# Z-score Alignment Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `scripts/audit_zscore_alignment.R`, which checks per trait that the pre-imputed bNMF z-score matrix is oriented to the T2D risk allele by comparing signs with GWAS Catalog (REST API v2) associations.

**Architecture:** One script of small, sourceable functions plus a command-line block. Catalog data are fetched per EFO term through an injectable `get_page(url)` function (so tests use canned responses), cached as RDS, aligned to the matrix by rsID, reduced to one call per SNP–trait pair, and summarised as a per-trait verdict. Output is a comparisons CSV, a per-trait CSV, and a self-contained HTML report built with the helpers in `scripts/qc_report.R`.

**Tech Stack:** R ≥ 4.1 (native pipe), `httr2` (HTTP), base R, `testthat` (tests, already used by `tests/run_tests.R`).

Spec: `docs/superpowers/specs/2026-11-13-zscore-alignment-audit-design.md`

> **Executed with a deviation:** Tasks 1–6 were built as written. In Task 7's
> live run the REST API took about 45 s per page, so the loader was replaced
> with the bulk associations file (see the revision note in the spec).
> Tasks 1–2's API code, `httr2` and the cache no longer exist; see commit
> 799872b for the replacement and its tests.

## Global Constraints

- Do not use `gwasrapidd`: it calls the retired v1 API (HTTP 410). Use v2 at `https://www.ebi.ac.uk/gwas/rest/api/v2`.
- Never compare against the user's own summary-statistics files; the Catalog is the only reference.
- Dependencies: `httr2` and base R only in the script (no `svglite`, `dplyr`, `data.table`, `jsonlite` calls). The HTML uses `qc_html_escape()`, `qc_html_table()`, `qc_name_list()` from `scripts/qc_report.R`.
- v2 facts (verified 2026-11-13): `associations?efo_id=X&size=N&page=P` with N ≤ 500 (1000 returns HTTP 400); response has `page$totalPages` (0 for an unknown term, with no `_embedded`) and `_embedded$associations`; each association has `association_id`, `accession_id`, `p_value`, `beta` (text such as `"0.018 unit decrease"` or `"-"`), `or_per_copy_num` (absent when null), `snp_allele` (list of `rs_id`, `effect_allele`), `reported_trait` (list), `multi_snp_haplotype`, `snp_interaction`. `efo-traits?trait=Q&size=N` returns `_embedded$efo_traits` with `efo_id`, `efo_trait`.
- Thresholds: genome-wide P ≤ 5e-8; informative when |matrix z| ≥ 2 (CLI `--min-abs-z`); verdicts need n ≥ 5 informative; PASS c ≥ 0.8; FLIPPED c ≤ 0.2; REVIEW otherwise.
- Rscript encodes spaces in `--file=` as `~+~`; the CLI must decode them (the user's repository path contains spaces).
- HTTP user agent is the generic string `bnmf_clustering z-score audit`; no personal details in requests.
- Git: work on branch `zscore-alignment-audit`. Stage explicit paths only — never `git add -A` or `git add .` (the user keeps private files in `scripts/`).
- Tests run with `Rscript tests/run_tests.R` from the repository root and must need no network access.

---

### Task 1: Catalog record parsing

**Files:**
- Create: `scripts/audit_zscore_alignment.R`
- Create: `tests/test-zscore-audit.R`

**Interfaces:**
- Produces: constants `CATALOG_API`, `AUDIT_GWS_P`, `AUDIT_BASES`; `audit_complement(x)`; `catalog_direction(beta_text, or_value) -> numeric (+1, -1, NA)`; `parse_catalog_associations(associations) -> data.frame(association_id chr, accession chr, rsID chr, effect_allele chr, p_value num, beta_text chr, or_value num, reported_trait chr, single_snp lgl)`. Test helper `fake_assoc(...)` returns one parsed-JSON association list.

- [ ] **Step 1: Create the branch**

```bash
git checkout -b zscore-alignment-audit
```

- [ ] **Step 2: Write the failing tests**

Create `tests/test-zscore-audit.R`:

```r
# Z-score alignment audit: scripts/audit_zscore_alignment.R
# No network access: Catalog responses are canned lists shaped like the
# parsed JSON of REST API v2.

suppressWarnings(suppressPackageStartupMessages({
  source(file.path(scripts_dir, "qc_report.R"))
  source(file.path(scripts_dir, "audit_zscore_alignment.R"))
}))

fake_assoc <- function(id, rs, ea, p = 1e-10, beta = "0.1 unit increase",
                       or = NULL, accession = "GCST000001", reported = "Trait X",
                       haplotype = FALSE, n_alleles = 1L) {
  list(association_id = id, accession_id = accession,
       snp_allele = rep(list(list(rs_id = rs, effect_allele = ea)), n_alleles),
       p_value = p, beta = beta, or_per_copy_num = or,
       reported_trait = list(reported), multi_snp_haplotype = haplotype,
       snp_interaction = FALSE)
}

test_that("catalog_direction reads beta text and odds ratios", {
  expect_equal(
    catalog_direction(c("0.1 unit increase", "0.2 unit decrease", "-", NA, "-", "-"),
                      c(NA, NA, 1.2, NA, 0.8, 1)),
    c(1, -1, 1, NA, -1, NA)
  )
  expect_equal(catalog_direction("0.1 unit decrease", 1.3), 1)
})

test_that("parse_catalog_associations flattens v2 records", {
  d <- parse_catalog_associations(list(
    fake_assoc(1, "rs1", "a", p = 2e-9, reported = "BMI"),
    fake_assoc(228396747, "rs2", "G", beta = "-", or = 1.1),
    fake_assoc(3, "rs3", "T", haplotype = TRUE),
    fake_assoc(4, "rs4", "C", n_alleles = 2L)
  ))
  expect_equal(d$association_id, c("1", "228396747", "3", "4"))
  expect_equal(d$rsID, c("rs1", "rs2", "rs3", "rs4"))
  expect_equal(d$effect_allele[1], "A")
  expect_equal(d$p_value[1], 2e-9)
  expect_equal(d$or_value, c(NA, 1.1, NA, NA))
  expect_equal(d$reported_trait[1], "BMI")
  expect_equal(d$single_snp, c(TRUE, TRUE, FALSE, FALSE))
})

test_that("parse_catalog_associations handles nulls and empty input", {
  d <- parse_catalog_associations(list(fake_assoc(5, NULL, NULL, p = NULL, beta = NULL)))
  expect_true(is.na(d$rsID))
  expect_true(is.na(d$effect_allele))
  expect_true(is.na(d$p_value))
  empty <- parse_catalog_associations(list())
  expect_equal(nrow(empty), 0L)
  expect_named(empty, names(d))
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `Rscript tests/run_tests.R`
Expected: `test-zscore-audit.R` errors with `cannot open file '.../scripts/audit_zscore_alignment.R'`.

- [ ] **Step 4: Write the implementation**

Create `scripts/audit_zscore_alignment.R`:

```r
# =============================================================================
# audit_zscore_alignment.R
# =============================================================================
# Checks that each trait column of the pre-imputed z-score matrix
# (zmat_preImputed.rds) is oriented to the T2D risk allele, by comparing signs
# with genome-wide-significant associations in the NHGRI-EBI GWAS Catalog.
#
# The check is independent of your own summary-statistics files, so it can
# catch an allele mistake made while preprocessing them (for example, the
# wrong effect-allele column), which flips or scrambles a whole trait.
#
# Usage:
#   1. Draft the audit manifest (one EFO term suggested per matrix trait):
#        Rscript scripts/audit_zscore_alignment.R --draft \
#          results/zmat_preImputed.rds audit_traits.tsv
#   2. Review audit_traits.tsv: correct efo_id (see the suggestions column;
#      separate several IDs with ";"), blank it to skip a trait, and
#      optionally add the GCST accession of the study the clustering used
#      and a reported_trait regular expression.
#   3. Run the audit:
#        Rscript scripts/audit_zscore_alignment.R \
#          results/zmat_preImputed.rds results/alignment_GWAS_summStats.csv \
#          audit_traits.tsv results/zscore_audit [--min-abs-z 2]
#      Writes results/zscore_audit.csv (one row per SNP-trait comparison),
#      results/zscore_audit_traits.csv (one row per trait),
#      results/zscore_audit.html, and results/zscore_audit_catalog_cache/.
#
# Verdicts: PASS (>= 80% of informative SNPs agree in sign), FLIPPED (<= 20%),
# REVIEW (in between), INSUFFICIENT (< 5 informative SNPs), LOOKUP_ERROR,
# NOT_CHECKED (no EFO ID).
#
# From R, source scripts/qc_report.R before this file.
#
# Dependencies: httr2.
# =============================================================================

CATALOG_API <- "https://www.ebi.ac.uk/gwas/rest/api/v2"
AUDIT_GWS_P <- 5e-8
AUDIT_BASES <- c("A", "C", "G", "T")

audit_complement <- function(x) chartr("ACGT", "TGCA", x)

catalog_direction <- function(beta_text, or_value) {
  # +1 / -1 for the effect allele, NA when the record gives no direction.
  # An odds ratio, when present, takes precedence over the beta text.
  beta_text <- tolower(ifelse(is.na(beta_text), "", beta_text))
  direction <- ifelse(grepl("increase", beta_text, fixed = TRUE), 1,
                      ifelse(grepl("decrease", beta_text, fixed = TRUE), -1, NA_real_))
  use_or <- is.finite(or_value) & or_value > 0 & or_value != 1
  direction[use_or] <- sign(log(or_value[use_or]))
  direction
}

catalog_scalar <- function(x, type) {
  if (is.null(x) || !length(x)) {
    return(if (type == "character") NA_character_ else NA_real_)
  }
  x <- x[[1]]
  if (type == "numeric") return(suppressWarnings(as.numeric(x)))
  if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE) else as.character(x)
}

parse_catalog_associations <- function(associations) {
  # Flatten v2 association records (parsed JSON) into one row each.
  chr <- function(f) vapply(associations, function(a) catalog_scalar(f(a), "character"),
                            character(1))
  num <- function(f) vapply(associations, function(a) catalog_scalar(f(a), "numeric"),
                            numeric(1))
  first_allele <- function(a) if (length(a$snp_allele)) a$snp_allele[[1]] else list()
  data.frame(
    association_id = chr(function(a) a$association_id),
    accession = chr(function(a) a$accession_id),
    rsID = chr(function(a) first_allele(a)$rs_id),
    effect_allele = toupper(chr(function(a) first_allele(a)$effect_allele)),
    p_value = num(function(a) a$p_value),
    beta_text = chr(function(a) a$beta),
    or_value = num(function(a) a$or_per_copy_num),
    reported_trait = vapply(associations, function(a) {
      paste(unlist(a$reported_trait), collapse = "; ")
    }, character(1)),
    single_snp = vapply(associations, function(a) {
      length(a$snp_allele) == 1L && !isTRUE(a$multi_snp_haplotype) &&
        !isTRUE(a$snp_interaction)
    }, logical(1)),
    stringsAsFactors = FALSE
  )
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript tests/run_tests.R`
Expected: all blocks pass, including the 3 new ones; last line `... 0 failed or errored.`

- [ ] **Step 6: Commit**

```bash
git add scripts/audit_zscore_alignment.R tests/test-zscore-audit.R
git commit -m "Add GWAS Catalog record parsing for the z-score alignment audit"
```

---

### Task 2: Catalog fetching, caching and EFO search

**Files:**
- Modify: `scripts/audit_zscore_alignment.R` (append)
- Modify: `tests/test-zscore-audit.R` (append)

**Interfaces:**
- Consumes: `parse_catalog_associations()`, `CATALOG_API`.
- Produces: `catalog_get_page(url) -> list` (parsed JSON); `catalog_associations_url(efo_id, page, page_size = 500L) -> chr`; `catalog_fetch_efo(efo_id, cache_dir = NULL, get_page = catalog_get_page, page_size = 500L) -> data.frame` (columns of `parse_catalog_associations`, unique `association_id`); `catalog_search_efo(query, get_page = catalog_get_page, n = 3L) -> data.frame(efo_id, efo_trait)`. Test helpers `fake_api(pages)` (environment with `$get(url)` and `$calls`) and `assoc_page(assocs, total)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test-zscore-audit.R`:

```r
fake_api <- function(pages) {
  env <- new.env()
  env$calls <- character(0)
  env$get <- function(url) {
    env$calls <- c(env$calls, url)
    if (!url %in% names(pages)) stop("unexpected URL: ", url)
    pages[[url]]
  }
  env
}

assoc_page <- function(assocs, total) {
  list(`_embedded` = list(associations = assocs), page = list(totalPages = total))
}

test_that("catalog_fetch_efo reads every page and caches the result", {
  api <- fake_api(setNames(
    list(assoc_page(list(fake_assoc(1, "rs1", "A"), fake_assoc(2, "rs2", "G")), 2),
         assoc_page(list(fake_assoc(2, "rs2", "G"), fake_assoc(3, "rs3", "C")), 2)),
    c(catalog_associations_url("EFO_0000001", 0), catalog_associations_url("EFO_0000001", 1))
  ))
  cache <- file.path(tempfile(), "cache")
  d <- catalog_fetch_efo("EFO_0000001", cache_dir = cache, get_page = api$get)
  expect_equal(d$association_id, c("1", "2", "3"))
  expect_length(api$calls, 2L)
  again <- catalog_fetch_efo("EFO_0000001", cache_dir = cache, get_page = api$get)
  expect_identical(again, d)
  expect_length(api$calls, 2L)
})

test_that("catalog_fetch_efo returns no rows for an unknown term", {
  api <- fake_api(setNames(list(list(page = list(totalPages = 0))),
                           catalog_associations_url("EFO_0000002", 0)))
  expect_equal(nrow(catalog_fetch_efo("EFO_0000002", get_page = api$get)), 0L)
})

test_that("catalog_fetch_efo rejects malformed IDs and responses", {
  expect_error(catalog_fetch_efo("EFO 1&x=2", get_page = function(url) stop("no call")),
               "Invalid EFO ID")
  api <- fake_api(setNames(list(list(oops = TRUE)), catalog_associations_url("EFO_0000003", 0)))
  expect_error(catalog_fetch_efo("EFO_0000003", get_page = api$get),
               "Unexpected GWAS Catalog response")
})

test_that("catalog_search_efo ranks exact and measurement labels first", {
  hits <- list(`_embedded` = list(efo_traits = list(
    list(efo_id = "HP_1", efo_trait = "Decreased HDL concentration"),
    list(efo_id = "EFO_2", efo_trait = "high density lipoprotein cholesterol measurement"),
    list(efo_id = "EFO_3", efo_trait = "HDL")
  )))
  d <- catalog_search_efo("HDL", get_page = function(url) {
    expect_match(url, "efo-traits?trait=HDL&", fixed = TRUE)
    hits
  })
  expect_equal(d$efo_id, c("EFO_3", "EFO_2", "HP_1"))
  none <- catalog_search_efo("zzz", get_page = function(url) list(page = list(totalPages = 0)))
  expect_equal(nrow(none), 0L)
  expect_named(none, c("efo_id", "efo_trait"))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript tests/run_tests.R`
Expected: the 4 new blocks fail with `could not find function "catalog_associations_url"` / `"catalog_search_efo"`.

- [ ] **Step 3: Write the implementation**

Append to `scripts/audit_zscore_alignment.R`:

```r
catalog_get_page <- function(url) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Install the httr2 package to query the GWAS Catalog.")
  }
  httr2::request(url) |>
    httr2::req_user_agent("bnmf_clustering z-score audit") |>
    httr2::req_timeout(60) |>
    httr2::req_retry(max_tries = 4) |>
    httr2::req_perform() |>
    httr2::resp_body_json()
}

catalog_associations_url <- function(efo_id, page, page_size = 500L) {
  sprintf("%s/associations?efo_id=%s&size=%d&page=%d",
          CATALOG_API, efo_id, as.integer(page_size), as.integer(page))
}

catalog_fetch_efo <- function(efo_id, cache_dir = NULL, get_page = catalog_get_page,
                              page_size = 500L) {
  # Every Catalog association annotated with one EFO term, cached as RDS.
  if (!grepl("^[A-Za-z]+_[0-9]+$", efo_id)) stop("Invalid EFO ID: ", efo_id)
  cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(efo_id, ".rds"))
  if (!is.null(cache_file) && file.exists(cache_file)) return(readRDS(cache_file))

  pages <- list()
  page <- 0L
  repeat {
    response <- get_page(catalog_associations_url(efo_id, page, page_size))
    total <- response$page$totalPages
    if (is.null(total)) stop("Unexpected GWAS Catalog response for ", efo_id)
    pages[[length(pages) + 1L]] <-
      parse_catalog_associations(response[["_embedded"]]$associations)
    page <- page + 1L
    if (total > 1L) message(sprintf("  %s: page %d of %d", efo_id, page, total))
    if (page >= total) break
  }
  out <- do.call(rbind, pages)
  out <- out[!duplicated(out$association_id), , drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(cache_file)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_file)
  }
  out
}

catalog_search_efo <- function(query, get_page = catalog_get_page, n = 3L) {
  # Up to n EFO terms whose label matches query, best first: exact label,
  # "<query> measurement", labels starting with the query, other
  # "measurement" terms, then the rest; shorter labels first within a tier.
  url <- sprintf("%s/efo-traits?trait=%s&size=50", CATALOG_API,
                 utils::URLencode(query, reserved = TRUE))
  hits <- get_page(url)[["_embedded"]]$efo_traits
  if (!length(hits)) {
    return(data.frame(efo_id = character(0), efo_trait = character(0)))
  }
  id <- vapply(hits, function(h) catalog_scalar(h$efo_id, "character"), character(1))
  label <- vapply(hits, function(h) catalog_scalar(h$efo_trait, "character"), character(1))
  q <- tolower(query)
  l <- tolower(label)
  tier <- ifelse(l == q, 0L,
          ifelse(l == paste(q, "measurement"), 1L,
          ifelse(startsWith(l, q), 2L,
          ifelse(grepl("measurement", l, fixed = TRUE), 3L, 4L))))
  ord <- order(tier, nchar(label))
  utils::head(data.frame(efo_id = id[ord], efo_trait = label[ord]), n)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript tests/run_tests.R`
Expected: `0 failed or errored.`

- [ ] **Step 5: Commit**

```bash
git add scripts/audit_zscore_alignment.R tests/test-zscore-audit.R
git commit -m "Fetch and cache GWAS Catalog v2 associations by EFO term"
```

---

### Task 3: Alignment to the matrix and per-SNP calls

**Files:**
- Modify: `scripts/audit_zscore_alignment.R` (append)
- Modify: `tests/test-zscore-audit.R` (append)

**Interfaces:**
- Consumes: `parse_catalog_associations()` output, `catalog_direction()`, `audit_complement()`, `AUDIT_BASES`, `AUDIT_GWS_P`.
- Produces:
  - `prepare_alignment(alignment, z) -> data.frame(ChrPos, rsID, REF, ALT, Risk_Allele, palindromic)` in `rownames(z)` order; invalid rsIDs set to `NA`.
  - `align_catalog_to_matrix(catalog, alignment, accession = "", reported_trait = "", p_threshold = AUDIT_GWS_P) -> data.frame` with the alignment and catalog columns plus `catalog_direction`, `strand_complemented`, `aligned_allele`, `expected_sign`, `same_study`; rows in matrix order.
  - `call_snp_trait(aligned, z_values, min_abs_z = 2) -> data.frame | NULL`, one row per `ChrPos`, adding `n_catalog_records`, `catalog_conflict`, `matrix_z`, `status` (`concordant`, `discordant`, `weak_matrix_z`, `missing_in_matrix`, `catalog_conflict`), `informative`, `concordant`. `z_values` is a numeric vector named by `ChrPos`.
  - Test helper `audit_fixture()` returning `list(alignment, records, catalog, z)`; `z` has columns `good`, `flipped`, `scrambled`, `lookup`, `blank`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test-zscore-audit.R`:

```r
audit_fixture <- function() {
  # Eight matrix variants. rs7 is A/T (palindromic). Every Catalog record
  # below implies +1 for the T2D risk allele, reached in different ways.
  alignment <- data.frame(
    ChrPos = paste0("1:", 1:8 * 100),
    rsID = paste0("rs", 1:8),
    REF = c("A", "C", "A", "G", "A", "C", "A", "G"),
    ALT = c("G", "T", "C", "T", "G", "A", "T", "A"),
    Risk_Allele = c("G", "C", "A", "T", "A", "A", "T", "G"),
    stringsAsFactors = FALSE
  )
  records <- list(
    fake_assoc(1, "rs1", "G"),                                # risk allele, increase
    fake_assoc(2, "rs2", "T", beta = "0.2 unit decrease"),    # non-risk, decrease
    fake_assoc(3, "rs3", "T"),                                # complement of risk A
    fake_assoc(4, "rs4", "T", beta = "-", or = 1.2),          # risk, OR > 1
    fake_assoc(5, "rs5", "G", beta = "-", or = 0.9),          # non-risk, OR < 1
    fake_assoc(6, "rs6", "A"),
    fake_assoc(7, "rs7", "T"),                                # palindromic: dropped
    fake_assoc(8, "rs8", "G")
  )
  z <- cbind(good = rep(5, 8), flipped = rep(-5, 8),
             scrambled = c(5, -5, 5, -5, 5, -5, 5, -5),
             lookup = rep(5, 8), blank = rep(5, 8))
  rownames(z) <- alignment$ChrPos
  list(alignment = alignment, records = records,
       catalog = parse_catalog_associations(records), z = z)
}

test_that("prepare_alignment orders rows like the matrix and flags palindromes", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment[8:1, ], fx$z)
  expect_equal(a$ChrPos, rownames(fx$z))
  expect_equal(a$palindromic, c(rep(FALSE, 6), TRUE, FALSE))
})

test_that("prepare_alignment rejects missing rows, bad risk alleles and absent rsIDs", {
  fx <- audit_fixture()
  expect_error(prepare_alignment(fx$alignment[-1, ], fx$z), "no alignment row")
  bad <- fx$alignment
  bad$Risk_Allele[1] <- "C"
  expect_error(prepare_alignment(bad, fx$z), "Risk_Allele")
  no_rsid <- fx$alignment
  no_rsid$rsID <- ""
  expect_error(prepare_alignment(no_rsid, fx$z), "rsID")
})

test_that("align_catalog_to_matrix orients every record to the T2D risk allele", {
  fx <- audit_fixture()
  d <- align_catalog_to_matrix(fx$catalog, prepare_alignment(fx$alignment, fx$z))
  expect_equal(d$rsID, c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs8"))
  expect_true(all(d$expected_sign == 1))
  expect_equal(d$strand_complemented, c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("align_catalog_to_matrix applies the P, allele and reported-trait filters", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment, fx$z)
  catalog <- parse_catalog_associations(list(
    fake_assoc(1, "rs1", "G", p = 1e-6),                    # not genome-wide
    fake_assoc(2, "rs2", "?"),                               # unknown allele
    fake_assoc(3, "rs3", "AT"),                              # not a single base
    fake_assoc(4, "rs4", "T", beta = "-"),                   # no direction
    fake_assoc(5, "rs5", "G", reported = "BMI-adjusted thing"),
    fake_assoc(6, "rs6", "A", reported = "Unadjusted thing"),
    fake_assoc(7, "rs99", "A")                               # not in the matrix
  ))
  d <- align_catalog_to_matrix(catalog, a, reported_trait = "unadjusted")
  expect_equal(d$rsID, "rs6")
  expect_equal(nrow(align_catalog_to_matrix(catalog[0, ], a)), 0L)
})

test_that("call_snp_trait prefers the clustering study and flags conflicts", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment, fx$z)
  catalog <- parse_catalog_associations(list(
    fake_assoc(1, "rs1", "G", p = 1e-30, accession = "GCST000008"),
    fake_assoc(2, "rs1", "G", p = 1e-9, accession = "GCST000009"),
    fake_assoc(3, "rs2", "C", p = 1e-12),                    # risk allele: +1
    fake_assoc(4, "rs2", "T", p = 1e-10),                    # non-risk: -1
    fake_assoc(5, "rs3", "A", p = 1e-10)
  ))
  d <- align_catalog_to_matrix(catalog, a, accession = "GCST000009")
  z <- setNames(c(5, 5, 1, rep(5, 5)), rownames(fx$z))
  calls <- call_snp_trait(d, z, min_abs_z = 2)
  expect_equal(calls$rsID, c("rs1", "rs2", "rs3"))
  expect_equal(calls$association_id, c("2", "3", "5"))
  expect_equal(calls$same_study, c(TRUE, FALSE, FALSE))
  expect_equal(calls$n_catalog_records, c(2L, 2L, 1L))
  expect_equal(calls$status, c("concordant", "catalog_conflict", "weak_matrix_z"))
  expect_equal(calls$informative, c(TRUE, FALSE, FALSE))
  z[3] <- NA
  expect_equal(call_snp_trait(d, z)$status[3], "missing_in_matrix")
  expect_null(call_snp_trait(d[0, ], z))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript tests/run_tests.R`
Expected: the 5 new blocks fail with `could not find function "prepare_alignment"`.

- [ ] **Step 3: Write the implementation**

Append to `scripts/audit_zscore_alignment.R`:

```r
prepare_alignment <- function(alignment, z) {
  # One alignment row per matrix row, in matrix order.
  needed <- c("ChrPos", "rsID", "REF", "ALT", "Risk_Allele")
  missing_cols <- setdiff(needed, names(alignment))
  if (length(missing_cols)) {
    stop("The alignment file lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  idx <- match(rownames(z), alignment$ChrPos)
  if (anyNA(idx)) {
    stop(sum(is.na(idx)), " matrix row(s) have no alignment row (first: ",
         rownames(z)[is.na(idx)][1], ").")
  }
  a <- as.data.frame(alignment, stringsAsFactors = FALSE)[idx, needed]
  a[] <- lapply(a, as.character)
  a[c("REF", "ALT", "Risk_Allele")] <- lapply(a[c("REF", "ALT", "Risk_Allele")], toupper)
  bad_risk <- is.na(a$Risk_Allele) | !(a$Risk_Allele == a$REF | a$Risk_Allele == a$ALT)
  if (any(bad_risk)) {
    stop(sum(bad_risk), " variant(s) have a Risk_Allele that is neither REF nor ALT (first: ",
         a$ChrPos[bad_risk][1], ").")
  }
  a$rsID[is.na(a$rsID) | !grepl("^rs[0-9]+$", a$rsID)] <- NA_character_
  if (all(is.na(a$rsID))) {
    stop("The alignment file has no rsIDs; the audit matches the Catalog by rsID. ",
         "rsIDs come from LDlink, so they are blank when LD pruning is off.")
  }
  a$palindromic <- paste0(a$REF, a$ALT) %in% c("AT", "TA", "CG", "GC")
  rownames(a) <- NULL
  a
}

align_catalog_to_matrix <- function(catalog, alignment, accession = "",
                                    reported_trait = "", p_threshold = AUDIT_GWS_P) {
  # Keep usable genome-wide-significant records, match them to matrix rows by
  # rsID, and express each record's direction for the T2D risk allele.
  catalog$catalog_direction <- catalog_direction(catalog$beta_text, catalog$or_value)
  keep <- catalog$single_snp & !is.na(catalog$rsID) & is.finite(catalog$p_value) &
    catalog$p_value <= p_threshold & !is.na(catalog$catalog_direction)
  if (nzchar(reported_trait)) {
    keep <- keep & grepl(reported_trait, catalog$reported_trait, ignore.case = TRUE,
                         perl = TRUE)
  }
  usable <- alignment[!is.na(alignment$rsID) & !alignment$palindromic, , drop = FALSE]
  d <- merge(usable, catalog[keep, , drop = FALSE], by = "rsID")

  ea <- d$effect_allele
  base <- !is.na(ea) & ea %in% AUDIT_BASES
  direct <- base & (ea == d$REF | ea == d$ALT)
  flipped_strand <- base & !direct &
    (audit_complement(ea) == d$REF | audit_complement(ea) == d$ALT)
  matched <- direct | flipped_strand
  d <- d[matched, , drop = FALSE]
  d$strand_complemented <- flipped_strand[matched]
  d$aligned_allele <- ifelse(d$strand_complemented, audit_complement(d$effect_allele),
                             d$effect_allele)
  d$expected_sign <- d$catalog_direction * ifelse(d$aligned_allele == d$Risk_Allele, 1, -1)
  d$same_study <- nzchar(accession) & !is.na(d$accession) & d$accession == accession
  d <- d[order(match(d$ChrPos, alignment$ChrPos), d$p_value), , drop = FALSE]
  rownames(d) <- NULL
  d
}

call_snp_trait <- function(aligned, z_values, min_abs_z = 2) {
  # One call per matrix variant: the smallest-P record from the clustering
  # study when there is one, otherwise from any study. Records in that pool
  # that disagree in direction make the call a catalog_conflict.
  if (!nrow(aligned)) return(NULL)
  rows <- lapply(split(seq_len(nrow(aligned)), aligned$ChrPos), function(i) {
    g <- aligned[i, , drop = FALSE]
    pool <- if (any(g$same_study)) g[g$same_study, , drop = FALSE] else g
    best <- pool[which.min(pool$p_value), , drop = FALSE]
    best$n_catalog_records <- nrow(g)
    best$catalog_conflict <- length(unique(pool$expected_sign)) > 1L
    best
  })
  calls <- do.call(rbind, rows)
  calls <- calls[order(match(calls$ChrPos, names(z_values))), , drop = FALSE]
  calls$matrix_z <- unname(z_values[calls$ChrPos])
  calls$status <- ifelse(calls$catalog_conflict, "catalog_conflict",
                  ifelse(!is.finite(calls$matrix_z), "missing_in_matrix",
                  ifelse(abs(calls$matrix_z) < min_abs_z, "weak_matrix_z",
                  ifelse(sign(calls$matrix_z) == calls$expected_sign,
                         "concordant", "discordant"))))
  calls$informative <- calls$status %in% c("concordant", "discordant")
  calls$concordant <- ifelse(calls$informative, calls$status == "concordant", NA)
  rownames(calls) <- NULL
  calls
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript tests/run_tests.R`
Expected: `0 failed or errored.`

- [ ] **Step 5: Commit**

```bash
git add scripts/audit_zscore_alignment.R tests/test-zscore-audit.R
git commit -m "Align Catalog records to the T2D risk allele and call each SNP"
```

---

### Task 4: Per-trait verdicts

**Files:**
- Modify: `scripts/audit_zscore_alignment.R` (append)
- Modify: `tests/test-zscore-audit.R` (append)

**Interfaces:**
- Consumes: `call_snp_trait()` output (uses only `informative`, `concordant`, `same_study`).
- Produces: `AUDIT_VERDICT_ORDER`; `trait_verdict(calls, accession = "", min_n = 5L, pass_at = 0.8, flip_at = 0.2) -> one-row data.frame(verdict, n_catalog_snps, n_informative, n_concordant, concordance, binom_p, n_same_study, same_study_concordance, study_differs)`. `calls = NULL` gives `INSUFFICIENT` with zero counts; `study_differs` is `NA` when `accession` is blank.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test-zscore-audit.R`:

```r
test_that("trait_verdict separates aligned, flipped and scrambled traits", {
  fx <- audit_fixture()
  d <- align_catalog_to_matrix(fx$catalog, prepare_alignment(fx$alignment, fx$z))
  verdict <- function(trait) trait_verdict(call_snp_trait(d, fx$z[, trait]))
  good <- verdict("good")
  flipped <- verdict("flipped")
  scrambled <- verdict("scrambled")
  expect_equal(c(good$verdict, flipped$verdict, scrambled$verdict),
               c("PASS", "FLIPPED", "REVIEW"))
  expect_equal(good$n_informative, 7L)
  expect_equal(good$concordance, 1)
  expect_equal(flipped$concordance, 0)
  expect_equal(scrambled$n_concordant, 3L)
  expect_lt(flipped$binom_p, 0.05)
  expect_true(is.na(good$study_differs))
})

test_that("trait_verdict needs enough informative SNPs and reports same-study agreement", {
  expect_equal(trait_verdict(NULL)$verdict, "INSUFFICIENT")
  expect_equal(trait_verdict(NULL)$n_informative, 0L)
  calls <- data.frame(informative = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
                      concordant = c(TRUE, TRUE, TRUE, TRUE, FALSE, NA),
                      same_study = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
  v <- trait_verdict(calls, accession = "GCST000009")
  expect_equal(v$verdict, "PASS")
  expect_equal(v$concordance, 0.8)
  expect_equal(v$n_same_study, 2L)
  expect_equal(v$same_study_concordance, 1)
  expect_false(v$study_differs)
  expect_equal(trait_verdict(calls[1:4, ])$verdict, "INSUFFICIENT")
  calls$same_study <- FALSE
  expect_true(trait_verdict(calls, accession = "GCST000009")$study_differs)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript tests/run_tests.R`
Expected: the 2 new blocks fail with `could not find function "trait_verdict"`.

- [ ] **Step 3: Write the implementation**

Append to `scripts/audit_zscore_alignment.R`:

```r
AUDIT_VERDICT_ORDER <- c("FLIPPED", "REVIEW", "LOOKUP_ERROR", "INSUFFICIENT",
                         "PASS", "NOT_CHECKED")

trait_verdict <- function(calls, accession = "", min_n = 5L, pass_at = 0.8,
                          flip_at = 0.2) {
  informative <- if (is.null(calls)) logical(0) else calls$informative
  n <- sum(informative)
  k <- if (n) sum(calls$concordant[informative]) else 0L
  concordance <- if (n) k / n else NA_real_
  same <- if (n) informative & calls$same_study else logical(0)
  n_same <- sum(same)
  verdict <- if (n < min_n) {
    "INSUFFICIENT"
  } else if (concordance >= pass_at) {
    "PASS"
  } else if (concordance <= flip_at) {
    "FLIPPED"
  } else {
    "REVIEW"
  }
  data.frame(
    verdict = verdict,
    n_catalog_snps = if (is.null(calls)) 0L else nrow(calls),
    n_informative = n,
    n_concordant = k,
    concordance = concordance,
    binom_p = if (n) stats::binom.test(k, n)$p.value else NA_real_,
    n_same_study = n_same,
    same_study_concordance = if (n_same) mean(calls$concordant[same]) else NA_real_,
    study_differs = if (nzchar(accession)) n_same == 0L else NA,
    stringsAsFactors = FALSE
  )
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript tests/run_tests.R`
Expected: `0 failed or errored.`

- [ ] **Step 5: Commit**

```bash
git add scripts/audit_zscore_alignment.R tests/test-zscore-audit.R
git commit -m "Add per-trait sign-concordance verdicts"
```

---

### Task 5: HTML report

**Files:**
- Modify: `scripts/audit_zscore_alignment.R` (append)
- Modify: `tests/test-zscore-audit.R` (append)

**Interfaces:**
- Consumes: `qc_html_escape()`, `qc_html_table()`, `qc_name_list()` (from `scripts/qc_report.R`); `AUDIT_VERDICT_ORDER`; summary rows = `trait, efo_id, accession, reported_trait` + `trait_verdict()` columns + `note`; comparisons = `trait` + `call_snp_trait()` columns.
- Produces: `concordance_plot_svg(summary) -> chr` (one string); `write_audit_report(summary, comparisons, file, settings = list()) -> invisible(flags chr)`. Recognised settings: `min_abs_z`, `p_threshold`, `n_variants`, `n_traits`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test-zscore-audit.R`:

```r
test_that("concordance_plot_svg escapes trait names and skips traits without data", {
  s <- data.frame(trait = c("A<B", "none"), verdict = c("PASS", "INSUFFICIENT"),
                  n_informative = c(10L, 0L), n_concordant = c(9L, 0L),
                  concordance = c(0.9, NA))
  svg <- concordance_plot_svg(s)
  expect_match(svg, "A&lt;B", fixed = TRUE)
  expect_false(grepl(">none<", svg, fixed = TRUE))
  expect_match(concordance_plot_svg(s[2, ]), "No trait", fixed = TRUE)
})

test_that("write_audit_report writes a self-contained page with every trait", {
  fx <- audit_fixture()
  d <- align_catalog_to_matrix(fx$catalog, prepare_alignment(fx$alignment, fx$z))
  traits <- c("good", "flipped", "scrambled")
  calls <- lapply(traits, function(t) cbind(trait = t, call_snp_trait(d, fx$z[, t])))
  summary <- do.call(rbind, lapply(seq_along(traits), function(i) {
    cbind(trait = traits[i], efo_id = "EFO_0000001", accession = "", reported_trait = "",
          trait_verdict(calls[[i]]), note = "")
  }))
  blank <- cbind(trait = "blank", efo_id = "", accession = "", reported_trait = "",
                 trait_verdict(NULL), note = "No EFO ID given.")
  blank$verdict <- "NOT_CHECKED"
  summary <- rbind(summary, blank)
  f <- file.path(tempfile(), "report.html")
  flags <- quietly(write_audit_report(summary, do.call(rbind, calls), f,
                                      settings = list(min_abs_z = 2)))
  html <- paste(readLines(f), collapse = "\n")
  for (t in c(traits, "blank")) expect_match(html, paste0(">", t, "<"), fixed = TRUE)
  expect_match(html, "<svg", fixed = TRUE)
  expect_match(html, "<details open>", fixed = TRUE)
  expect_false(grepl("<script", html, fixed = TRUE))
  expect_true(any(grepl("flipped", flags, fixed = TRUE)))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript tests/run_tests.R`
Expected: the 2 new blocks fail with `could not find function "concordance_plot_svg"` / `"write_audit_report"`.

- [ ] **Step 3: Write the implementation**

Append to `scripts/audit_zscore_alignment.R`:

```r
AUDIT_COLORS <- c(PASS = "#2f7d4f", REVIEW = "#b26b00", FLIPPED = "#b3261e")

audit_setting <- function(settings, name, default) {
  if (is.null(settings[[name]])) default else settings[[name]]
}

concordance_plot_svg <- function(summary) {
  # Dot plot: one row per trait with informative SNPs, x = share concordant,
  # dot area grows with the number of informative SNPs.
  s <- summary[summary$n_informative > 0 & is.finite(summary$concordance), , drop = FALSE]
  if (!nrow(s)) return("<p class=\"muted\">No trait had informative SNPs to plot.</p>")
  s <- s[order(s$concordance, s$trait), , drop = FALSE]
  left <- 240; plot_w <- 480; right <- 60; top <- 16; row_h <- 18
  width <- left + plot_w + right
  bottom <- top + nrow(s) * row_h
  height <- bottom + 40
  f <- function(v) formatC(v, format = "f", digits = 1)
  x <- function(v) left + v * plot_w
  y <- top + (seq_len(nrow(s)) - 0.5) * row_h
  r <- pmin(8, 2.5 + sqrt(s$n_informative) / 1.5)
  col <- unname(AUDIT_COLORS[s$verdict])
  col[is.na(col)] <- "#7a8794"
  ticks <- c(0, 0.2, 0.5, 0.8, 1)
  paste(c(
    sprintf(paste0('<svg viewBox="0 0 %d %s" width="100%%" role="img" ',
                   'aria-label="Sign concordance by trait" ',
                   'style="max-width:%dpx;font-size:11px;font-family:inherit">'),
            width, f(height), width),
    sprintf('<rect x="%s" y="%d" width="%s" height="%s" fill="#b3261e" fill-opacity="0.08"/>',
            f(x(0)), top, f(0.2 * plot_w), f(bottom - top)),
    sprintf('<rect x="%s" y="%d" width="%s" height="%s" fill="#b26b00" fill-opacity="0.07"/>',
            f(x(0.2)), top, f(0.6 * plot_w), f(bottom - top)),
    sprintf(paste0('<line x1="%1$s" x2="%1$s" y1="%2$d" y2="%3$s" stroke="#c9ced8"/>',
                   '<text x="%1$s" y="%4$s" text-anchor="middle" fill="#6a7486">%5$s</text>'),
            f(x(ticks)), top, f(bottom), f(bottom + 14), format(ticks)),
    sprintf(paste0('<text x="%s" y="%s" text-anchor="middle" fill="#3b4658">',
                   'Share of informative SNPs whose matrix sign matches the Catalog</text>'),
            f(x(0.5)), f(bottom + 32)),
    sprintf('<text x="%d" y="%s" text-anchor="end" dominant-baseline="middle" fill="#1c2330">%s</text>',
            left - 8, f(y), qc_html_escape(s$trait)),
    sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s"><title>%s: %d of %d concordant (%s)</title></circle>',
            f(x(s$concordance)), f(y), f(r), col, qc_html_escape(s$trait),
            as.integer(s$n_concordant), as.integer(s$n_informative), s$verdict),
    sprintf('<text x="%s" y="%s" dominant-baseline="middle" fill="#6a7486">n=%d</text>',
            f(x(1) + 12), f(y), as.integer(s$n_informative)),
    "</svg>"
  ), collapse = "\n")
}

write_audit_report <- function(summary, comparisons, file, settings = list()) {
  if (!exists("qc_html_table", mode = "function")) {
    stop("Source scripts/qc_report.R before writing the audit report.")
  }
  min_abs_z <- audit_setting(settings, "min_abs_z", 2)
  p_threshold <- audit_setting(settings, "p_threshold", AUDIT_GWS_P)
  summary <- summary[order(match(summary$verdict, AUDIT_VERDICT_ORDER), summary$trait), ,
                     drop = FALSE]
  count <- table(factor(summary$verdict, levels = AUDIT_VERDICT_ORDER))
  traits_with <- function(v) qc_name_list(summary$trait[summary$verdict == v])

  flags <- character(0)
  if (count[["FLIPPED"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) look flipped: the matrix sign disagrees with the Catalog at most informative SNPs (%s). Check which column was used as the effect allele when these GWAS were formatted.",
    count[["FLIPPED"]], traits_with("FLIPPED")))
  if (count[["REVIEW"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) have mixed agreement (%s). Possible causes: a partial allele mix-up, or an EFO term that pools a differently defined phenotype (for example adjusted and unadjusted). Their SNP tables are open below.",
    count[["REVIEW"]], traits_with("REVIEW")))
  if (count[["LOOKUP_ERROR"]] > 0) flags <- c(flags, sprintf(
    "The Catalog lookup failed for %d trait(s) (%s); rerun to retry.",
    count[["LOOKUP_ERROR"]], traits_with("LOOKUP_ERROR")))
  differs <- summary$trait[summary$verdict %in% c("PASS", "REVIEW", "FLIPPED") &
                             summary$study_differs %in% TRUE]
  if (length(differs)) flags <- c(flags, sprintf(
    "No informative Catalog record came from the clustering study for: %s. Their verdicts rest on other studies of the trait.",
    qc_name_list(differs)))

  check_block <- if (length(flags)) {
    paste0("<div class=\"checks\"><h2>Things to check</h2><ul>",
           paste0("<li>", qc_html_escape(flags), "</li>", collapse = ""), "</ul></div>")
  } else {
    "<div class=\"checks ok\"><h2>Things to check</h2><p>No flags were raised.</p></div>"
  }

  summary_cols <- c("trait", "verdict", "n_informative", "n_concordant", "concordance",
                    "binom_p", "n_catalog_snps", "n_same_study", "same_study_concordance",
                    "study_differs", "efo_id", "accession", "reported_trait", "note")
  detail_cols <- c("rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "effect_allele",
                   "strand_complemented", "catalog_direction", "expected_sign", "matrix_z",
                   "status", "p_value", "accession", "same_study", "association_id",
                   "reported_trait")
  details <- vapply(seq_len(nrow(summary)), function(i) {
    s <- summary[i, ]
    rows <- if (is.null(comparisons)) NULL else comparisons[comparisons$trait == s$trait, , drop = FALSE]
    if (is.null(rows) || !nrow(rows)) return("")
    sprintf("<details%s><summary><strong>%s</strong> &middot; %s (%d of %d informative SNPs concordant)</summary>%s</details>",
            if (s$verdict %in% c("FLIPPED", "REVIEW")) " open" else "",
            qc_html_escape(s$trait), s$verdict, as.integer(s$n_concordant),
            as.integer(s$n_informative), qc_html_table(rows[detail_cols], max_rows = 500L))
  }, character(1))

  meta <- sprintf("%s variants &times; %s traits &middot; generated %s",
                  audit_setting(settings, "n_variants", "?"),
                  audit_setting(settings, "n_traits", nrow(summary)),
                  format(Sys.time(), "%Y-%m-%d %H:%M"))
  counts_line <- paste(sprintf("%s: %d", names(count), as.integer(count)), collapse = " &middot; ")
  methods <- sprintf(paste(
    "<p>For each trait, every GWAS Catalog association annotated with the trait's EFO term(s) was taken from any study,",
    "keeping single-SNP records with P &le; %s and a stated direction (an odds ratio, or a beta described as an increase or decrease).",
    "Records were matched to matrix variants by rsID. A/T and C/G variants were excluded; an effect allele matching only the",
    "complement of REF/ALT was read on the opposite strand. Each record's direction was expressed for the T2D risk allele.</p>",
    "<p>Each variant gets one call: the smallest-P record from the clustering study when its accession is given and present,",
    "otherwise from any study. When those records disagree in direction the call is a Catalog conflict and is not counted.",
    "A call is informative when the matrix |z| &ge; %s. Verdicts need at least 5 informative variants: PASS when at least 80%%",
    "agree in sign, FLIPPED when at most 20%% do, REVIEW in between. The binomial P tests agreement against 50%%.</p>",
    "<p>Limitations: an EFO term can pool related phenotypes whose directions differ (adjusted versus unadjusted, sex-specific);",
    "the reported_trait filter narrows a term. Catalog records are curated top hits, so traits with few significant loci among",
    "these variants may be INSUFFICIENT. Only the sign is checked, not the magnitude.</p>"),
    format(p_threshold), format(min_abs_z))

  css <- paste(
    "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f6f7f9;color:#1c2330;}",
    "main{max-width:1100px;margin:0 auto;padding:24px 16px 64px;}",
    "h1{font-size:24px;margin:8px 0 4px;} h2{font-size:18px;margin:32px 0 8px;}",
    "p{line-height:1.5;margin:6px 0;} .muted{color:#6a7486;font-size:13px;}",
    ".meta{color:#6a7486;font-size:13px;margin-bottom:16px;}",
    ".checks{background:#fff7e6;border:1px solid #f0c36d;border-radius:8px;padding:4px 16px 12px;}",
    ".checks.ok{background:#eef8f0;border-color:#9fd3aa;} .checks h2{margin-top:12px;} .checks li{margin:6px 0;line-height:1.45;}",
    ".tablewrap{overflow-x:auto;background:#fff;border:1px solid #dfe3ea;border-radius:8px;}",
    "table{border-collapse:collapse;width:100%;font-size:13px;}",
    "th,td{padding:6px 10px;border-bottom:1px solid #eceff4;text-align:left;white-space:nowrap;}",
    "th{background:#f0f2f6;font-weight:600;} tr:last-child td{border-bottom:none;}",
    ".num{text-align:right;font-variant-numeric:tabular-nums;}",
    "details{margin:8px 0;} summary{cursor:pointer;padding:4px 0;}",
    ".plot{background:#fff;border:1px solid #dfe3ea;border-radius:8px;padding:12px;}",
    sep = "\n"
  )
  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Z-score alignment audit</title>",
    paste0("<style>", css, "</style></head><body><main>"),
    "<h1>Z-score alignment audit</h1>",
    paste0("<div class=\"meta\">", meta, "</div>"),
    check_block,
    "<h2>Results by trait</h2>",
    paste0("<p class=\"muted\">", counts_line, "</p>"),
    qc_html_table(summary[summary_cols], max_rows = 1000L),
    "<h2>Sign concordance by trait</h2>",
    paste0("<div class=\"plot\">", concordance_plot_svg(summary), "</div>"),
    "<p class=\"muted\">Red band: at most 20% agree (FLIPPED). Amber band: 20&ndash;80% (REVIEW).</p>",
    "<h2>Variant calls by trait</h2>",
    details[nzchar(details)],
    "<h2>Method</h2>",
    methods,
    "</main></body></html>"
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, file, useBytes = TRUE)
  message(sprintf("Audit report written to %s (%d item(s) to check).", file, length(flags)))
  invisible(flags)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript tests/run_tests.R`
Expected: `0 failed or errored.`

- [ ] **Step 5: Commit**

```bash
git add scripts/audit_zscore_alignment.R tests/test-zscore-audit.R
git commit -m "Write the z-score audit HTML report with a concordance plot"
```

---

### Task 6: Manifest, draft mode and the audit driver

**Files:**
- Modify: `scripts/audit_zscore_alignment.R` (append)
- Modify: `tests/test-zscore-audit.R` (append)

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `audit_read_matrix(path) -> numeric matrix` (validated).
  - `read_audit_manifest(path, traits) -> data.frame(trait, efo_id, accession, reported_trait)` (all character, blanks as `""`).
  - `draft_audit_manifest(traits, out_file, search = catalog_search_efo) -> invisible(data.frame)`; columns `trait, efo_id, accession, reported_trait, suggestions`.
  - `audit_zscore_alignment(z, alignment, manifest, output_prefix = NULL, min_abs_z = 2, p_threshold = AUDIT_GWS_P, get_page = catalog_get_page) -> invisible(list(summary, comparisons))`; writes `<prefix>.csv`, `<prefix>_traits.csv`, `<prefix>.html`, `<prefix>_catalog_cache/` when `output_prefix` is given.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test-zscore-audit.R`:

```r
write_lines_tsv <- function(lines) {
  f <- tempfile(fileext = ".tsv")
  writeLines(lines, f)
  f
}

test_that("read_audit_manifest fills optional columns and validates entries", {
  traits <- c("good", "flipped")
  m <- read_audit_manifest(write_lines_tsv(c("trait\tefo_id", "good\tEFO_0000001", "flipped\t")),
                           traits)
  expect_equal(m$efo_id, c("EFO_0000001", ""))
  expect_equal(m$accession, c("", ""))
  expect_equal(m$reported_trait, c("", ""))
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tefo_id", "nope\tEFO_1")), traits),
               "not in the matrix")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tefo_id", "good\tbad id")), traits),
               "EFO ID")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tefo_id\taccession",
                                                     "good\tEFO_1\tXYZ")), traits),
               "accession")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tefo_id\treported_trait",
                                                     "good\tEFO_1\t(")), traits),
               "reported_trait")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait", "good")), traits), "efo_id")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tefo_id", "good\tEFO_1",
                                                     "good\tEFO_2")), traits), "Duplicate")
})

test_that("draft_audit_manifest writes suggestions and never overwrites", {
  out <- tempfile(fileext = ".tsv")
  search <- function(q) {
    if (q == "Waist circumference male") {
      data.frame(efo_id = c("EFO_1", "EFO_2"), efo_trait = c("waist circumference", "waist"))
    } else {
      data.frame(efo_id = character(0), efo_trait = character(0))
    }
  }
  quietly(draft_audit_manifest(c("Waist circumference_male", "XINSG30"), out, search = search))
  m <- utils::read.delim(out, colClasses = "character", na.strings = character(0))
  expect_equal(names(m), c("trait", "efo_id", "accession", "reported_trait", "suggestions"))
  expect_equal(m$efo_id, c("EFO_1", ""))
  expect_equal(m$suggestions[1], "EFO_1 (waist circumference); EFO_2 (waist)")
  expect_error(draft_audit_manifest("x", out, search = search), "already exists")
  expect_equal(nrow(read_audit_manifest(out, c("Waist circumference_male", "XINSG30"))), 2L)
})

test_that("audit_zscore_alignment labels every trait and writes its outputs", {
  fx <- audit_fixture()
  api <- fake_api(setNames(list(assoc_page(fx$records, 1)),
                           catalog_associations_url("EFO_0000001", 0)))
  get_page <- function(url) {
    if (grepl("EFO_0000666", url, fixed = TRUE)) stop("HTTP 500")
    api$get(url)
  }
  manifest <- data.frame(
    trait = c("good", "flipped", "scrambled", "lookup", "blank"),
    efo_id = c("EFO_0000001", "EFO_0000001", "EFO_0000001", "EFO_0000001; EFO_0000666", ""),
    accession = "", reported_trait = "", stringsAsFactors = FALSE
  )
  prefix <- file.path(tempfile(), "audit")
  res <- quietly(audit_zscore_alignment(fx$z, fx$alignment, manifest, prefix,
                                        get_page = get_page))
  s <- res$summary
  expect_equal(s$verdict[match(manifest$trait, s$trait)],
               c("PASS", "FLIPPED", "REVIEW", "LOOKUP_ERROR", "NOT_CHECKED"))
  expect_match(s$note[s$trait == "lookup"], "HTTP 500")
  expect_length(api$calls, 1L)
  for (suffix in c(".csv", "_traits.csv", ".html")) {
    expect_true(file.exists(paste0(prefix, suffix)))
  }
  expect_true(file.exists(file.path(paste0(prefix, "_catalog_cache"), "EFO_0000001.rds")))
  comp <- utils::read.csv(paste0(prefix, ".csv"))
  expect_equal(sort(unique(comp$trait)), c("flipped", "good", "scrambled"))
  expect_equal(nrow(comp), 21L)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript tests/run_tests.R`
Expected: the 3 new blocks fail with `could not find function "read_audit_manifest"` / `"draft_audit_manifest"` / `"audit_zscore_alignment"`.

- [ ] **Step 3: Write the implementation**

Append to `scripts/audit_zscore_alignment.R`:

```r
AUDIT_COMPARISON_COLS <- c(
  "trait", "rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "association_id", "accession",
  "same_study", "reported_trait", "p_value", "effect_allele", "strand_complemented",
  "catalog_direction", "expected_sign", "matrix_z", "n_catalog_records", "informative",
  "concordant", "status"
)

audit_read_matrix <- function(path) {
  z <- readRDS(path)
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z) || !is.numeric(z) || is.null(rownames(z)) || is.null(colnames(z)) ||
      anyDuplicated(rownames(z)) || anyDuplicated(colnames(z))) {
    stop(path, " is not a numeric matrix with unique row and column names.")
  }
  z
}

audit_split_ids <- function(x) {
  ids <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  ids[nzchar(ids)]
}

read_audit_manifest <- function(path, traits) {
  m <- utils::read.delim(path, colClasses = "character", na.strings = character(0),
                         quote = "", check.names = FALSE)
  missing_cols <- setdiff(c("trait", "efo_id"), names(m))
  if (length(missing_cols)) {
    stop(path, " lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  for (col in c("accession", "reported_trait")) if (!col %in% names(m)) m[[col]] <- ""
  m <- m[c("trait", "efo_id", "accession", "reported_trait")]
  m[] <- lapply(m, function(x) trimws(ifelse(is.na(x), "", x)))

  if (anyDuplicated(m$trait)) {
    stop("Duplicate traits in ", path, ": ", qc_name_list(unique(m$trait[duplicated(m$trait)])))
  }
  unknown <- setdiff(m$trait, traits)
  if (length(unknown)) stop("Traits not in the matrix: ", paste(unknown, collapse = ", "))
  ids <- unlist(lapply(m$efo_id, audit_split_ids))
  bad_ids <- ids[!grepl("^[A-Za-z]+_[0-9]+$", ids)]
  if (length(bad_ids)) stop("Invalid EFO ID(s): ", paste(bad_ids, collapse = ", "))
  bad_acc <- m$accession[nzchar(m$accession) & !grepl("^GCST[0-9]+$", m$accession)]
  if (length(bad_acc)) {
    stop("Invalid accession(s) (expected GCST followed by digits): ", paste(bad_acc, collapse = ", "))
  }
  for (re in m$reported_trait[nzchar(m$reported_trait)]) {
    ok <- tryCatch({ grepl(re, "", perl = TRUE); TRUE }, error = function(e) FALSE)
    if (!ok) stop("Invalid reported_trait regular expression: ", re)
  }
  m
}

draft_audit_manifest <- function(traits, out_file, search = catalog_search_efo) {
  # Suggest EFO terms from each trait name. The result needs human review.
  if (file.exists(out_file)) {
    stop(out_file, " already exists; edit it, or choose another file name.")
  }
  rows <- lapply(traits, function(trait) {
    query <- trimws(gsub("[_.]+", " ", trait))
    hits <- tryCatch(search(query), error = function(e) {
      warning(trait, ": EFO search failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    found <- !is.null(hits) && nrow(hits) > 0L
    message(sprintf("%s: %s", trait, if (found) hits$efo_id[1] else "no EFO match"))
    data.frame(
      trait = trait,
      efo_id = if (found) hits$efo_id[1] else "",
      accession = "", reported_trait = "",
      suggestions = if (found) paste(sprintf("%s (%s)", hits$efo_id, hits$efo_trait),
                                     collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(out)
}

audit_zscore_alignment <- function(z, alignment, manifest, output_prefix = NULL,
                                   min_abs_z = 2, p_threshold = AUDIT_GWS_P,
                                   get_page = catalog_get_page) {
  a <- prepare_alignment(alignment, z)
  cache_dir <- if (!is.null(output_prefix)) paste0(output_prefix, "_catalog_cache")

  # Fetch each EFO term once, whichever traits share it.
  efo_ids <- unique(unlist(lapply(manifest$efo_id, audit_split_ids)))
  catalog <- list()
  errors <- character(0)
  for (id in efo_ids) {
    message("Fetching GWAS Catalog associations for ", id)
    got <- tryCatch(catalog_fetch_efo(id, cache_dir, get_page), error = function(e) e)
    if (inherits(got, "error")) {
      errors[[id]] <- conditionMessage(got)
      warning(id, ": ", conditionMessage(got), call. = FALSE)
    } else {
      catalog[[id]] <- got
    }
  }

  summaries <- list()
  comparisons <- list()
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, ]
    ids <- audit_split_ids(row$efo_id)
    calls <- NULL
    note <- ""
    if (!length(ids)) {
      v <- trait_verdict(NULL)
      v$verdict <- "NOT_CHECKED"
      note <- "No EFO ID given."
    } else if (any(ids %in% names(errors))) {
      v <- trait_verdict(NULL)
      v$verdict <- "LOOKUP_ERROR"
      note <- paste(unname(errors[intersect(ids, names(errors))]), collapse = "; ")
    } else {
      records <- do.call(rbind, unname(catalog[ids]))
      records <- records[!duplicated(records$association_id), , drop = FALSE]
      aligned <- align_catalog_to_matrix(records, a, row$accession, row$reported_trait,
                                         p_threshold)
      calls <- call_snp_trait(aligned, z[, row$trait], min_abs_z)
      v <- trait_verdict(calls, row$accession)
      if (!is.null(calls)) {
        comparisons[[length(comparisons) + 1L]] <- cbind(trait = row$trait, calls)
      }
    }
    message(sprintf("%s: %s (%d of %d informative SNPs concordant)", row$trait, v$verdict,
                    as.integer(v$n_concordant), as.integer(v$n_informative)))
    summaries[[i]] <- cbind(row, v, note = note, stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  comparisons <- if (length(comparisons)) do.call(rbind, comparisons)[AUDIT_COMPARISON_COLS]

  if (!is.null(output_prefix)) {
    dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(if (is.null(comparisons)) data.frame() else comparisons,
                     paste0(output_prefix, ".csv"), row.names = FALSE)
    utils::write.csv(summary, paste0(output_prefix, "_traits.csv"), row.names = FALSE)
    write_audit_report(summary, comparisons, paste0(output_prefix, ".html"),
                       settings = list(min_abs_z = min_abs_z, p_threshold = p_threshold,
                                       n_variants = nrow(z), n_traits = nrow(manifest)))
  }
  invisible(list(summary = summary, comparisons = comparisons))
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript tests/run_tests.R`
Expected: `0 failed or errored.`

- [ ] **Step 5: Commit**

```bash
git add scripts/audit_zscore_alignment.R tests/test-zscore-audit.R
git commit -m "Add the audit manifest, draft mode and audit driver"
```

---

### Task 7: Command line, README and a live run on a real matrix

**Files:**
- Modify: `scripts/audit_zscore_alignment.R` (append the command-line block)
- Modify: `README.md` (pipeline scripts list, outputs table, new section after "Stability analysis")

**Interfaces:**
- Consumes: `audit_read_matrix()`, `read_audit_manifest()`, `draft_audit_manifest()`, `audit_zscore_alignment()`, `AUDIT_VERDICT_ORDER`.
- Produces: the CLI documented in the script header.

- [ ] **Step 1: Append the command-line block**

Append to `scripts/audit_zscore_alignment.R`:

```r
if (sys.nframe() == 0L) {
  # Rscript writes spaces in the script path as "~+~".
  script_path <- gsub("~+~", " ", sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), fixed = TRUE)
  source(file.path(dirname(normalizePath(script_path)), "qc_report.R"))

  usage <- paste(
    "Usage:",
    "  Rscript audit_zscore_alignment.R --draft zmat_preImputed.rds audit_traits.tsv",
    "  Rscript audit_zscore_alignment.R zmat_preImputed.rds alignment_GWAS_summStats.csv \\",
    "    audit_traits.tsv output_prefix [--min-abs-z 2]",
    "See the header of this script for the manifest columns.", sep = "\n")
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[1] %in% c("-h", "--help")) {
    cat(usage, "\n")
    quit(status = 0)
  }
  if (args[1] == "--draft") {
    if (length(args) != 3L) { cat(usage, "\n"); quit(status = 1) }
    z <- audit_read_matrix(args[2])
    draft_audit_manifest(colnames(z), args[3])
    message("Wrote ", args[3], ". Review efo_id for every trait before running the audit.")
  } else {
    min_abs_z <- 2
    flag <- which(args == "--min-abs-z")
    if (length(flag)) {
      min_abs_z <- suppressWarnings(as.numeric(args[flag[1] + 1L]))
      args <- args[-c(flag[1], flag[1] + 1L)]
    }
    if (length(args) != 4L || !is.finite(min_abs_z) || min_abs_z < 0) {
      cat(usage, "\n")
      quit(status = 1)
    }
    z <- audit_read_matrix(args[1])
    alignment <- utils::read.csv(args[2], colClasses = "character", check.names = FALSE)
    manifest <- read_audit_manifest(args[3], colnames(z))
    result <- audit_zscore_alignment(z, alignment, manifest, args[4], min_abs_z = min_abs_z)
    print(table(factor(result$summary$verdict, levels = AUDIT_VERDICT_ORDER)))
  }
}
```

- [ ] **Step 2: Check help output from a path with spaces**

Run: `Rscript "$PWD/scripts/audit_zscore_alignment.R" --help`
Expected: the usage text, exit status 0 (the repository path contains spaces, so this also checks the `~+~` decoding).

Run: `Rscript tests/run_tests.R`
Expected: `0 failed or errored.` (sourcing the file in tests must not run the CLI block).

- [ ] **Step 3: Live run on the T2DGGI matrix (network; not committed)**

Outputs go to the session scratchpad, not the repository. Set, in the same shell for every command below:

```bash
export RES="/Users/kirk/Partners HealthCare Dropbox/Kirk Smith/MGH/Collaborations/T2DGGI_2025/bNMF/pheno_only_v3_results"
export SCRATCH="<the session scratchpad directory>/zscore_audit_live"
mkdir -p "$SCRATCH"
```

First confirm the alignment file has rsIDs:

Run: `Rscript -e 'a <- read.csv(file.path(Sys.getenv("RES"), "alignment_GWAS_summStats.csv")); print(names(a)); print(mean(grepl("^rs[0-9]+$", a$rsID)))'`
Expected: includes `ChrPos, rsID, REF, ALT, Risk_Allele`; share of valid rsIDs near 1.

Write a six-trait manifest `$SCRATCH/audit_six.tsv`:

```
trait	efo_id	accession	reported_trait
BMI	EFO_0004340		
HDL	EFO_0004612		
LDL	EFO_0004611		
HEIGHT	EFO_0004339		
Urate	EFO_0004531		
SBP	EFO_0006335		
```

Run: `Rscript scripts/audit_zscore_alignment.R "$RES/zmat_preImputed.rds" "$RES/alignment_GWAS_summStats.csv" "$SCRATCH/audit_six.tsv" "$SCRATCH/live/audit"`
Expected: six verdicts printed; BMI, HDL, LDL, HEIGHT, Urate and SBP should be PASS. If any is REVIEW or FLIPPED, open `$SCRATCH/live/audit.html` and inspect its SNP table before changing code — the finding may be real.

Negative control — flip one column and confirm it is caught:

Run: `Rscript -e 'z <- readRDS(file.path(Sys.getenv("RES"), "zmat_preImputed.rds")); z[, "HDL"] <- -z[, "HDL"]; saveRDS(z, file.path(Sys.getenv("SCRATCH"), "zmat_flipHDL.rds"))'`
Run: `Rscript scripts/audit_zscore_alignment.R "$SCRATCH/zmat_flipHDL.rds" "$RES/alignment_GWAS_summStats.csv" "$SCRATCH/audit_six.tsv" "$SCRATCH/flip/audit"`
Expected: HDL is FLIPPED; the other five unchanged. (A second run reuses nothing — different prefix — so this also exercises a fresh fetch; to save time, copy `$SCRATCH/live/audit_catalog_cache` to `$SCRATCH/flip/audit_catalog_cache` first.)

Draft-mode smoke test:

Run: `Rscript scripts/audit_zscore_alignment.R --draft "$RES/zmat_preImputed.rds" "$SCRATCH/draft.tsv"`
Expected: one line per trait (94), a TSV with 94 rows; spot-check that `BMI` and `HDL` got sensible suggestions.

- [ ] **Step 4: Update the README**

In `README.md`, under "### Pipeline scripts", after the `qc_report.R` entry, add:

```markdown
#### `audit_zscore_alignment.R` — Allele-alignment audit against the GWAS Catalog
* Command-line script; see "Checking z-score alignment" below. Checks each trait column of `zmat_preImputed.rds` against genome-wide-significant GWAS Catalog associations, independently of your own summary-statistics files.
```

In the "### Outputs" table, add rows:

```markdown
| `zscore_audit.html`, `zscore_audit.csv`, `zscore_audit_traits.csv` | Only when you run the alignment audit: per-trait sign agreement with the GWAS Catalog, and every variant comparison behind it |
```

After the "### Stability analysis" section (before "### Contributors"), add:

````markdown
### Checking z-score alignment

The most likely place for an allele error is the preprocessing of trait summary statistics (for example, using the wrong column as the effect allele), which flips or scrambles a whole trait. `audit_zscore_alignment.R` checks every trait column of `zmat_preImputed.rds` against the NHGRI-EBI GWAS Catalog, so it does not depend on the files that may contain the error. It needs `httr2` and network access; the rsIDs in `alignment_GWAS_summStats.csv` come from LDlink, so LD pruning must have been on.

1. Draft a manifest with one suggested EFO term per trait:

   ```sh
   Rscript scripts/audit_zscore_alignment.R --draft my_results/zmat_preImputed.rds audit_traits.tsv
   ```

2. Review `audit_traits.tsv`. Correct `efo_id` using the `suggestions` column or the [Catalog trait search](https://www.ebi.ac.uk/gwas/) (separate several IDs with `;`; leave blank to skip a trait). Optionally add the `accession` (GCST…) of the study you used, and a `reported_trait` regular expression to narrow a broad term, for example to BMI-adjusted analyses.

3. Run the audit:

   ```sh
   Rscript scripts/audit_zscore_alignment.R my_results/zmat_preImputed.rds \
     my_results/alignment_GWAS_summStats.csv audit_traits.tsv my_results/zscore_audit
   ```

For each trait, the variants with a genome-wide-significant Catalog association (any study; the clustering study is preferred when its accession is given) are compared in sign with the matrix, after orienting the Catalog effect allele to the T2D risk allele. A/T and C/G variants are skipped, and so are variants with matrix |z| < 2 (`--min-abs-z`). A trait is PASS when at least 80% of at least 5 informative variants agree, FLIPPED when at most 20% do, and REVIEW in between. A REVIEW can also mean the EFO term pools a differently defined phenotype; the report's per-trait tables show which variants disagree. Catalog downloads are cached in `zscore_audit_catalog_cache/`; delete it to refresh.
````

- [ ] **Step 5: Run the full test suite**

Run: `Rscript tests/run_tests.R`
Expected: `0 failed or errored.`

- [ ] **Step 6: Commit**

```bash
git add scripts/audit_zscore_alignment.R README.md
git commit -m "Add the audit command line and document the alignment check"
```
