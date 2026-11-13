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
  d <- quietly(catalog_fetch_efo("EFO_0000001", cache_dir = cache, get_page = api$get))
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
