# Z-score alignment audit: scripts/audit_zscore_alignment.R
# No network access: Catalog data are small tables in the layout of the bulk
# "All associations with added ontology annotations" file.

suppressWarnings(suppressPackageStartupMessages({
  source(file.path(scripts_dir, "qc_report.R"))
  source(file.path(scripts_dir, "audit_zscore_alignment.R"))
}))

efo_uri <- function(id) paste0("http://www.ebi.ac.uk/efo/", id)

cat_row <- function(rs, ea, p = "1E-10", value = "0.1", ci = "[0.05-0.15] unit increase",
                    uri = efo_uri("EFO_0000001"), mapped = "trait x",
                    accession = "GCST000001", reported = "Trait X", snps = rs,
                    current = "", risk = paste0(rs, "-", ea)) {
  row <- list(risk, snps, current, p, value, ci, mapped, uri, accession, reported)
  names(row) <- unname(CATALOG_COLUMNS)
  as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
}

fake_catalog <- function(...) parse_catalog_rows(do.call(rbind, list(...)))

test_that("catalog_direction reads CI text first, then odds ratios above 1", {
  expect_equal(
    catalog_direction(c("[0.05-0.15] unit increase", "[0.02-0.02] mmol/L decrease",
                        "[1.05-1.13]", "", "[0.8-0.9]", NA, "[1.8-2.14] unit decrease"),
                      c("0.1", "0.02", "1.09", "", "0.85", "2", "1.97")),
    c(1, -1, 1, NA, NA, 1, -1)
  )
})

test_that("parse_catalog_rows reads alleles, rsIDs and ontology terms", {
  d <- parse_catalog_rows(rbind(
    cat_row("rs1", "a", p = "2E-9"),
    cat_row("rs2", "T", risk = "rs2 -T"),
    cat_row("rs3", "", risk = "rs3--"),
    cat_row("rs4", "?"),
    cat_row("rs5; rs6", "A", risk = "rs5-A"),
    cat_row("chr6:31892560", "G"),
    cat_row("rs7", "C", current = "777"),
    cat_row("rs8", "G", uri = paste(efo_uri("EFO_0000002"),
                                    "http://purl.obolibrary.org/obo/OBA_VT0005416", sep = ", ")),
    cat_row("rs9", "G", risk = "rs10-G")
  ))
  expect_equal(d$catalog_row, as.character(1:9))
  expect_equal(d$effect_allele, c("A", "T", NA, "?", NA, "G", "C", "G", NA))
  expect_equal(d$single_snp, c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_true(is.na(d$rsID[5]) && is.na(d$rsID[6]))
  expect_equal(d$rsID_current[7], "rs777")
  expect_equal(d$term_ids[c(1, 8)], c("EFO_0000001", "EFO_0000002,OBA_VT0005416"))
  expect_equal(d$p_value[1], 2e-9)
  expect_error(parse_catalog_rows(data.frame(SNPS = "rs1")), "lacks column")
})

test_that("read_catalog_bulk reads the TSV and the zip", {
  raw <- rbind(cat_row("rs1", "A"), cat_row("rs2", "G", value = "", ci = ""))
  raw$extra <- "ignored"
  dir <- tempfile()
  dir.create(dir)
  tsv <- file.path(dir, "assoc.tsv")
  utils::write.table(raw, tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  from_tsv <- quietly(read_catalog_bulk(tsv))
  expect_equal(from_tsv, parse_catalog_rows(raw[unname(CATALOG_COLUMNS)]),
               ignore_attr = TRUE)
  expect_true(is.na(from_tsv$effect_value[2]))
  skip_if(!nzchar(Sys.which("zip")) || !nzchar(Sys.which("unzip")), "zip/unzip not available")
  zip <- file.path(dir, "assoc.zip")
  utils::zip(zip, tsv, flags = "-jq")
  expect_equal(quietly(read_catalog_bulk(zip)), from_tsv)
  expect_error(read_catalog_bulk(file.path(dir, "missing.zip")), "not found")
})

test_that("catalog_for_terms matches whole IDs only", {
  d <- fake_catalog(cat_row("rs1", "A", uri = efo_uri("EFO_1")),
                    cat_row("rs2", "A", uri = efo_uri("EFO_10")),
                    cat_row("rs3", "A", uri = paste(efo_uri("EFO_2"), efo_uri("EFO_1"), sep = ", ")))
  expect_equal(catalog_for_terms(d, "EFO_1")$rsID, c("rs1", "rs3"))
  expect_equal(catalog_for_terms(d, c("EFO_10", "EFO_2"))$rsID, c("rs2", "rs3"))
})

test_that("suggest_efo_terms ranks label matches and counts associations", {
  d <- fake_catalog(
    cat_row("rs1", "A", uri = efo_uri("EFO_1"), mapped = "body mass index"),
    cat_row("rs2", "A", uri = efo_uri("EFO_1"), mapped = "body mass index"),
    cat_row("rs3", "A", uri = efo_uri("EFO_2"), mapped = "BMI-adjusted waist-hip ratio"),
    cat_row("rs4", "A", uri = paste(efo_uri("EFO_3"), efo_uri("EFO_1"), sep = ", "),
            mapped = "urate measurement, body mass index"),
    cat_row("rs5", "A", uri = efo_uri("EFO_4"), mapped = "glucose measurement"),
    cat_row("rs6", "A", uri = efo_uri("EFO_5"), mapped = "glucose tolerance test")
  )
  terms <- catalog_terms(d)
  expect_equal(terms$n[terms$efo_id == "EFO_1"], 3L)
  expect_equal(suggest_efo_terms(terms, "Body mass index")$efo_id, "EFO_1")
  expect_equal(suggest_efo_terms(terms, "Urate")$efo_id, "EFO_3")
  expect_equal(suggest_efo_terms(terms, "Glucose")$efo_id, c("EFO_4", "EFO_5"))
  # Abbreviation expanded; the literal "BMI" match ranks after it.
  expect_equal(suggest_efo_terms(terms, "BMI")$efo_id, c("EFO_1", "EFO_2"))
  # A label contained in the trait name ("urate" in "Urate_male").
  expect_equal(suggest_efo_terms(terms, "Urate male")$efo_id, "EFO_3")
  expect_equal(nrow(suggest_efo_terms(terms, "XINSG30")), 0L)
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
  raw <- rbind(
    cat_row("rs1", "G"),                                       # risk allele, increase
    cat_row("rs2", "T", ci = "[0.1-0.3] unit decrease"),       # non-risk, decrease
    cat_row("rs3", "T"),                                       # complement of risk A
    cat_row("rs4", "T", value = "1.2", ci = "[1.1-1.3]"),      # risk, odds ratio
    cat_row("rs5", "A", value = "1.3", ci = ""),               # risk, odds ratio, no CI
    cat_row("rs6", "A"),
    cat_row("rs7", "T"),                                       # palindromic: dropped
    cat_row("rs8", "G")
  )
  z <- cbind(good = rep(5, 8), flipped = rep(-5, 8),
             scrambled = c(5, -5, 5, -5, 5, -5, 5, -5),
             unknown = rep(5, 8), blank = rep(5, 8))
  rownames(z) <- alignment$ChrPos
  list(alignment = alignment, catalog = parse_catalog_rows(raw), z = z)
}

test_that("prepare_alignment orders rows like the matrix and flags palindromes", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment[8:1, ], fx$z)
  expect_equal(a$ChrPos, rownames(fx$z))
  expect_equal(a$palindromic, c(rep(FALSE, 6), TRUE, FALSE))
  z_chr <- fx$z
  rownames(z_chr) <- paste0("chr", rownames(z_chr))
  expect_equal(prepare_alignment(fx$alignment, z_chr)$ChrPos, rownames(z_chr))
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
  catalog <- fake_catalog(
    cat_row("rs1", "G", p = "1E-6"),                         # not genome-wide
    cat_row("rs2", "?"),                                      # unknown allele
    cat_row("rs3", "AT"),                                     # not a single base
    cat_row("rs4", "T", value = "", ci = ""),                 # no direction
    cat_row("rs5", "G", reported = "BMI-adjusted thing"),
    cat_row("rs6", "A", reported = "Unadjusted thing"),
    cat_row("rs99", "A")                                      # not in the matrix
  )
  d <- align_catalog_to_matrix(catalog, a, reported_trait = "unadjusted")
  expect_equal(d$rsID, "rs6")
  expect_equal(nrow(align_catalog_to_matrix(catalog[0, ], a)), 0L)
})

test_that("align_catalog_to_matrix matches a merged rsID by its current ID", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment, fx$z)
  d <- align_catalog_to_matrix(fake_catalog(cat_row("rs555", "A", current = "6")), a)
  expect_equal(d$rsID, "rs6")
  expect_equal(d$expected_sign, 1)
})

test_that("call_snp_trait prefers the clustering study and flags conflicts", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment, fx$z)
  catalog <- fake_catalog(
    cat_row("rs1", "G", p = "1E-30", accession = "GCST000008"),
    cat_row("rs1", "G", p = "1E-9", accession = "GCST000009"),
    cat_row("rs2", "C", p = "1E-12"),                         # risk allele: +1
    cat_row("rs2", "T", p = "1E-10"),                         # non-risk: -1
    cat_row("rs3", "A", p = "1E-10")
  )
  d <- align_catalog_to_matrix(catalog, a, accession = "GCST000009")
  z <- setNames(c(5, 5, 1, rep(5, 5)), rownames(fx$z))
  calls <- call_snp_trait(d, z, min_abs_z = 2)
  expect_equal(calls$rsID, c("rs1", "rs2", "rs3"))
  expect_equal(calls$catalog_row, c("2", "3", "5"))
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
                                      settings = list(min_abs_z = 2, catalog_label = "test.zip")))
  html <- paste(readLines(f), collapse = "\n")
  for (t in c(traits, "blank")) expect_match(html, paste0(">", t, "<"), fixed = TRUE)
  expect_match(html, "<svg", fixed = TRUE)
  expect_match(html, "<details open>", fixed = TRUE)
  expect_match(html, "test.zip", fixed = TRUE)
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
  m <- read_audit_manifest(write_lines_tsv(c("trait\tefo_id", "good\tOBA_VT0005416", "flipped\t")),
                           traits)
  expect_equal(m$efo_id, c("OBA_VT0005416", ""))
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
  catalog <- fake_catalog(
    cat_row("rs1", "A", uri = efo_uri("EFO_1"), mapped = "waist circumference male"),
    cat_row("rs2", "A", uri = efo_uri("EFO_2"), mapped = "waist circumference male adjusted")
  )
  quietly(draft_audit_manifest(c("Waist circumference_male", "XINSG30"), out, catalog))
  m <- utils::read.delim(out, colClasses = "character", na.strings = character(0))
  expect_equal(names(m), c("trait", "efo_id", "accession", "reported_trait", "suggestions"))
  expect_equal(m$efo_id, c("EFO_1", ""))
  expect_equal(m$suggestions[1],
               "EFO_1 (waist circumference male, n=1); EFO_2 (waist circumference male adjusted, n=1)")
  expect_error(draft_audit_manifest("x", out, catalog), "already exists")
  expect_equal(nrow(read_audit_manifest(out, c("Waist circumference_male", "XINSG30"))), 2L)
})

test_that("audit_zscore_alignment labels every trait and writes its outputs", {
  fx <- audit_fixture()
  manifest <- data.frame(
    trait = c("good", "flipped", "scrambled", "unknown", "blank"),
    efo_id = c("EFO_0000001", "EFO_0000001", "EFO_0000001", "EFO_0000666", ""),
    accession = "", reported_trait = "", stringsAsFactors = FALSE
  )
  prefix <- file.path(tempfile(), "audit")
  res <- quietly(audit_zscore_alignment(fx$z, fx$alignment, manifest, fx$catalog, prefix,
                                        catalog_label = "test.zip"))
  s <- res$summary
  expect_equal(s$verdict[match(manifest$trait, s$trait)],
               c("PASS", "FLIPPED", "REVIEW", "INSUFFICIENT", "NOT_CHECKED"))
  expect_match(s$note[s$trait == "unknown"], "not in the Catalog file")
  for (suffix in c(".csv", "_traits.csv", ".html")) {
    expect_true(file.exists(paste0(prefix, suffix)))
  }
  comp <- utils::read.csv(paste0(prefix, ".csv"))
  expect_equal(sort(unique(comp$trait)), c("flipped", "good", "scrambled"))
  expect_equal(nrow(comp), 21L)
})
