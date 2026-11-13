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
