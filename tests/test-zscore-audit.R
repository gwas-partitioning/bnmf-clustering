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
