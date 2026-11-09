# Per-run QC report: write_qc_report().

suppressWarnings(suppressPackageStartupMessages(
  source(file.path(scripts_dir, "qc_report.R"))
))

fx_qc <- make_discovery_fixtures()
vars_sig_qc <- quietly(get_sig_snps(gwas = fx_qc$gwas, PVCUTOFF = 1e-6))
vars_main_qc <- quietly(get_biggest_gwas(fx_qc$paths[["primary"]], vars_sig_qc))
tf_qc <- make_trait_fixtures()
fetched_qc <- quietly(fetch_summary_stats(
  df_input = tf_qc$df_input, trait_ss_files = tf_qc$files, trait_ss_size = tf_qc$sizes
))

full_qc <- list(
  run_name = "unit_test_run",
  settings = list(PVCUTOFF = 5e-8, MAIN_PVCUTOFF = 0.05, PROXY_MIN_R2 = 0.8),
  input_qc = attr(vars_sig_qc, "input_qc"),
  concordance = attr(vars_main_qc, "concordance"),
  dropped = attr(vars_main_qc, "dropped"),
  drop_summary = attr(vars_main_qc, "drop_summary"),
  sentinel_pvalue_cutoff = 5e-8,
  funnel = data.frame(
    stage = c("Discovery variants in the pool", "Kept by the primary-GWAS filters",
              "Sentinel candidates", "Final variants"),
    n = c(12L, 5L, 5L, 4L)
  ),
  trait_qc = fetched_qc$trait_qc,
  trait_filters = data.frame(
    trait = c("trait1", "trait2", "trait3"),
    result = c("trait kept", "removed (correlation)", "trait kept"),
    note = c(NA, "correlated w/ trait1", NA)
  ),
  sentinel_nonmissingness = c("1_1000_A_G" = 1, "1_2000_C_T" = 0.5, "2_5000_A_C" = 0),
  nonmissing_cutoff = 0.8,
  proxies_needed = data.frame(VAR_ID = c("1_2000_C_T", "2_5000_A_C"),
                              reason = c("high_missingness", "high_missingness")),
  proxies_found = data.frame(original_SNP = "1_2000_C_T", VAR_ID = "1_2100_A_G", r2 = 0.93),
  final_variants = data.frame(VAR_ID = c("1_1000_A_G", "1_2100_A_G"),
                              Variant_Type = c("original", "proxy")),
  bnmf_run_summary = data.frame(run = 1:4, K = c(5, 5, 6, 5), evid = c(-1, -2, -3, -4))
)

report_file <- tempfile("qc_report_", fileext = ".html")
flags <- quietly(write_qc_report(full_qc, report_file))
html <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

test_that("the report is written as one self-contained HTML file", {
  expect_true(file.exists(report_file))
  expect_match(html, "<html", fixed = TRUE)
  expect_match(html, "unit_test_run", fixed = TRUE)
  expect_false(grepl("<script src=|<link rel=\"stylesheet\"", html))
})

test_that("each supplied part of the run appears as a section", {
  for (heading in c("Run settings", "Discovery GWAS inputs",
                    "Agreement with the primary GWAS",
                    "Variants removed by the primary-GWAS filters",
                    "Variant counts by stage", "Trait GWAS inputs",
                    "Trait filters", "Sentinel trait coverage",
                    "Proxy search", "Final variant set", "bNMF")) {
    expect_match(html, heading, fixed = TRUE)
  }
})

test_that("the tables carry the run's actual content", {
  expect_match(html, "discB_EAS", fixed = TRUE)
  expect_match(html, "risk_allele_mismatch", fixed = TRUE)
  expect_match(html, "assumed ALT", fixed = TRUE)
  expect_match(html, "removed (correlation)", fixed = TRUE)
})

test_that("removed variants that were sentinel-level are listed by name", {
  # 1:4000 (discovery P 1e-8 < 5e-8) was removed for a risk-allele mismatch.
  expect_match(html, "1:4000:G:T", fixed = TRUE)
})

test_that("the report opens with the things that need checking", {
  expect_true(is.character(flags))
  expect_true(any(grepl("assumed", flags)))                 # effect allele assumed ALT
  expect_true(any(grepl("strand complement", flags)))       # trait2 matched by complement
  expect_true(any(grepl("no data in any trait", flags)))    # 2_5000_A_C has coverage 0
  expect_true(any(grepl("no proxy", flags)))                # 2_5000_A_C got no proxy
  expect_match(html, "Things to check", fixed = TRUE)
})

test_that("low agreement with the primary GWAS is flagged", {
  low <- full_qc
  low$concordance$n_testable <- 100L
  low$concordance$n_concordant <- c(100L, 4L, 100L, 100L)[seq_len(nrow(low$concordance))]
  low$concordance$concordance <- low$concordance$n_concordant / low$concordance$n_testable
  low_flags <- quietly(write_qc_report(low, tempfile(fileext = ".html")))
  expect_true(any(grepl("agreement", low_flags)))
})

test_that("a report with only some parts still renders, without the missing sections", {
  partial_file <- tempfile("qc_partial_", fileext = ".html")
  partial_flags <- quietly(write_qc_report(
    list(run_name = "partial", input_qc = full_qc$input_qc), partial_file
  ))
  partial <- paste(readLines(partial_file, warn = FALSE), collapse = "\n")
  expect_match(partial, "Discovery GWAS inputs", fixed = TRUE)
  expect_false(grepl("Proxy search", partial, fixed = TRUE))
  expect_false(grepl("Trait GWAS inputs", partial, fixed = TRUE))
})

test_that("text from the data is escaped", {
  odd <- full_qc
  odd$trait_filters$trait[2] <- "HDL<C> & friends"   # the removed trait, listed by name
  odd_file <- tempfile(fileext = ".html")
  quietly(write_qc_report(odd, odd_file))
  odd_html <- paste(readLines(odd_file, warn = FALSE), collapse = "\n")
  expect_match(odd_html, "HDL&lt;C&gt; &amp; friends", fixed = TRUE)
})

test_that("proxy candidate outcomes and unanswered LDlink queries appear in the report", {
  with_proxy_audit <- full_qc
  with_proxy_audit$proxy_candidates <- data.frame(
    VAR_ID = c("1_2000_C_T", "1_2000_C_T", "2_5000_A_C"),
    proxy_ChrPos = c("1:2100", "1:2200", "2:5100"),
    status = c("chosen", "main_pvalue_above_cutoff", "not_in_main_gwas")
  )
  with_proxy_audit$proxies_not_found <- data.frame(
    VAR_ID = "2_5000_A_C", reason = "no_candidate_passed_primary_filters"
  )
  with_proxy_audit$proxy_unanswered <- "3_7000_A_G"
  audit_file <- tempfile(fileext = ".html")
  audit_flags <- quietly(write_qc_report(with_proxy_audit, audit_file))
  audit_html <- paste(readLines(audit_file, warn = FALSE), collapse = "\n")
  expect_match(audit_html, "Proxy candidates by outcome", fixed = TRUE)
  expect_match(audit_html, "main_pvalue_above_cutoff", fixed = TRUE)
  expect_match(audit_html, "no_candidate_passed_primary_filters", fixed = TRUE)
  expect_true(any(grepl("LDlink returned no answer", audit_flags)))
})
