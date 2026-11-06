#!/usr/bin/env Rscript

# Run the pipeline unit tests:
#   Rscript tests/run_tests.R
#
# The tests build small synthetic summary-statistics files in a temporary
# directory; they need no network access and no LDlink token.

suppressPackageStartupMessages(library(testthat))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
tests_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
} else {
  normalizePath("tests", mustWork = TRUE)
}
Sys.setenv(BNMF_REPO_DIR = normalizePath(file.path(tests_dir, ".."), mustWork = TRUE))

results <- as.data.frame(
  testthat::test_dir(tests_dir, reporter = "summary", stop_on_failure = FALSE)
)
n_failed <- sum(results$failed) + sum(results$error)
cat(sprintf("\n%d test blocks, %d failed or errored.\n", nrow(results), n_failed))
if (n_failed > 0L) quit(status = 1L)
