# Post-hoc functions that do not need annotation databases: the W-only weight
# cutoff. post_bNMF_2025.R is read as text and only the cutoff function is
# evaluated, so the Bioconductor annotation packages are not needed.

post_src <- readLines(file.path(scripts_dir, "post_bNMF_2025.R"))
cutoff_start <- grep("^calculate_cutoff <- function", post_src)[1]
cutoff_end <- grep("^calculate_h_core_traits <- function", post_src)[1] - 1L
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
eval(parse(text = post_src[cutoff_start:cutoff_end]))

make_w <- function(n_variants = 300L, n_clusters = 10L, n_strong = 25L) {
  set.seed(1)
  w <- as.data.frame(matrix(rexp(n_variants * n_clusters, rate = 4),
                            n_variants, n_clusters))
  names(w) <- paste0("X", seq_len(n_clusters))
  w[seq_len(n_strong), 1:3] <- w[seq_len(n_strong), 1:3] + 3
  w
}

test_that("the weight cutoff uses W only and a top-0.2% anchor with a minimum of 25", {
  w <- make_w()
  out_dir <- tempfile("cutoff_")
  cutoff <- suppressWarnings(calculate_cutoff(w, output_dir = out_dir))
  diag <- attr(cutoff, "diagnostics")
  expect_true(is.numeric(cutoff) && length(cutoff) == 1L)
  expect_gt(cutoff, 0)
  expect_lt(cutoff, max(as.matrix(w)))
  expect_equal(diag$n_pooled_weights, 3000L)
  expect_equal(diag$upper_anchor_n, 25L)          # 0.2% of 3000 is 6, so the minimum applies
  expect_true(file.exists(file.path(out_dir, "bNMF_cutoff_diagnostics.csv")))
  # A larger matrix: 0.2% of 30000 values is 60, above the minimum.
  big <- suppressWarnings(calculate_cutoff(make_w(3000L), output_dir = tempfile("cutoff_")))
  expect_equal(attr(big, "diagnostics")$upper_anchor_n, 60L)
})

test_that("H is not used for the cutoff, and negative weights are refused", {
  w <- make_w()
  a <- suppressWarnings(calculate_cutoff(w, h = NULL, output_dir = tempfile("cutoff_")))
  b <- suppressWarnings(calculate_cutoff(w, h = data.frame(x = 1:3), output_dir = tempfile("cutoff_")))
  expect_identical(as.numeric(a), as.numeric(b))
  w_bad <- w; w_bad[1, 1] <- -1
  expect_error(calculate_cutoff(w_bad, output_dir = tempfile("cutoff_")), "nonnegative")
})

test_that("the public post-hoc code has no Open Targets or enrichment dependencies", {
  post_text <- paste(post_src, collapse = "\n")
  report_files <- list.files(file.path(scripts_dir, "make_bNMF_report"),
                             full.names = TRUE)
  report_text <- paste(unlist(lapply(report_files, readLines, warn = FALSE)), collapse = "\n")
  for (txt in list(post_text, report_text)) {
    expect_false(grepl("OpenTargets|v2g_index|enrichR|STRINGdb|use_string_api", txt))
  }
})
