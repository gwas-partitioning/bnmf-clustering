#!/usr/bin/env Rscript

# Render the versioned bNMF interpretation report.
#
# Usage:
#   Rscript render_bNMF_report_INTERPRETATION_V2.R RESULTS_DIR [PREVIOUS_RESULTS_DIR]

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: Rscript render_bNMF_report_INTERPRETATION_V2.R RESULTS_DIR [PREVIOUS_RESULTS_DIR]")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
previous_main_dir <- if (length(args) >= 2 && nzchar(args[2])) {
  normalizePath(args[2], mustWork = TRUE)
} else {
  NULL
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_arg)) stop("Could not determine this script's directory.")
report_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))

report_file <- file.path(report_dir, "format_bNMF_results_2025_INTERPRETATION_V2.Rmd")
helper_file <- file.path(report_dir, "make_cluster_figures_INTERPRETATION_V2.R")
posthoc_file <- file.path(main_dir, "cluster_post_hoc.rds")

if (!file.exists(report_file) || !file.exists(helper_file)) {
  stop("The V2 report and helper must remain together in: ", report_dir)
}
if (!file.exists(posthoc_file)) stop("Missing: ", posthoc_file)

posthoc <- readRDS(posthoc_file)
run_name <- basename(main_dir)
output_file <- sprintf("%s_K%i_interpretation_v2.html", run_name, posthoc$k)

rmarkdown::render(
  input = report_file,
  output_file = output_file,
  output_dir = main_dir,
  params = list(
    main_dir = main_dir,
    k = NULL,
    previous_main_dir = previous_main_dir,
    helper_script = helper_file
  ),
  envir = new.env(parent = globalenv())
)

message("Report written to: ", file.path(main_dir, output_file))
