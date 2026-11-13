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
