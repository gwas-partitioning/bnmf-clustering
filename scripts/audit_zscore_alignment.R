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
  bad_risk[is.na(bad_risk)] <- TRUE
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
