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
