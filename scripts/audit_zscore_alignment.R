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

AUDIT_VERDICT_ORDER <- c("FLIPPED", "REVIEW", "LOOKUP_ERROR", "INSUFFICIENT",
                         "PASS", "NOT_CHECKED")

trait_verdict <- function(calls, accession = "", min_n = 5L, pass_at = 0.8,
                          flip_at = 0.2) {
  informative <- if (is.null(calls)) logical(0) else calls$informative
  n <- sum(informative)
  k <- if (n) sum(calls$concordant[informative]) else 0L
  concordance <- if (n) k / n else NA_real_
  same <- if (n) informative & calls$same_study else logical(0)
  n_same <- sum(same)
  verdict <- if (n < min_n) {
    "INSUFFICIENT"
  } else if (concordance >= pass_at) {
    "PASS"
  } else if (concordance <= flip_at) {
    "FLIPPED"
  } else {
    "REVIEW"
  }
  data.frame(
    verdict = verdict,
    n_catalog_snps = if (is.null(calls)) 0L else nrow(calls),
    n_informative = n,
    n_concordant = k,
    concordance = concordance,
    binom_p = if (n) stats::binom.test(k, n)$p.value else NA_real_,
    n_same_study = n_same,
    same_study_concordance = if (n_same) mean(calls$concordant[same]) else NA_real_,
    study_differs = if (nzchar(accession)) n_same == 0L else NA,
    stringsAsFactors = FALSE
  )
}

AUDIT_COLORS <- c(PASS = "#2f7d4f", REVIEW = "#b26b00", FLIPPED = "#b3261e")

audit_setting <- function(settings, name, default) {
  if (is.null(settings[[name]])) default else settings[[name]]
}

concordance_plot_svg <- function(summary) {
  # Dot plot: one row per trait with informative SNPs, x = share concordant,
  # dot area grows with the number of informative SNPs.
  s <- summary[summary$n_informative > 0 & is.finite(summary$concordance), , drop = FALSE]
  if (!nrow(s)) return("<p class=\"muted\">No trait had informative SNPs to plot.</p>")
  s <- s[order(s$concordance, s$trait), , drop = FALSE]
  left <- 240; plot_w <- 480; right <- 60; top <- 16; row_h <- 18
  width <- left + plot_w + right
  bottom <- top + nrow(s) * row_h
  height <- bottom + 40
  f <- function(v) formatC(v, format = "f", digits = 1)
  x <- function(v) left + v * plot_w
  y <- top + (seq_len(nrow(s)) - 0.5) * row_h
  r <- pmin(8, 2.5 + sqrt(s$n_informative) / 1.5)
  col <- unname(AUDIT_COLORS[s$verdict])
  col[is.na(col)] <- "#7a8794"
  ticks <- c(0, 0.2, 0.5, 0.8, 1)
  paste(c(
    sprintf(paste0('<svg viewBox="0 0 %d %s" width="100%%" role="img" ',
                   'aria-label="Sign concordance by trait" ',
                   'style="max-width:%dpx;font-size:11px;font-family:inherit">'),
            width, f(height), width),
    sprintf('<rect x="%s" y="%d" width="%s" height="%s" fill="#b3261e" fill-opacity="0.08"/>',
            f(x(0)), top, f(0.2 * plot_w), f(bottom - top)),
    sprintf('<rect x="%s" y="%d" width="%s" height="%s" fill="#b26b00" fill-opacity="0.07"/>',
            f(x(0.2)), top, f(0.6 * plot_w), f(bottom - top)),
    sprintf(paste0('<line x1="%1$s" x2="%1$s" y1="%2$d" y2="%3$s" stroke="#c9ced8"/>',
                   '<text x="%1$s" y="%4$s" text-anchor="middle" fill="#6a7486">%5$s</text>'),
            f(x(ticks)), top, f(bottom), f(bottom + 14), format(ticks)),
    sprintf(paste0('<text x="%s" y="%s" text-anchor="middle" fill="#3b4658">',
                   'Share of informative SNPs whose matrix sign matches the Catalog</text>'),
            f(x(0.5)), f(bottom + 32)),
    sprintf('<text x="%d" y="%s" text-anchor="end" dominant-baseline="middle" fill="#1c2330">%s</text>',
            left - 8, f(y), qc_html_escape(s$trait)),
    sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s"><title>%s: %d of %d concordant (%s)</title></circle>',
            f(x(s$concordance)), f(y), f(r), col, qc_html_escape(s$trait),
            as.integer(s$n_concordant), as.integer(s$n_informative), s$verdict),
    sprintf('<text x="%s" y="%s" dominant-baseline="middle" fill="#6a7486">n=%d</text>',
            f(x(1) + 12), f(y), as.integer(s$n_informative)),
    "</svg>"
  ), collapse = "\n")
}

write_audit_report <- function(summary, comparisons, file, settings = list()) {
  if (!exists("qc_html_table", mode = "function")) {
    stop("Source scripts/qc_report.R before writing the audit report.")
  }
  min_abs_z <- audit_setting(settings, "min_abs_z", 2)
  p_threshold <- audit_setting(settings, "p_threshold", AUDIT_GWS_P)
  summary <- summary[order(match(summary$verdict, AUDIT_VERDICT_ORDER), summary$trait), ,
                     drop = FALSE]
  count <- table(factor(summary$verdict, levels = AUDIT_VERDICT_ORDER))
  traits_with <- function(v) qc_name_list(summary$trait[summary$verdict == v])

  flags <- character(0)
  if (count[["FLIPPED"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) look flipped: the matrix sign disagrees with the Catalog at most informative SNPs (%s). Check which column was used as the effect allele when these GWAS were formatted.",
    count[["FLIPPED"]], traits_with("FLIPPED")))
  if (count[["REVIEW"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) have mixed agreement (%s). Possible causes: a partial allele mix-up, or an EFO term that pools a differently defined phenotype (for example adjusted and unadjusted). Their SNP tables are open below.",
    count[["REVIEW"]], traits_with("REVIEW")))
  if (count[["LOOKUP_ERROR"]] > 0) flags <- c(flags, sprintf(
    "The Catalog lookup failed for %d trait(s) (%s); rerun to retry.",
    count[["LOOKUP_ERROR"]], traits_with("LOOKUP_ERROR")))
  differs <- summary$trait[summary$verdict %in% c("PASS", "REVIEW", "FLIPPED") &
                             summary$study_differs %in% TRUE]
  if (length(differs)) flags <- c(flags, sprintf(
    "No informative Catalog record came from the clustering study for: %s. Their verdicts rest on other studies of the trait.",
    qc_name_list(differs)))

  check_block <- if (length(flags)) {
    paste0("<div class=\"checks\"><h2>Things to check</h2><ul>",
           paste0("<li>", qc_html_escape(flags), "</li>", collapse = ""), "</ul></div>")
  } else {
    "<div class=\"checks ok\"><h2>Things to check</h2><p>No flags were raised.</p></div>"
  }

  summary_cols <- c("trait", "verdict", "n_informative", "n_concordant", "concordance",
                    "binom_p", "n_catalog_snps", "n_same_study", "same_study_concordance",
                    "study_differs", "efo_id", "accession", "reported_trait", "note")
  detail_cols <- c("rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "effect_allele",
                   "strand_complemented", "catalog_direction", "expected_sign", "matrix_z",
                   "status", "p_value", "accession", "same_study", "association_id",
                   "reported_trait")
  details <- vapply(seq_len(nrow(summary)), function(i) {
    s <- summary[i, ]
    rows <- if (is.null(comparisons)) NULL else comparisons[comparisons$trait == s$trait, , drop = FALSE]
    if (is.null(rows) || !nrow(rows)) return("")
    sprintf("<details%s><summary><strong>%s</strong> &middot; %s (%d of %d informative SNPs concordant)</summary>%s</details>",
            if (s$verdict %in% c("FLIPPED", "REVIEW")) " open" else "",
            qc_html_escape(s$trait), s$verdict, as.integer(s$n_concordant),
            as.integer(s$n_informative), qc_html_table(rows[detail_cols], max_rows = 500L))
  }, character(1))

  meta <- sprintf("%s variants &times; %s traits &middot; generated %s",
                  audit_setting(settings, "n_variants", "?"),
                  audit_setting(settings, "n_traits", nrow(summary)),
                  format(Sys.time(), "%Y-%m-%d %H:%M"))
  counts_line <- paste(sprintf("%s: %d", names(count), as.integer(count)), collapse = " &middot; ")
  methods <- sprintf(paste(
    "<p>For each trait, every GWAS Catalog association annotated with the trait's EFO term(s) was taken from any study,",
    "keeping single-SNP records with P &le; %s and a stated direction (an odds ratio, or a beta described as an increase or decrease).",
    "Records were matched to matrix variants by rsID. A/T and C/G variants were excluded; an effect allele matching only the",
    "complement of REF/ALT was read on the opposite strand. Each record's direction was expressed for the T2D risk allele.</p>",
    "<p>Each variant gets one call: the smallest-P record from the clustering study when its accession is given and present,",
    "otherwise from any study. When those records disagree in direction the call is a Catalog conflict and is not counted.",
    "A call is informative when the matrix |z| &ge; %s. Verdicts need at least 5 informative variants: PASS when at least 80%%",
    "agree in sign, FLIPPED when at most 20%% do, REVIEW in between. The binomial P tests agreement against 50%%.</p>",
    "<p>Limitations: an EFO term can pool related phenotypes whose directions differ (adjusted versus unadjusted, sex-specific);",
    "the reported_trait filter narrows a term. Catalog records are curated top hits, so traits with few significant loci among",
    "these variants may be INSUFFICIENT. Only the sign is checked, not the magnitude.</p>"),
    format(p_threshold), format(min_abs_z))

  css <- paste(
    "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f6f7f9;color:#1c2330;}",
    "main{max-width:1100px;margin:0 auto;padding:24px 16px 64px;}",
    "h1{font-size:24px;margin:8px 0 4px;} h2{font-size:18px;margin:32px 0 8px;}",
    "p{line-height:1.5;margin:6px 0;} .muted{color:#6a7486;font-size:13px;}",
    ".meta{color:#6a7486;font-size:13px;margin-bottom:16px;}",
    ".checks{background:#fff7e6;border:1px solid #f0c36d;border-radius:8px;padding:4px 16px 12px;}",
    ".checks.ok{background:#eef8f0;border-color:#9fd3aa;} .checks h2{margin-top:12px;} .checks li{margin:6px 0;line-height:1.45;}",
    ".tablewrap{overflow-x:auto;background:#fff;border:1px solid #dfe3ea;border-radius:8px;}",
    "table{border-collapse:collapse;width:100%;font-size:13px;}",
    "th,td{padding:6px 10px;border-bottom:1px solid #eceff4;text-align:left;white-space:nowrap;}",
    "th{background:#f0f2f6;font-weight:600;} tr:last-child td{border-bottom:none;}",
    ".num{text-align:right;font-variant-numeric:tabular-nums;}",
    "details{margin:8px 0;} summary{cursor:pointer;padding:4px 0;}",
    ".plot{background:#fff;border:1px solid #dfe3ea;border-radius:8px;padding:12px;}",
    sep = "\n"
  )
  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Z-score alignment audit</title>",
    paste0("<style>", css, "</style></head><body><main>"),
    "<h1>Z-score alignment audit</h1>",
    paste0("<div class=\"meta\">", meta, "</div>"),
    check_block,
    "<h2>Results by trait</h2>",
    paste0("<p class=\"muted\">", counts_line, "</p>"),
    qc_html_table(summary[summary_cols], max_rows = 1000L),
    "<h2>Sign concordance by trait</h2>",
    paste0("<div class=\"plot\">", concordance_plot_svg(summary), "</div>"),
    "<p class=\"muted\">Red band: at most 20% agree (FLIPPED). Amber band: 20&ndash;80% (REVIEW).</p>",
    "<h2>Variant calls by trait</h2>",
    details[nzchar(details)],
    "<h2>Method</h2>",
    methods,
    "</main></body></html>"
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, file, useBytes = TRUE)
  message(sprintf("Audit report written to %s (%d item(s) to check).", file, length(flags)))
  invisible(flags)
}

AUDIT_COMPARISON_COLS <- c(
  "trait", "rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "association_id", "accession",
  "same_study", "reported_trait", "p_value", "effect_allele", "strand_complemented",
  "catalog_direction", "expected_sign", "matrix_z", "n_catalog_records", "informative",
  "concordant", "status"
)

audit_read_matrix <- function(path) {
  z <- readRDS(path)
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z) || !is.numeric(z) || is.null(rownames(z)) || is.null(colnames(z)) ||
      anyDuplicated(rownames(z)) || anyDuplicated(colnames(z))) {
    stop(path, " is not a numeric matrix with unique row and column names.")
  }
  z
}

audit_split_ids <- function(x) {
  ids <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  ids[nzchar(ids)]
}

read_audit_manifest <- function(path, traits) {
  m <- utils::read.delim(path, colClasses = "character", na.strings = character(0),
                         quote = "", check.names = FALSE)
  missing_cols <- setdiff(c("trait", "efo_id"), names(m))
  if (length(missing_cols)) {
    stop(path, " lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  for (col in c("accession", "reported_trait")) if (!col %in% names(m)) m[[col]] <- ""
  m <- m[c("trait", "efo_id", "accession", "reported_trait")]
  m[] <- lapply(m, function(x) trimws(ifelse(is.na(x), "", x)))

  if (anyDuplicated(m$trait)) {
    stop("Duplicate traits in ", path, ": ", qc_name_list(unique(m$trait[duplicated(m$trait)])))
  }
  unknown <- setdiff(m$trait, traits)
  if (length(unknown)) stop("Traits not in the matrix: ", paste(unknown, collapse = ", "))
  ids <- unlist(lapply(m$efo_id, audit_split_ids))
  bad_ids <- ids[!grepl("^[A-Za-z]+_[0-9]+$", ids)]
  if (length(bad_ids)) stop("Invalid EFO ID(s): ", paste(bad_ids, collapse = ", "))
  bad_acc <- m$accession[nzchar(m$accession) & !grepl("^GCST[0-9]+$", m$accession)]
  if (length(bad_acc)) {
    stop("Invalid accession(s) (expected GCST followed by digits): ", paste(bad_acc, collapse = ", "))
  }
  for (re in m$reported_trait[nzchar(m$reported_trait)]) {
    ok <- tryCatch({ suppressWarnings(grepl(re, "", perl = TRUE)); TRUE },
                   error = function(e) FALSE)
    if (!ok) stop("Invalid reported_trait regular expression: ", re)
  }
  m
}

draft_audit_manifest <- function(traits, out_file, search = catalog_search_efo) {
  # Suggest EFO terms from each trait name. The result needs human review.
  if (file.exists(out_file)) {
    stop(out_file, " already exists; edit it, or choose another file name.")
  }
  rows <- lapply(traits, function(trait) {
    query <- trimws(gsub("[_.]+", " ", trait))
    hits <- tryCatch(search(query), error = function(e) {
      warning(trait, ": EFO search failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    found <- !is.null(hits) && nrow(hits) > 0L
    message(sprintf("%s: %s", trait, if (found) hits$efo_id[1] else "no EFO match"))
    data.frame(
      trait = trait,
      efo_id = if (found) hits$efo_id[1] else "",
      accession = "", reported_trait = "",
      suggestions = if (found) paste(sprintf("%s (%s)", hits$efo_id, hits$efo_trait),
                                     collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(out)
}

audit_zscore_alignment <- function(z, alignment, manifest, output_prefix = NULL,
                                   min_abs_z = 2, p_threshold = AUDIT_GWS_P,
                                   get_page = catalog_get_page) {
  a <- prepare_alignment(alignment, z)
  cache_dir <- if (!is.null(output_prefix)) paste0(output_prefix, "_catalog_cache")

  # Fetch each EFO term once, whichever traits share it.
  efo_ids <- unique(unlist(lapply(manifest$efo_id, audit_split_ids)))
  catalog <- list()
  errors <- character(0)
  for (id in efo_ids) {
    message("Fetching GWAS Catalog associations for ", id)
    got <- tryCatch(catalog_fetch_efo(id, cache_dir, get_page), error = function(e) e)
    if (inherits(got, "error")) {
      errors[[id]] <- conditionMessage(got)
      warning(id, ": ", conditionMessage(got), call. = FALSE)
    } else {
      catalog[[id]] <- got
    }
  }

  summaries <- list()
  comparisons <- list()
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, ]
    ids <- audit_split_ids(row$efo_id)
    calls <- NULL
    note <- ""
    if (!length(ids)) {
      v <- trait_verdict(NULL)
      v$verdict <- "NOT_CHECKED"
      note <- "No EFO ID given."
    } else if (any(ids %in% names(errors))) {
      v <- trait_verdict(NULL)
      v$verdict <- "LOOKUP_ERROR"
      note <- paste(unname(errors[intersect(ids, names(errors))]), collapse = "; ")
    } else {
      records <- do.call(rbind, unname(catalog[ids]))
      records <- records[!duplicated(records$association_id), , drop = FALSE]
      aligned <- align_catalog_to_matrix(records, a, row$accession, row$reported_trait,
                                         p_threshold)
      calls <- call_snp_trait(aligned, z[, row$trait], min_abs_z)
      v <- trait_verdict(calls, row$accession)
      if (!is.null(calls)) {
        comparisons[[length(comparisons) + 1L]] <- cbind(trait = row$trait, calls)
      }
    }
    message(sprintf("%s: %s (%d of %d informative SNPs concordant)", row$trait, v$verdict,
                    as.integer(v$n_concordant), as.integer(v$n_informative)))
    summaries[[i]] <- cbind(row, v, note = note, stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  comparisons <- if (length(comparisons)) do.call(rbind, comparisons)[AUDIT_COMPARISON_COLS]

  if (!is.null(output_prefix)) {
    dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(if (is.null(comparisons)) data.frame() else comparisons,
                     paste0(output_prefix, ".csv"), row.names = FALSE)
    utils::write.csv(summary, paste0(output_prefix, "_traits.csv"), row.names = FALSE)
    write_audit_report(summary, comparisons, paste0(output_prefix, ".html"),
                       settings = list(min_abs_z = min_abs_z, p_threshold = p_threshold,
                                       n_variants = nrow(z), n_traits = nrow(manifest)))
  }
  invisible(list(summary = summary, comparisons = comparisons))
}
