# =============================================================================
# audit_zscore_alignment.R
# =============================================================================
# Checks that each trait column of the pre-imputed z-score matrix
# (zmat_preImputed.rds) is oriented to the T2D risk allele, by comparing it
# with the same variants in the Common Metabolic Diseases Knowledge Portal
# (https://hugeamp.org), which holds full summary statistics (every variant,
# not only significant hits) for about 1,700 phenotypes.
#
# The check is independent of your own summary-statistics files, so it can
# catch an allele mistake made while preprocessing them (for example, the
# wrong effect-allele column), which flips or scrambles a whole trait.
#
# Usage:
#   1. Draft the audit manifest (a Portal phenotype suggested per trait):
#        Rscript scripts/audit_zscore_alignment.R --draft \
#          results/zmat_preImputed.rds audit_traits.tsv
#   2. Review audit_traits.tsv: correct kp_phenotype (see the suggestions
#      column, or https://hugeamp.org for the phenotype list) or blank it to
#      skip a trait.
#   3. Run the audit:
#        Rscript scripts/audit_zscore_alignment.R \
#          results/zmat_preImputed.rds results/alignment_GWAS_summStats.csv \
#          audit_traits.tsv results/zscore_audit \
#          [--ancestry Mixed] [--min-abs-z 2] [--workers 4]
#      Writes results/zscore_audit.csv (one row per variant-trait comparison),
#      results/zscore_audit_traits.csv (one row per trait),
#      results/zscore_audit.html, and results/zscore_audit_portal_cache/.
#      The first run queries the Portal once per variant (all phenotypes come
#      back together), about 10-20 minutes for 1,000 variants; reruns use the
#      cache. --ancestry selects the Portal's ancestry-specific analysis
#      (for example EU, AF, EA, SA, HS); the default is its multi-ancestry
#      meta-analysis.
#
# Verdicts: PASS (>= 80% of informative variants agree in sign), FLIPPED
# (<= 20%), REVIEW (in between), INSUFFICIENT (< 5 informative variants),
# NOT_CHECKED (no Portal phenotype). A variant is informative when both the
# matrix and the Portal z-score have |z| >= 2.
#
# From R, source scripts/qc_report.R before this file.
#
# Dependencies: httr2 (and jsonlite, which it uses); parallel (base R).
# =============================================================================

KP_API <- "https://bioindex.hugeamp.org/api"
AUDIT_USER_AGENT <- "bnmf_clustering z-score audit"

kp_get_json <- function(url) {
  # The Portal refuses requests without a user agent.
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Install the httr2 package to query the Knowledge Portal.")
  }
  httr2::request(url) |>
    httr2::req_user_agent(AUDIT_USER_AGENT) |>
    httr2::req_timeout(180) |>
    httr2::req_retry(max_tries = 4) |>
    httr2::req_perform() |>
    httr2::resp_body_json()
}

kp_variant_url <- function(var_id, ancestry = "Mixed") {
  if (identical(ancestry, "Mixed")) {
    sprintf("%s/bio/query/phewas-associations?q=%s", KP_API, var_id)
  } else {
    sprintf("%s/bio/query/ancestry-phewas-associations?q=%s,%s", KP_API, ancestry, var_id)
  }
}

audit_or <- function(x, y) if (is.null(x) || !length(x)) y else x

kp_parse_rows <- function(rows) {
  # Portal rows (parsed JSON) to phenotype, beta, se, p, n.
  num <- function(f) vapply(rows, function(r) {
    suppressWarnings(as.numeric(audit_or(r[[f]], NA)[[1]]))
  }, numeric(1))
  data.frame(
    phenotype = vapply(rows, function(r) as.character(audit_or(r$phenotype, NA)[[1]]),
                       character(1)),
    beta = num("beta"), se = num("stdErr"), p = num("pValue"), n = num("n"),
    stringsAsFactors = FALSE
  )
}

kp_empty <- function() {
  data.frame(phenotype = character(0), beta = numeric(0), se = numeric(0), p = numeric(0),
             n = numeric(0), effect_allele = character(0), other_allele = character(0),
             kp_var_id = character(0), stringsAsFactors = FALSE)
}

kp_fetch_variant <- function(chrom, pos, ref, alt, ancestry = "Mixed", cache_dir = NULL,
                             get_json = kp_get_json) {
  # All Portal phenotypes at one variant. The Portal names a variant
  # chr:pos:A1:A2 (hg19) and reports beta for A2; its allele order need not
  # match ours, so the other order is tried when the first returns nothing.
  cache_file <- if (!is.null(cache_dir)) {
    file.path(cache_dir, sprintf("%s_%s_%s_%s_%s.rds", chrom, pos, ref, alt, ancestry))
  }
  if (!is.null(cache_file) && file.exists(cache_file)) return(readRDS(cache_file))
  out <- NULL
  for (alleles in list(c(ref, alt), c(alt, ref))) {
    var_id <- paste(chrom, pos, alleles[1], alleles[2], sep = ":")
    response <- get_json(kp_variant_url(var_id, ancestry))
    if (is.null(response$count)) stop("Unexpected Knowledge Portal response for ", var_id)
    if (response$count > 0) {
      out <- kp_parse_rows(response$data)
      out$effect_allele <- alleles[2]
      out$other_allele <- alleles[1]
      out$kp_var_id <- var_id
      break
    }
  }
  if (is.null(out)) out <- kp_empty()
  if (!is.null(cache_file)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_file)
  }
  out
}

kp_fetch_matrix <- function(alignment, ancestry = "Mixed", cache_dir = NULL,
                            get_json = kp_get_json, workers = 4L) {
  # Portal results for every matrix variant, fetched in parallel where the
  # platform allows forking. Variants whose query fails are listed in
  # attr(, "failed") and left out; rerunning retries them.
  pos <- sub("^.*:", "", alignment$ChrPos)
  chrom <- sub("^chr", "", sub(":.*$", "", alignment$ChrPos), ignore.case = TRUE)
  one <- function(i) {
    tryCatch({
      d <- kp_fetch_variant(chrom[i], pos[i], alignment$REF[i], alignment$ALT[i], ancestry,
                            cache_dir, get_json)
      if (nrow(d)) d$ChrPos <- alignment$ChrPos[i] else d$ChrPos <- character(0)
      d
    }, error = function(e) structure(conditionMessage(e), class = "kp_failure"))
  }
  idx <- seq_len(nrow(alignment))
  workers <- if (.Platform$OS.type == "windows") 1L else max(1L, as.integer(workers))
  message(sprintf("Querying the Knowledge Portal for %d variants (%s, %d worker%s)...",
                  length(idx), ancestry, workers, if (workers == 1L) "" else "s"))
  results <- if (workers > 1L) {
    parallel::mclapply(idx, one, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    lapply(idx, one)
  }
  failed <- vapply(results, function(r) inherits(r, "kp_failure") || inherits(r, "try-error"),
                   logical(1))
  out <- do.call(rbind, results[!failed])
  if (is.null(out)) out <- cbind(kp_empty(), ChrPos = character(0))
  rownames(out) <- NULL
  attr(out, "failed") <- alignment$ChrPos[failed]
  attr(out, "not_found") <- setdiff(alignment$ChrPos[!failed], out$ChrPos)
  if (any(failed)) {
    warning(sum(failed), " variant(s) could not be fetched (first: ", alignment$ChrPos[failed][1],
            ": ", unclass(results[failed][[1]]), "). Rerun to retry them.", call. = FALSE)
  }
  out
}

kp_phenotypes <- function(get_json = kp_get_json) {
  rows <- get_json(sprintf("%s/portal/phenotypes", KP_API))$data
  data.frame(
    name = vapply(rows, function(r) as.character(audit_or(r$name, NA)[[1]]), character(1)),
    description = vapply(rows, function(r) as.character(audit_or(r$description, "")[[1]]), character(1)),
    group = vapply(rows, function(r) as.character(audit_or(r$group, "")[[1]]), character(1)),
    stringsAsFactors = FALSE
  )
}

audit_norm_label <- function(x) trimws(gsub("[^a-z0-9]+", " ", tolower(x)))

suggest_kp_phenotypes <- function(phenotypes, trait, n = 3L) {
  # Up to n Portal phenotypes for a matrix trait name, best first:
  #   0 the phenotype code equals the trait name (ignoring case and punctuation)
  #   1 the description has the same words as the trait name
  #   2 the description has every word of the trait name (and more)
  #   3 the trait name has every word of the description
  # Shorter descriptions first within a tier.
  empty <- phenotypes[0, c("name", "description")]
  q <- audit_norm_label(trait)
  if (!nzchar(q) || !nrow(phenotypes)) return(empty)
  squash <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  qw <- unique(strsplit(q, " ", fixed = TRUE)[[1]])
  desc <- audit_norm_label(sub("\\s*\\([^)]*\\)\\s*$", "", phenotypes$description))
  dw <- strsplit(desc, " ", fixed = TRUE)
  has_query <- vapply(dw, function(w) all(qw %in% w), logical(1))
  in_query <- vapply(dw, function(w) length(w) > 0L && all(w %in% qw), logical(1)) &
    nchar(desc) >= 4L
  tier <- ifelse(squash(phenotypes$name) == squash(trait), 0L,
          ifelse(has_query & in_query, 1L,
          ifelse(has_query, 2L,
          ifelse(in_query, 3L, NA_integer_))))
  keep <- which(!is.na(tier))
  if (!length(keep)) return(empty)
  keep <- keep[order(tier[keep], nchar(desc[keep]))]
  out <- utils::head(phenotypes[keep, c("name", "description")], n)
  rownames(out) <- NULL
  out
}

prepare_alignment <- function(alignment, z) {
  # One alignment row per matrix row, in matrix order.
  needed <- c("ChrPos", "REF", "ALT", "Risk_Allele")
  missing_cols <- setdiff(needed, names(alignment))
  if (length(missing_cols)) {
    stop("The alignment file lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  # Runs from earlier pipeline versions name matrix rows "chr1:123" while the
  # alignment file has "1:123"; match on the position without the prefix.
  strip_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
  idx <- match(strip_chr(rownames(z)), strip_chr(alignment$ChrPos))
  if (anyNA(idx)) {
    stop(sum(is.na(idx)), " matrix row(s) have no alignment row (first: ",
         rownames(z)[is.na(idx)][1], ").")
  }
  keep_cols <- intersect(c(needed, "rsID"), names(alignment))
  a <- as.data.frame(alignment, stringsAsFactors = FALSE)[idx, keep_cols, drop = FALSE]
  a[] <- lapply(a, as.character)
  if (!"rsID" %in% names(a)) a$rsID <- NA_character_
  a$ChrPos <- rownames(z)
  a[c("REF", "ALT", "Risk_Allele")] <- lapply(a[c("REF", "ALT", "Risk_Allele")], toupper)
  bad_risk <- is.na(a$Risk_Allele) | !(a$Risk_Allele == a$REF | a$Risk_Allele == a$ALT)
  bad_risk[is.na(bad_risk)] <- TRUE
  if (any(bad_risk)) {
    stop(sum(bad_risk), " variant(s) have a Risk_Allele that is neither REF nor ALT (first: ",
         a$ChrPos[bad_risk][1], ").")
  }
  rownames(a) <- NULL
  a[c("ChrPos", "rsID", "REF", "ALT", "Risk_Allele")]
}

compare_trait <- function(portal, alignment, z_values, phenotype, min_abs_z = 2) {
  # One row per matrix variant that the Portal has for this phenotype. The
  # Portal z-score (beta / se for its effect allele) is re-expressed for the
  # T2D risk allele and compared in sign with the matrix z-score.
  p <- portal[portal$phenotype == phenotype & is.finite(portal$beta) &
                is.finite(portal$se) & portal$se > 0, , drop = FALSE]
  p <- p[!duplicated(p$ChrPos), , drop = FALSE]
  if (!nrow(p)) return(NULL)
  d <- merge(alignment, p, by = "ChrPos")
  if (!nrow(d)) return(NULL)
  aligned_to_risk <- d$effect_allele == d$Risk_Allele
  d$portal_z <- ifelse(aligned_to_risk, 1, -1) * d$beta / d$se
  d$matrix_z <- unname(z_values[d$ChrPos])
  d$status <- ifelse(!is.finite(d$matrix_z), "missing_in_matrix",
              ifelse(abs(d$matrix_z) < min_abs_z | abs(d$portal_z) < min_abs_z, "weak",
              ifelse(sign(d$matrix_z) == sign(d$portal_z), "concordant", "discordant")))
  d$informative <- d$status %in% c("concordant", "discordant")
  d$concordant <- ifelse(d$informative, d$status == "concordant", NA)
  d <- d[order(match(d$ChrPos, names(z_values))), , drop = FALSE]
  rownames(d) <- NULL
  d
}

AUDIT_VERDICT_ORDER <- c("FLIPPED", "REVIEW", "INSUFFICIENT", "PASS", "NOT_CHECKED")

trait_verdict <- function(calls, min_n = 5L, pass_at = 0.8, flip_at = 0.2) {
  informative <- if (is.null(calls)) logical(0) else calls$informative
  n <- sum(informative)
  k <- if (n) sum(calls$concordant[informative]) else 0L
  concordance <- if (n) k / n else NA_real_
  both <- if (is.null(calls)) logical(0) else is.finite(calls$matrix_z) & is.finite(calls$portal_z)
  z_cor <- if (sum(both) >= 3L) suppressWarnings(stats::cor(calls$matrix_z[both], calls$portal_z[both])) else NA_real_
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
    n_portal_variants = if (is.null(calls)) 0L else nrow(calls),
    n_informative = n,
    n_concordant = k,
    concordance = concordance,
    binom_p = if (n) stats::binom.test(k, n)$p.value else NA_real_,
    z_correlation = z_cor,
    stringsAsFactors = FALSE
  )
}

AUDIT_COLORS <- c(PASS = "#2f7d4f", REVIEW = "#b26b00", FLIPPED = "#b3261e")

audit_setting <- function(settings, name, default) {
  if (is.null(settings[[name]])) default else settings[[name]]
}

concordance_plot_svg <- function(summary) {
  # Dot plot: one row per trait with informative variants, x = share
  # concordant, dot area grows with the number of informative variants.
  s <- summary[summary$n_informative > 0 & is.finite(summary$concordance), , drop = FALSE]
  if (!nrow(s)) return("<p class=\"muted\">No trait had informative variants to plot.</p>")
  s <- s[order(s$concordance, s$trait), , drop = FALSE]
  left <- 240; plot_w <- 480; right <- 60; top <- 16; row_h <- 18
  width <- left + plot_w + right
  bottom <- top + nrow(s) * row_h
  height <- bottom + 40
  f <- function(v) formatC(v, format = "f", digits = 1)
  x <- function(v) left + v * plot_w
  y <- top + (seq_len(nrow(s)) - 0.5) * row_h
  r <- pmin(8, 2.5 + sqrt(s$n_informative) / 2)
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
                   'Share of informative variants whose matrix sign matches the Portal</text>'),
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
  summary <- summary[order(match(summary$verdict, AUDIT_VERDICT_ORDER), summary$trait), ,
                     drop = FALSE]
  count <- table(factor(summary$verdict, levels = AUDIT_VERDICT_ORDER))
  traits_with <- function(v) qc_name_list(summary$trait[summary$verdict == v])

  flags <- character(0)
  if (count[["FLIPPED"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) look flipped: the matrix sign disagrees with the Portal at most informative variants (%s). Check which column was used as the effect allele when these GWAS were formatted.",
    count[["FLIPPED"]], traits_with("FLIPPED")))
  if (count[["REVIEW"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) have mixed agreement (%s). Possible causes: a partial allele mix-up, or a Portal phenotype defined differently from your GWAS (for example adjusted versus unadjusted). Their variant tables are open below.",
    count[["REVIEW"]], traits_with("REVIEW")))
  negative <- summary$trait[summary$verdict %in% c("PASS", "INSUFFICIENT") &
                              is.finite(summary$z_correlation) & summary$z_correlation < 0]
  if (length(negative)) flags <- c(flags, sprintf(
    "Matrix and Portal z-scores are negatively correlated for: %s. Check these traits even though their verdict is not FLIPPED.",
    qc_name_list(negative)))
  unknown <- summary$trait[grepl("not a Portal phenotype", summary$note, fixed = TRUE)]
  if (length(unknown)) flags <- c(flags, sprintf(
    "Phenotype codes not found in the Portal for: %s.", qc_name_list(unknown)))
  n_failed <- length(audit_setting(settings, "failed", character(0)))
  if (n_failed) flags <- c(flags, sprintf(
    "%d variant(s) could not be fetched from the Portal and were left out; rerun to retry them.",
    n_failed))

  check_block <- if (length(flags)) {
    paste0("<div class=\"checks\"><h2>Things to check</h2><ul>",
           paste0("<li>", qc_html_escape(flags), "</li>", collapse = ""), "</ul></div>")
  } else {
    "<div class=\"checks ok\"><h2>Things to check</h2><p>No flags were raised.</p></div>"
  }

  summary_cols <- c("trait", "verdict", "n_informative", "n_concordant", "concordance",
                    "binom_p", "z_correlation", "n_portal_variants", "kp_phenotype", "note")
  detail_cols <- c("rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "kp_var_id", "effect_allele",
                   "matrix_z", "portal_z", "status", "p", "n")
  details <- vapply(seq_len(nrow(summary)), function(i) {
    s <- summary[i, ]
    rows <- if (is.null(comparisons)) NULL else comparisons[comparisons$trait == s$trait, , drop = FALSE]
    if (is.null(rows) || !nrow(rows)) return("")
    sprintf("<details%s><summary><strong>%s</strong> &middot; %s (%d of %d informative variants concordant)</summary>%s</details>",
            if (s$verdict %in% c("FLIPPED", "REVIEW")) " open" else "",
            qc_html_escape(s$trait), s$verdict, as.integer(s$n_concordant),
            as.integer(s$n_informative), qc_html_table(rows[detail_cols], max_rows = 1500L))
  }, character(1))

  meta <- sprintf("%s variants &times; %s traits &middot; Knowledge Portal, %s analysis &middot; generated %s",
                  audit_setting(settings, "n_variants", "?"),
                  audit_setting(settings, "n_traits", nrow(summary)),
                  qc_html_escape(audit_setting(settings, "ancestry", "Mixed")),
                  format(Sys.time(), "%Y-%m-%d %H:%M"))
  counts_line <- paste(sprintf("%s: %d", names(count), as.integer(count)), collapse = " &middot; ")
  methods <- sprintf(paste(
    "<p>Each matrix variant was looked up in the Common Metabolic Diseases Knowledge Portal (hugeamp.org) by hg19",
    "position and alleles, which returns that variant's beta and standard error for every Portal phenotype. For each",
    "trait, the Portal z-score (beta / SE) for the phenotype named in the manifest was expressed for the T2D risk allele",
    "and compared with the matrix z-score. Positions and alleles are matched exactly on the same build, so strand-ambiguous",
    "(A/T, C/G) variants are included.</p>",
    "<p>A variant is informative when both z-scores have |z| &ge; %s. Verdicts need at least 5 informative variants:",
    "PASS when at least 80%% agree in sign, FLIPPED when at most 20%% do, REVIEW in between. The binomial P tests",
    "agreement against 50%%; z_correlation is the Pearson correlation over all variants the Portal has for the trait.</p>",
    "<p>Limitations: a Portal phenotype is usually a meta-analysis across the Portal's datasets and may differ from your",
    "GWAS in definition (adjustment, sex) or ancestry; such differences lower agreement without any error on your side.",
    "If your trait GWAS was itself downloaded from the Portal, the check is not independent of it.</p>"),
    format(min_abs_z))

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
    "<h2>Variant comparisons by trait</h2>",
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
  "trait", "kp_phenotype", "rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "kp_var_id",
  "effect_allele", "beta", "se", "p", "n", "portal_z", "matrix_z", "informative",
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

read_audit_manifest <- function(path, traits) {
  m <- utils::read.delim(path, colClasses = "character", na.strings = character(0),
                         quote = "", check.names = FALSE)
  missing_cols <- setdiff(c("trait", "kp_phenotype"), names(m))
  if (length(missing_cols)) {
    stop(path, " lacks column(s): ", paste(missing_cols, collapse = ", "),
         ". Draft a new manifest with --draft.")
  }
  m <- m[c("trait", "kp_phenotype")]
  m[] <- lapply(m, function(x) trimws(ifelse(is.na(x), "", x)))
  if (anyDuplicated(m$trait)) {
    stop("Duplicate traits in ", path, ": ", qc_name_list(unique(m$trait[duplicated(m$trait)])))
  }
  unknown <- setdiff(m$trait, traits)
  if (length(unknown)) stop("Traits not in the matrix: ", paste(unknown, collapse = ", "))
  bad <- m$kp_phenotype[nzchar(m$kp_phenotype) & grepl("[[:space:],;]", m$kp_phenotype)]
  if (length(bad)) stop("Give one Portal phenotype code per trait (got: ", paste(bad, collapse = ", "), ").")
  m
}

draft_audit_manifest <- function(traits, out_file, phenotypes) {
  # Suggest a Portal phenotype for each trait name. The result needs review.
  if (file.exists(out_file)) {
    stop(out_file, " already exists; edit it, or choose another file name.")
  }
  rows <- lapply(traits, function(trait) {
    hits <- suggest_kp_phenotypes(phenotypes, gsub("[_.]+", " ", trait))
    if (!nrow(hits)) hits <- suggest_kp_phenotypes(phenotypes, trait)
    found <- nrow(hits) > 0L
    message(sprintf("%s: %s", trait,
                    if (found) sprintf("%s (%s)", hits$name[1], hits$description[1]) else "no match"))
    data.frame(
      trait = trait,
      kp_phenotype = if (found) hits$name[1] else "",
      suggestions = if (found) paste(sprintf("%s (%s)", hits$name, hits$description),
                                     collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(out)
}

audit_zscore_alignment <- function(z, alignment, manifest, portal, output_prefix = NULL,
                                   min_abs_z = 2, phenotype_names = NULL,
                                   settings = list()) {
  a <- prepare_alignment(alignment, z)
  summaries <- list()
  comparisons <- list()
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, ]
    note <- ""
    calls <- NULL
    if (!nzchar(row$kp_phenotype)) {
      v <- trait_verdict(NULL)
      v$verdict <- "NOT_CHECKED"
      note <- "No Portal phenotype given."
    } else {
      if (!is.null(phenotype_names) && !row$kp_phenotype %in% phenotype_names) {
        note <- sprintf("%s is not a Portal phenotype code.", row$kp_phenotype)
      }
      calls <- compare_trait(portal, a, z[, row$trait], row$kp_phenotype, min_abs_z)
      v <- trait_verdict(calls)
      if (!is.null(calls)) {
        calls$trait <- row$trait
        calls$kp_phenotype <- row$kp_phenotype
        comparisons[[length(comparisons) + 1L]] <- calls[AUDIT_COMPARISON_COLS]
      }
    }
    message(sprintf("%s: %s (%d of %d informative variants concordant)", row$trait, v$verdict,
                    as.integer(v$n_concordant), as.integer(v$n_informative)))
    summaries[[i]] <- cbind(row, v, note = note, stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  comparisons <- if (length(comparisons)) do.call(rbind, comparisons)

  if (!is.null(output_prefix)) {
    dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(if (is.null(comparisons)) data.frame() else comparisons,
                     paste0(output_prefix, ".csv"), row.names = FALSE)
    utils::write.csv(summary, paste0(output_prefix, "_traits.csv"), row.names = FALSE)
    settings$min_abs_z <- min_abs_z
    settings$n_variants <- nrow(z)
    settings$n_traits <- nrow(manifest)
    write_audit_report(summary, comparisons, paste0(output_prefix, ".html"), settings)
  }
  invisible(list(summary = summary, comparisons = comparisons))
}

if (sys.nframe() == 0L) {
  # Rscript writes spaces in the script path as "~+~".
  script_path <- gsub("~+~", " ", sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), fixed = TRUE)
  source(file.path(dirname(normalizePath(script_path)), "qc_report.R"))

  usage <- paste(
    "Usage:",
    "  Rscript audit_zscore_alignment.R --draft zmat_preImputed.rds audit_traits.tsv",
    "  Rscript audit_zscore_alignment.R zmat_preImputed.rds alignment_GWAS_summStats.csv \\",
    "    audit_traits.tsv output_prefix [--ancestry Mixed] [--min-abs-z 2] [--workers 4]",
    "See the header of this script for the manifest columns.", sep = "\n")
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat(usage, "\n")
    quit(status = 0)
  }
  take_value <- function(flag, default) {
    at <- which(args == flag)
    if (!length(at)) return(default)
    if (at[1] == length(args)) { cat(usage, "\n"); quit(status = 1) }
    value <- args[at[1] + 1L]
    args <<- args[-c(at[1], at[1] + 1L)]
    value
  }
  ancestry <- take_value("--ancestry", "Mixed")
  min_abs_z <- suppressWarnings(as.numeric(take_value("--min-abs-z", "2")))
  workers <- suppressWarnings(as.integer(take_value("--workers", "4")))
  draft <- "--draft" %in% args
  args <- setdiff(args, "--draft")
  if (!is.finite(min_abs_z) || min_abs_z < 0 || is.na(workers) ||
      !grepl("^[A-Za-z]+$", ancestry) || length(args) != if (draft) 2L else 4L) {
    cat(usage, "\n")
    quit(status = 1)
  }
  phenotypes <- kp_phenotypes()
  if (draft) {
    z <- audit_read_matrix(args[1])
    draft_audit_manifest(colnames(z), args[2], phenotypes)
    message("Wrote ", args[2], ". Review kp_phenotype for every trait before running the audit.")
  } else {
    z <- audit_read_matrix(args[1])
    alignment <- utils::read.csv(args[2], colClasses = "character", check.names = FALSE)
    manifest <- read_audit_manifest(args[3], colnames(z))
    a <- prepare_alignment(alignment, z)
    portal <- kp_fetch_matrix(a, ancestry, paste0(args[4], "_portal_cache"), workers = workers)
    result <- audit_zscore_alignment(
      z, alignment, manifest, portal, args[4], min_abs_z = min_abs_z,
      phenotype_names = phenotypes$name,
      settings = list(ancestry = ancestry, failed = attr(portal, "failed"))
    )
    print(table(factor(result$summary$verdict, levels = AUDIT_VERDICT_ORDER)))
  }
}
