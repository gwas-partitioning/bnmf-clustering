# =============================================================================
# qc_report.R
# =============================================================================
# Per-run QC report for the bNMF pipeline.
#
# Functions:
#   write_qc_report()  — Write one self-contained HTML file summarising how the
#                        inputs were read, which variants and traits were
#                        removed and why, and what needs a second look.
#
# Usage:
#   Collect the pieces of a run in a named list as the run proceeds, then call
#   write_qc_report(qc, file). Every element is optional; a section is written
#   only when its data are present, so the report can be written part-way
#   through a run and again at the end.
#
#   qc <- list(run_name = version, settings = list(PVCUTOFF = PVCUTOFF, ...))
#   qc$input_qc     <- attr(vars_sig,  "input_qc")
#   qc$concordance  <- attr(vars_main, "concordance")
#   qc$dropped      <- attr(vars_main, "dropped")
#   qc$drop_summary <- attr(vars_main, "drop_summary")
#   qc$trait_qc     <- z_n_mats$trait_qc
#   write_qc_report(qc, file.path(main_dir, "QC_report.html"))
#
# Recognised elements of `qc`:
#   run_name                 character
#   settings                 named list or vector of thresholds and options
#   input_qc                 attr(get_sig_snps(), "input_qc")
#   concordance              attr(get_biggest_gwas(), "concordance")
#   sentinels, main_pvalue_cutoff
#                            final sentinel table (VAR_ID, PVALUE_discovery, PVALUE_main,
#                            Risk_Allele, Risk_Allele_discovery) and the MAIN_PVCUTOFF used
#   dropped, drop_summary    attr(get_biggest_gwas(), "dropped" / "drop_summary")
#                            (or fetch_summary_stats()$dropped on the legacy path)
#   sentinel_pvalue_cutoff   discovery P that defines a sentinel; used to list
#                            removed variants that would have been sentinels
#   funnel                   data frame: stage, n  (variant counts by stage)
#   trait_qc                 fetch_summary_stats()$trait_qc
#   trait_filters            data frame: trait, result, note
#   sentinel_nonmissingness  named numeric vector: fraction of traits with data
#   nonmissing_cutoff        threshold below which a sentinel needs a proxy
#   proxies_needed           data frame with VAR_ID and reason
#   proxies_found            data frame with original_SNP (and the proxy chosen)
#   proxies_not_found        select_proxies()$no_proxy: VAR_ID and why no proxy
#   proxy_candidates         select_proxies()$candidates: every LD candidate
#                            with its final status
#   proxy_unanswered         sentinels LDlink gave no answer for
#   final_variants           data frame with VAR_ID and Variant_Type
#   bnmf_run_summary         data frame with K per repetition
#   comparison               list(label, previous, current): character vectors
#                            of variant IDs from an earlier run and this one
#   notes                    character vector of free-text remarks
#
# Dependencies: base R only.
# =============================================================================

qc_html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

qc_format_column <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.logical(v)) return(ifelse(is.na(v), "", ifelse(v, "yes", "no")))
  if (is.numeric(v)) {
    out <- rep("", length(v))
    ok <- !is.na(v)
    # Tiny non-zero values (P-values) must not be mistaken for whole numbers.
    whole <- ok & abs(v - round(v)) < 1e-9 & abs(v) < 1e15 & (v == 0 | abs(v) >= 0.5)
    out[whole] <- format(round(v[whole]), big.mark = ",", scientific = FALSE, trim = TRUE)
    frac <- ok & !whole
    out[frac] <- formatC(v[frac], digits = 3, format = "g")
    return(out)
  }
  out <- as.character(v)
  out[is.na(out)] <- ""
  qc_html_escape(out)
}

qc_html_table <- function(df, max_rows = 200L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0L) return("<p class=\"muted\">None.</p>")
  shown <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  cells <- lapply(shown, qc_format_column)
  numeric_col <- vapply(shown, is.numeric, logical(1))
  header <- paste0(
    "<tr>",
    paste0("<th", ifelse(numeric_col, " class=\"num\"", ""), ">",
           qc_html_escape(names(shown)), "</th>", collapse = ""),
    "</tr>"
  )
  rows <- vapply(seq_len(nrow(shown)), function(i) {
    paste0(
      "<tr>",
      paste0("<td", ifelse(numeric_col, " class=\"num\"", ""), ">",
             vapply(cells, function(col) col[i], character(1)), "</td>",
             collapse = ""),
      "</tr>"
    )
  }, character(1))
  note <- if (nrow(df) > max_rows) {
    sprintf("<p class=\"muted\">Showing the first %s of %s rows.</p>",
            format(max_rows, big.mark = ","), format(nrow(df), big.mark = ","))
  } else {
    ""
  }
  paste0("<div class=\"tablewrap\"><table>", header,
         paste(rows, collapse = "\n"), "</table></div>", note)
}

qc_name_list <- function(x, max_items = 8L) {
  x <- as.character(x)
  if (length(x) <= max_items) return(paste(x, collapse = ", "))
  paste0(paste(x[seq_len(max_items)], collapse = ", "),
         sprintf(", and %d more", length(x) - max_items))
}

write_qc_report <- function(qc, file,
                            title = "bNMF run QC report",
                            concordance_warn = 0.8,
                            concordance_min_n = 20L,
                            mismatch_warn = 0.05) {
  # Write the QC report for one run. Returns (invisibly) the character vector
  # of flags shown under "Things to check".

  has <- function(name) !is.null(qc[[name]]) && length(qc[[name]]) > 0L
  has_rows <- function(name) has(name) && is.data.frame(qc[[name]])
  flags <- character(0)
  sections <- character(0)
  add_section <- function(heading, ...) {
    sections <<- c(sections, paste0("<h2>", qc_html_escape(heading), "</h2>"), ...)
  }
  para <- function(...) paste0("<p>", paste0(...), "</p>")

  # ---- Run settings ----------------------------------------------------------
  if (has("settings")) {
    settings <- qc$settings
    settings_df <- data.frame(
      setting = names(settings),
      value = vapply(settings, function(v) paste(format(v, scientific = NA), collapse = ", "),
                     character(1)),
      stringsAsFactors = FALSE
    )
    add_section("Run settings", qc_html_table(settings_df))
  }

  # ---- Discovery GWAS inputs -------------------------------------------------
  if (has_rows("input_qc")) {
    input_qc <- qc$input_qc
    assumed <- input_qc$GWAS[grepl("^assumed", input_qc$effect_allele_rule)]
    if (length(assumed) > 0L) {
      flags <- c(flags, sprintf(
        "%d discovery GWAS file(s) have no effect-allele column, so the effect allele was assumed to be ALT: %s. Confirm this against each study's documentation.",
        length(assumed), qc_name_list(assumed)))
    }
    add_section(
      "Discovery GWAS inputs",
      para("How each discovery GWAS was read. <b>n_invalid_id</b> and ",
           "<b>n_invalid_alleles</b> count rows dropped because they are not ",
           "biallelic autosomal SNPs with A/C/G/T alleles. <b>n_no_direction</b> ",
           "counts rows whose effect is missing or zero."),
      qc_html_table(input_qc)
    )
  }

  # ---- Agreement with the primary GWAS --------------------------------------
  if (has_rows("concordance")) {
    concordance <- qc$concordance
    low <- concordance[!is.na(concordance$concordance) &
                         concordance$n_testable >= concordance_min_n &
                         concordance$concordance < concordance_warn, , drop = FALSE]
    if (nrow(low) > 0L) {
      flags <- c(flags, sprintf(
        "Low risk-allele agreement with the primary GWAS for %s. Agreement near 0%% usually means the file's effect allele is the opposite of what was assumed; agreement near 50%% points to a genome-build or allele-labelling problem.",
        qc_name_list(sprintf("%s (%.0f%%)", low$GWAS, 100 * low$concordance))))
    }
    shown <- concordance
    shown$concordance <- ifelse(is.na(shown$concordance), NA_character_,
                                sprintf("%.1f%%", 100 * shown$concordance))
    add_section(
      "Agreement with the primary GWAS",
      para("For each discovery GWAS: how many of its variants are in the primary ",
           "GWAS, how many it lists in the opposite allele order (harmless, ",
           "handled automatically), and how often the two GWAS agree on the risk ",
           "allele among variants with a direction in both and a primary-GWAS P ",
           "below the cutoff. Agreement should be close to 100%."),
      qc_html_table(shown)
    )
  }

  # ---- Variants removed by the primary-GWAS filters -------------------------
  if (has_rows("dropped") || has_rows("drop_summary")) {
    parts <- para("Variants removed before clumping and pruning because they are ",
                  "absent from the primary GWAS, lack a direction of effect, exceed ",
                  "the primary-GWAS P cutoff, or have a different risk allele in the ",
                  "discovery and primary GWAS.")
    if (has_rows("drop_summary")) {
      parts <- c(parts, "<h3>Counts by discovery GWAS and reason</h3>",
                 qc_html_table(qc$drop_summary))
    } else if (has_rows("dropped") && "reason" %in% names(qc$dropped)) {
      counts <- as.data.frame(table(reason = qc$dropped$reason), stringsAsFactors = FALSE)
      names(counts)[2] <- "n"
      parts <- c(parts, "<h3>Counts by reason</h3>", qc_html_table(counts))
    }
    if (has_rows("dropped") && has("sentinel_pvalue_cutoff") &&
        "PVALUE_discovery" %in% names(qc$dropped)) {
      dropped <- qc$dropped
      sentinel_level <- dropped[!is.na(dropped$PVALUE_discovery) &
                                  dropped$PVALUE_discovery < qc$sentinel_pvalue_cutoff, , drop = FALSE]
      sentinel_level <- sentinel_level[order(sentinel_level$PVALUE_discovery), , drop = FALSE]
      keep_cols <- intersect(c("VAR_KEY", "GWAS", "Population", "PVALUE_discovery",
                               "PVALUE_main", "Risk_Allele_discovery",
                               "Risk_Allele_main", "reason"), names(sentinel_level))
      if (nrow(sentinel_level) > 0L) {
        flags <- c(flags, sprintf(
          "%s variant(s) that met the sentinel P cutoff in a discovery GWAS were removed by the primary-GWAS filters (%s). They are listed in the report.",
          format(nrow(sentinel_level), big.mark = ","),
          qc_name_list(sprintf("%s: %d", names(table(sentinel_level$reason)),
                               as.integer(table(sentinel_level$reason))))))
      }
      parts <- c(parts,
                 sprintf("<h3>Removed variants that met the sentinel cutoff (discovery P &lt; %s)</h3>",
                         format(qc$sentinel_pvalue_cutoff, scientific = TRUE)),
                 para("These would have been sentinel candidates. Each locus can ",
                      "still be represented by another variant that passed."),
                 qc_html_table(sentinel_level[, keep_cols, drop = FALSE], max_rows = 100L))
    }
    add_section("Variants removed by the primary-GWAS filters", parts)
  }

  # ---- Sentinel checks --------------------------------------------------------
  # Independent re-check of the final sentinel list (whatever its source),
  # so a bad list or a regression in the selection code cannot pass silently.
  if (has_rows("sentinels") && has("sentinel_pvalue_cutoff") && has("main_pvalue_cutoff")) {
    sen <- qc$sentinels
    n_sen <- nrow(sen)
    col_or_na <- function(name, na) if (name %in% names(sen)) sen[[name]] else rep(na, n_sen)
    p_disc <- suppressWarnings(as.numeric(col_or_na("PVALUE_discovery", NA_real_)))
    p_main <- suppressWarnings(as.numeric(col_or_na("PVALUE_main", NA_real_)))
    ra_disc <- as.character(col_or_na("Risk_Allele_discovery", NA_character_))
    ra_main <- as.character(col_or_na("Risk_Allele", NA_character_))
    # Missing values count as failures, never as passes.
    ok_disc <- !is.na(p_disc) & p_disc < qc$sentinel_pvalue_cutoff
    ok_main <- !is.na(p_main) & p_main <= qc$main_pvalue_cutoff
    ok_allele <- !is.na(ra_disc) & !is.na(ra_main) & ra_disc == ra_main
    checks <- data.frame(
      check = c(sprintf("Discovery P < %s in at least one discovery GWAS",
                        format(qc$sentinel_pvalue_cutoff, scientific = TRUE)),
                sprintf("Primary-GWAS P <= %s", format(qc$main_pvalue_cutoff)),
                "Same risk allele in the discovery and primary GWAS"),
      n_sentinels = n_sen,
      n_pass = c(sum(ok_disc), sum(ok_main), sum(ok_allele)),
      stringsAsFactors = FALSE
    )
    checks$n_fail <- checks$n_sentinels - checks$n_pass
    checks$result <- ifelse(checks$n_fail == 0L, "all pass", "FAIL")
    parts <- c(
      para(sprintf("Re-checked on the final list of %s sentinels (before any proxy replacement). ",
                   format(n_sen, big.mark = ",")),
           "The discovery P is the lowest P across the discovery GWAS. A missing value counts as a failure."),
      qc_html_table(checks)
    )
    failed <- !(ok_disc & ok_main & ok_allele)
    if (any(failed)) {
      id <- if ("VAR_ID" %in% names(sen)) as.character(sen$VAR_ID) else seq_len(n_sen)
      flags <- c(flags, sprintf(
        "%d sentinel(s) fail at least one of the sentinel checks (discovery P, primary-GWAS P, shared risk allele): %s. See the Sentinel checks section.",
        sum(failed), qc_name_list(id[failed])))
      bad <- data.frame(VAR_ID = id[failed],
                        PVALUE_discovery = p_disc[failed], PVALUE_main = p_main[failed],
                        Risk_Allele_discovery = ra_disc[failed], Risk_Allele_main = ra_main[failed],
                        discovery_P_ok = ok_disc[failed], main_P_ok = ok_main[failed],
                        allele_ok = ok_allele[failed], stringsAsFactors = FALSE)
      parts <- c(parts, "<h3>Sentinels failing a check</h3>", qc_html_table(bad))
    } else {
      parts <- c(parts, para("Every sentinel passes all three checks."))
    }
    add_section("Sentinel checks", parts)
  }

  # ---- Variant counts by stage ----------------------------------------------
  if (has_rows("funnel")) {
    add_section("Variant counts by stage", qc_html_table(qc$funnel))
  }

  # ---- Trait GWAS inputs -----------------------------------------------------
  if (has_rows("trait_qc")) {
    trait_qc <- qc$trait_qc
    if ("effect_allele_rule" %in% names(trait_qc)) {
      assumed_traits <- trait_qc$trait[grepl("^assumed", trait_qc$effect_allele_rule)]
      if (length(assumed_traits) > 0L) {
        flags <- c(flags, sprintf(
          "%d trait file(s) have no effect-allele column, so the effect allele was assumed to be ALT: %s. A wrong assumption flips that trait's z-scores and nothing downstream can detect it.",
          length(assumed_traits), qc_name_list(assumed_traits)))
      }
    }
    if ("n_complemented" %in% names(trait_qc)) {
      comp <- trait_qc[!is.na(trait_qc$n_complemented) & trait_qc$n_complemented > 0, , drop = FALSE]
      if (nrow(comp) > 0L) {
        flags <- c(flags, sprintf(
          "%d trait file(s) had variants matched by strand complement, i.e. reported on the opposite strand: %s.",
          nrow(comp), qc_name_list(sprintf("%s (%s)", comp$trait,
                                           format(comp$n_complemented, big.mark = ",", trim = TRUE)))))
      }
    }
    if (all(c("n_allele_mismatch", "n_position_matched") %in% names(trait_qc))) {
      share <- trait_qc$n_allele_mismatch / pmax(trait_qc$n_position_matched, 1)
      high <- trait_qc[!is.na(share) & trait_qc$n_position_matched >= concordance_min_n &
                         share > mismatch_warn, , drop = FALSE]
      if (nrow(high) > 0L) {
        flags <- c(flags, sprintf(
          "%d trait file(s) had more than %.0f%% of position-matched variants rejected for allele mismatch: %s. This can indicate a different genome build or mislabelled alleles.",
          nrow(high), 100 * mismatch_warn, qc_name_list(high$trait)))
      }
    }
    add_section(
      "Trait GWAS inputs",
      para("How each trait file was read at the requested positions. ",
           "<b>n_direct</b>: alleles matched as reported. <b>n_complemented</b>: ",
           "matched after complementing (opposite strand). <b>n_allele_mismatch</b>: ",
           "position found but alleles matched on neither strand, so no z-score."),
      qc_html_table(trait_qc, max_rows = 2000L)
    )
  }

  # ---- Trait filters ---------------------------------------------------------
  if (has_rows("trait_filters")) {
    trait_filters <- qc$trait_filters
    counts <- as.data.frame(table(result = trait_filters$result), stringsAsFactors = FALSE)
    names(counts)[2] <- "n_traits"
    removed <- trait_filters[trait_filters$result != "trait kept", , drop = FALSE]
    add_section(
      "Trait filters",
      qc_html_table(counts),
      "<h3>Traits removed</h3>",
      qc_html_table(removed, max_rows = 2000L)
    )
  }

  # ---- Sentinel trait coverage ----------------------------------------------
  if (has("sentinel_nonmissingness")) {
    coverage <- qc$sentinel_nonmissingness
    cutoff <- if (has("nonmissing_cutoff")) qc$nonmissing_cutoff else 0.8
    n_zero <- sum(!is.na(coverage) & coverage == 0)
    if (n_zero > 0L) {
      flags <- c(flags, sprintf(
        "%d sentinel(s) have no data in any trait file: %s.",
        n_zero, qc_name_list(names(coverage)[!is.na(coverage) & coverage == 0])))
    }
    bins <- data.frame(
      coverage = c("No trait data", "Above 0 and below 50%",
                   sprintf("50%% to below %.0f%%", 100 * cutoff),
                   sprintf("%.0f%% or more", 100 * cutoff)),
      n_sentinels = c(
        n_zero,
        sum(coverage > 0 & coverage < 0.5, na.rm = TRUE),
        sum(coverage >= 0.5 & coverage < cutoff, na.rm = TRUE),
        sum(coverage >= cutoff, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
    add_section(
      "Sentinel trait coverage",
      para(sprintf("Fraction of the retained traits with data, per sentinel. Sentinels below %.0f%% are sent to the proxy search.",
                   100 * cutoff)),
      qc_html_table(bins)
    )
  }

  # ---- Proxy search ----------------------------------------------------------
  if (has_rows("proxies_needed") || has_rows("proxies_found")) {
    needed <- if (has_rows("proxies_needed")) qc$proxies_needed else data.frame()
    found <- if (has_rows("proxies_found")) qc$proxies_found else data.frame()
    parts <- character(0)
    if (nrow(needed) > 0L && "reason" %in% names(needed)) {
      reason_counts <- as.data.frame(table(reason = needed$reason), stringsAsFactors = FALSE)
      names(reason_counts)[2] <- "n_sentinels"
      parts <- c(parts, "<h3>Sentinels needing a proxy, by reason</h3>",
                 qc_html_table(reason_counts))
    }
    if (has("proxy_unanswered")) {
      flags <- c(flags, sprintf(
        "LDlink returned no answer for %d flagged sentinel(s): %s. They got no proxy in this run; rerun to ask again.",
        length(qc$proxy_unanswered), qc_name_list(qc$proxy_unanswered)))
    }
    if (nrow(needed) > 0L && "VAR_ID" %in% names(needed)) {
      found_ids <- if ("original_SNP" %in% names(found)) found$original_SNP else character(0)
      not_found <- setdiff(needed$VAR_ID, found_ids)
      parts <- c(parts, para(sprintf(
        "%d sentinel(s) needed a proxy; a proxy was found for %d and not for %d. A flagged sentinel with no proxy is dropped.",
        nrow(needed), length(intersect(needed$VAR_ID, found_ids)), length(not_found))))
      if (length(not_found) > 0L) {
        flags <- c(flags, sprintf(
          "%d sentinel(s) needed a proxy and no proxy was found, so they were dropped: %s.",
          length(not_found), qc_name_list(not_found)))
        # Prefer the proxy stage's own account of why (select_proxies()$no_proxy).
        without <- if (has_rows("proxies_not_found")) {
          names(qc$proxies_not_found)[names(qc$proxies_not_found) == "reason"] <- "why_no_proxy"
          merge(needed[needed$VAR_ID %in% not_found, , drop = FALSE],
                qc$proxies_not_found, by = "VAR_ID", all.x = TRUE, sort = FALSE)
        } else {
          needed[needed$VAR_ID %in% not_found, , drop = FALSE]
        }
        if ("why_no_proxy" %in% names(without)) {
          why_counts <- as.data.frame(table(why_no_proxy = without$why_no_proxy),
                                      stringsAsFactors = FALSE)
          names(why_counts)[2] <- "n_sentinels"
          parts <- c(parts, "<h3>Why no proxy was found</h3>", qc_html_table(why_counts))
        }
        parts <- c(parts, "<h3>Sentinels with no proxy found</h3>",
                   qc_html_table(without, max_rows = 500L))
      }
    }
    if (has_rows("proxy_candidates") && "status" %in% names(qc$proxy_candidates)) {
      outcome_counts <- as.data.frame(table(outcome = qc$proxy_candidates$status),
                                      stringsAsFactors = FALSE)
      names(outcome_counts)[2] <- "n_candidates"
      outcome_counts <- outcome_counts[order(-outcome_counts$n_candidates), , drop = FALSE]
      parts <- c(parts, "<h3>Proxy candidates by outcome</h3>",
                 para("Every LD partner LDlink returned at or above the r2 threshold, and what became of it."),
                 qc_html_table(outcome_counts))
    }
    if (nrow(found) > 0L) {
      parts <- c(parts, "<h3>Proxies chosen</h3>", qc_html_table(found, max_rows = 500L))
    }
    add_section("Proxy search", parts)
  }

  # ---- Final variant set -----------------------------------------------------
  if (has_rows("final_variants")) {
    final_variants <- qc$final_variants
    parts <- para(sprintf("%s variants in the final set.",
                          format(nrow(final_variants), big.mark = ",")))
    if ("Variant_Type" %in% names(final_variants)) {
      type_counts <- as.data.frame(table(Variant_Type = final_variants$Variant_Type),
                                   stringsAsFactors = FALSE)
      names(type_counts)[2] <- "n_variants"
      parts <- c(parts, qc_html_table(type_counts))
    }
    # Per-variant p-values. For a proxy, PVALUE_discovery is the replaced
    # sentinel's discovery P and PVALUE_main is the proxy's own main-GWAS P.
    p_cols <- intersect(c("PVALUE_discovery", "PVALUE_main"), names(final_variants))
    if (length(p_cols) > 0L) {
      range_lines <- vapply(p_cols, function(col) {
        p <- suppressWarnings(as.numeric(final_variants[[col]]))
        if (all(is.na(p))) return(sprintf("%s: no values", col))
        sprintf("%s: min %s, median %s, max %s", col,
                format(min(p, na.rm = TRUE), digits = 3, scientific = TRUE),
                format(stats::median(p, na.rm = TRUE), digits = 3, scientific = TRUE),
                format(max(p, na.rm = TRUE), digits = 3, scientific = TRUE))
      }, character(1))
      parts <- c(parts, para(paste(qc_html_escape(range_lines), collapse = "<br>")))
      id_cols <- intersect(c("VAR_ID", "rsID", "Variant_Type", "original_SNP"),
                           names(final_variants))
      p_table <- final_variants[, c(id_cols, p_cols), drop = FALSE]
      if ("PVALUE_main" %in% p_cols) {
        p_table <- p_table[order(p_table$PVALUE_main), , drop = FALSE]
      }
      parts <- c(parts, qc_html_table(p_table, max_rows = nrow(p_table)))
    }
    add_section("Final variant set", parts)
  }

  # ---- Comparison with an earlier run ---------------------------------------
  if (has("comparison") && is.list(qc$comparison)) {
    cmp <- qc$comparison
    previous <- unique(as.character(cmp$previous))
    current <- unique(as.character(cmp$current))
    only_previous <- setdiff(previous, current)
    only_current <- setdiff(current, previous)
    cmp_table <- data.frame(
      set = c("In both runs", "Only in the earlier run", "Only in this run"),
      n_variants = c(length(intersect(previous, current)),
                     length(only_previous), length(only_current)),
      stringsAsFactors = FALSE
    )
    parts <- c(para("Final variant sets compared with: <b>",
                    qc_html_escape(if (!is.null(cmp$label)) cmp$label else "an earlier run"),
                    "</b>."),
               qc_html_table(cmp_table))
    if (length(only_previous) > 0L) {
      parts <- c(parts, "<h3>Only in the earlier run</h3>",
                 qc_html_table(data.frame(variant = only_previous), max_rows = 500L))
    }
    if (length(only_current) > 0L) {
      parts <- c(parts, "<h3>Only in this run</h3>",
                 qc_html_table(data.frame(variant = only_current), max_rows = 500L))
    }
    add_section("Comparison with an earlier run", parts)
  }

  # ---- bNMF ------------------------------------------------------------------
  if (has_rows("bnmf_run_summary") && "K" %in% names(qc$bnmf_run_summary)) {
    k_counts <- as.data.frame(table(K = qc$bnmf_run_summary$K), stringsAsFactors = FALSE)
    names(k_counts)[2] <- "n_repetitions"
    k_counts$share <- sprintf("%.0f%%", 100 * k_counts$n_repetitions / sum(k_counts$n_repetitions))
    add_section(
      "bNMF",
      para("Number of clusters (K) chosen across repetitions."),
      qc_html_table(k_counts)
    )
  }

  # ---- Notes -----------------------------------------------------------------
  if (has("notes")) {
    add_section("Notes", paste0("<ul>", paste0("<li>", qc_html_escape(qc$notes), "</li>",
                                               collapse = ""), "</ul>"))
  }

  # ---- Assemble --------------------------------------------------------------
  run_name <- if (has("run_name")) as.character(qc$run_name)[1] else ""
  check_block <- if (length(flags) > 0L) {
    paste0("<div class=\"checks\"><h2>Things to check</h2><ul>",
           paste0("<li>", qc_html_escape(flags), "</li>", collapse = ""),
           "</ul></div>")
  } else {
    "<div class=\"checks ok\"><h2>Things to check</h2><p>No flags were raised.</p></div>"
  }
  css <- paste(
    "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f6f7f9;color:#1c2330;}",
    "main{max-width:1100px;margin:0 auto;padding:24px 16px 64px;}",
    "h1{font-size:24px;margin:8px 0 4px;} h2{font-size:18px;margin:32px 0 8px;} h3{font-size:14px;margin:18px 0 6px;color:#3b4658;}",
    "p{line-height:1.5;margin:6px 0;} .muted{color:#6a7486;font-size:13px;}",
    ".meta{color:#6a7486;font-size:13px;margin-bottom:16px;}",
    ".checks{background:#fff7e6;border:1px solid #f0c36d;border-radius:8px;padding:4px 16px 12px;}",
    ".checks.ok{background:#eef8f0;border-color:#9fd3aa;} .checks h2{margin-top:12px;} .checks li{margin:6px 0;line-height:1.45;}",
    ".tablewrap{overflow-x:auto;background:#fff;border:1px solid #dfe3ea;border-radius:8px;}",
    "table{border-collapse:collapse;width:100%;font-size:13px;}",
    "th,td{padding:6px 10px;border-bottom:1px solid #eceff4;text-align:left;white-space:nowrap;}",
    "th{background:#f0f2f6;font-weight:600;} tr:last-child td{border-bottom:none;}",
    ".num{text-align:right;font-variant-numeric:tabular-nums;}",
    sep = "\n"
  )
  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    paste0("<title>", qc_html_escape(title), "</title>"),
    paste0("<style>", css, "</style></head><body><main>"),
    paste0("<h1>", qc_html_escape(title), "</h1>"),
    paste0("<div class=\"meta\">", qc_html_escape(run_name),
           if (nzchar(run_name)) " &middot; " else "",
           "generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), "</div>"),
    check_block,
    sections,
    "</main></body></html>"
  )

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, file, useBytes = TRUE)
  message(sprintf("QC report written to %s (%d item(s) to check).", file, length(flags)))
  invisible(flags)
}
