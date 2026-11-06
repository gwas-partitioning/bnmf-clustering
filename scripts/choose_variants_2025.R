# =============================================================================
# choose_variants_2025.R
# =============================================================================
# Helper functions for GWAS variant selection, LD clumping, and proxy search.
#
# Functions:
#   standardize_gwas_variants()  — Put one GWAS table into the internal variant
#                                  representation (variant key, effect allele, risk allele)
#   get_sig_snps()               — Collect variants below a P cutoff from each discovery GWAS
#   get_biggest_gwas()           — Cross-reference them against the primary GWAS and apply
#                                  the primary-GWAS filters (P and risk-allele agreement)
#   snp_clump()                  — Position-based clumping (no external API required)
#   ld_pruning_SNP.clip()        — LD pruning via LDlink SNPclip (multi-population)
#   prune()                      — Single-chromosome LD pruning via LDlinkR or curl
#   ld_pruning()                 — Full-genome LD pruning via LDlinkR (legacy)
#   ld_pruning_topLD_api()       — LD pruning via TopLD local binary (legacy)
#   count_traits_per_variant()   — (deprecated) Count non-missing traits per variant
#   count_traits_per_variant_2025() — Parallel version of above
#   window_to_sentinels()        — Filter variant pool to within kb-window of sentinel variants
#   find_variants_needing_proxies() — Flag strand-ambiguous, multiallelic, or high-missingness variants
#   choose_proxies()             — Search for LD proxies via LDlinkR::LDproxy_batch
#
# Assumptions:
#   - Genome build: hg19/GRCh37
#   - Biallelic autosomal SNPs only
#   - Each GWAS file locates variants by VAR_ID (CHR_POS_REF_ALT), by CHR + POS,
#     or by SNP (chr:pos), and gives BETA (or ODDS_RATIO) and P_VALUE
#   - Files do NOT need to share one REF/ALT order; see the allele helpers below
#
# P-value columns:
#   PVALUE_discovery — lowest P across the discovery GWAS. Decides whether a
#                      variant is eligible to be a sentinel.
#   PVALUE_main      — P in the primary GWAS. Must be <= main_pval_cutoff for a
#                      variant to stay in, and orders variants within a locus.
#   There is deliberately no bare PVALUE column after get_sig_snps().
# =============================================================================

library(tidyverse)
library(data.table)
library(LDlinkR)
library(future)
library(future.apply)

# -----------------------------------------------------------------------------
# Allele and variant-key helpers
# -----------------------------------------------------------------------------
# The pipeline does not require input files to share one REF/ALT order. Every
# variant is keyed internally by position plus its allele pair in a fixed
# (alphabetical) order, and the effect allele is tracked explicitly, so two
# files that list the same variant as A_G and G_A meet on one key.

# Column names accepted as the effect allele / other allele, in priority order.
# (prep_bNMF_2025.R keeps an identical copy so each script can be sourced alone.)
bnmf_effect_allele_names <- function() {
  c("EA", "Effect_Allele", "Effect_Allele_PH", "effect_allele", "A1")
}
bnmf_other_allele_names <- function() {
  c("NEA", "Other_Allele", "other_allele", "Non_Effect_Allele", "A2")
}

canonical_variant_key <- function(chr, pos, a1, a2) {
  # "CHR:POS:A:B" with the two alleles in alphabetical order.
  chr <- sub("^chr", "", as.character(chr), ignore.case = TRUE)
  pos <- as.integer(pos)
  a1 <- toupper(as.character(a1))
  a2 <- toupper(as.character(a2))
  paste(chr, pos, pmin(a1, a2), pmax(a1, a2), sep = ":")
}

standardize_gwas_variants <- function(df, effect_allele = NULL) {
  # Put one GWAS table into the pipeline's internal variant representation.
  #
  # Accepted layouts:
  #   - VAR_ID (CHR_POS_REF_ALT), optionally with an effect-allele column
  #   - CHR + POS (or SNP as chr:pos) with effect- and other-allele columns,
  #     or with REF + ALT columns
  # plus BETA (or ODDS_RATIO).
  #
  # Effect-allele rule:
  #   1. `effect_allele` override, when given: "ALT", "REF", or a column name
  #   2. else the first column found among bnmf_effect_allele_names()
  #   3. else ALT (the second allele of CHR_POS_REF_ALT) is assumed
  #
  # Only biallelic autosomal SNPs with A/C/G/T alleles are kept.
  #
  # Returns the retained rows with standardized columns (VAR_KEY, VAR_ID, CHR,
  # POS, REF, ALT, EA, NEA, BETA, Risk_Allele, Nonrisk_Allele) followed by the
  # file's remaining columns. attr(, "qc") records what was dropped and which
  # effect-allele rule was used.

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  cols <- names(df)
  n_in <- nrow(df)

  # --- effect-allele rule (decided from the header, so it is known even for
  #     an empty table) -------------------------------------------------------
  ea_col <- NULL
  ea_rule <- NULL
  has_override <- !is.null(effect_allele) && length(effect_allele) == 1L &&
    !is.na(effect_allele) && nzchar(trimws(as.character(effect_allele)))
  if (has_override) {
    effect_allele <- trimws(as.character(effect_allele))
    if (toupper(effect_allele) %in% c("ALT", "REF")) {
      ea_rule <- toupper(effect_allele)
      ea_label <- paste0(ea_rule, " (manifest override)")
    } else if (effect_allele %in% cols) {
      ea_col <- effect_allele
      ea_label <- paste0("column ", ea_col, " (manifest override)")
    } else {
      stop("effect_allele override '", effect_allele,
           "' is neither ALT, REF, nor a column in the file.")
    }
  } else {
    hit <- intersect(bnmf_effect_allele_names(), cols)
    if (length(hit) > 0L) {
      ea_col <- hit[1]
      ea_label <- paste0("column ", ea_col)
    } else {
      ea_rule <- "ALT"
      ea_label <- "assumed ALT (no effect-allele column)"
    }
  }

  # --- position and the file's own allele order -----------------------------
  has_varid <- "VAR_ID" %in% cols
  has_ref_alt <- all(c("REF", "ALT") %in% cols)
  nea_col <- intersect(bnmf_other_allele_names(), cols)
  nea_col <- if (length(nea_col) > 0L) nea_col[1] else NULL

  if (has_varid) {
    var_id <- toupper(as.character(df$VAR_ID))
    valid_id <- !is.na(var_id) & grepl("^[0-9]+_[0-9]+_[ACGT]_[ACGT]$", var_id)
    parts <- data.table::tstrsplit(var_id[valid_id], "_", fixed = TRUE)
    if (length(parts) < 4L) parts <- rep(list(character(0)), 4L)
    chr <- parts[[1]]; pos <- parts[[2]]; ref <- parts[[3]]; alt <- parts[[4]]
  } else {
    if (all(c("CHR", "POS") %in% cols)) {
      chr_all <- sub("^chr", "", as.character(df$CHR), ignore.case = TRUE)
      pos_all <- as.character(df$POS)
    } else if ("SNP" %in% cols) {
      snp_parts <- data.table::tstrsplit(
        sub("^chr", "", as.character(df$SNP), ignore.case = TRUE), ":", fixed = TRUE
      )
      if (length(snp_parts) < 2L) snp_parts <- list(rep(NA_character_, n_in),
                                                    rep(NA_character_, n_in))
      chr_all <- snp_parts[[1]]
      pos_all <- snp_parts[[2]]
    } else {
      stop("A GWAS file needs VAR_ID (CHR_POS_REF_ALT), or CHR and POS columns, ",
           "or a SNP column in chr:pos form. Columns found: ",
           paste(cols, collapse = ", "))
    }
    valid_id <- !is.na(chr_all) & !is.na(pos_all) &
      grepl("^[0-9]+$", chr_all) & grepl("^[0-9]+$", pos_all)
    chr <- chr_all[valid_id]
    pos <- pos_all[valid_id]
    if (has_ref_alt) {
      ref <- toupper(trimws(as.character(df$REF[valid_id])))
      alt <- toupper(trimws(as.character(df$ALT[valid_id])))
    } else {
      ref <- NULL
      alt <- NULL
    }
  }
  df <- df[valid_id, , drop = FALSE]
  n_invalid_id <- n_in - nrow(df)

  # --- effect and other allele ----------------------------------------------
  if (!is.null(ea_col)) {
    ea <- toupper(trimws(as.character(df[[ea_col]])))
  } else {
    if (is.null(ref)) {
      stop("The file has no effect-allele column and no REF/ALT (or VAR_ID) to ",
           "fall back on. Accepted effect-allele columns: ",
           paste(bnmf_effect_allele_names(), collapse = ", "), ".")
    }
    ea <- if (ea_rule == "ALT") alt else ref
  }

  if (!is.null(ref)) {
    nea <- ifelse(ea == alt, ref, ifelse(ea == ref, alt, NA_character_))
  } else if (!is.null(nea_col)) {
    nea <- toupper(trimws(as.character(df[[nea_col]])))
    ref <- nea
    alt <- ea
  } else {
    stop("The file has an effect-allele column but no other-allele column and ",
         "no REF/ALT (or VAR_ID) to derive it from. Accepted other-allele columns: ",
         paste(bnmf_other_allele_names(), collapse = ", "), ".")
  }

  valid_alleles <- !is.na(ea) & !is.na(nea) &
    grepl("^[ACGT]$", ea) & grepl("^[ACGT]$", nea) & ea != nea &
    !is.na(ref) & !is.na(alt) & grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt)
  n_invalid_alleles <- sum(!valid_alleles)
  df <- df[valid_alleles, , drop = FALSE]
  chr <- chr[valid_alleles]; pos <- pos[valid_alleles]
  ref <- ref[valid_alleles]; alt <- alt[valid_alleles]
  ea <- ea[valid_alleles];   nea <- nea[valid_alleles]

  # --- effect size and risk allele ------------------------------------------
  if ("BETA" %in% cols) {
    beta <- suppressWarnings(as.numeric(df$BETA))
  } else if ("ODDS_RATIO" %in% cols) {
    beta <- suppressWarnings(log(as.numeric(df$ODDS_RATIO)))
  } else {
    stop("A GWAS file needs a BETA or ODDS_RATIO column. Columns found: ",
         paste(cols, collapse = ", "))
  }
  beta[!is.finite(beta)] <- NA_real_
  # A missing or exactly-zero effect gives no direction; the risk allele is NA
  # and the variant is removed when risk alleles are compared.
  has_direction <- !is.na(beta) & beta != 0
  risk <- ifelse(has_direction, ifelse(beta > 0, ea, nea), NA_character_)
  nonrisk <- ifelse(has_direction, ifelse(beta > 0, nea, ea), NA_character_)

  chr <- as.integer(chr)
  pos <- as.integer(pos)
  standardized <- data.frame(
    VAR_KEY = canonical_variant_key(chr, pos, ea, nea),
    VAR_ID = paste(chr, pos, ref, alt, sep = "_"),
    CHR = chr, POS = pos, REF = ref, ALT = alt, EA = ea, NEA = nea,
    BETA = beta, Risk_Allele = risk, Nonrisk_Allele = nonrisk,
    stringsAsFactors = FALSE
  )
  if (nrow(standardized) == 0L) {
    standardized$VAR_KEY <- character(0)
    standardized$VAR_ID <- character(0)
  }
  other_cols <- setdiff(names(df), names(standardized))
  out <- cbind(standardized, df[, other_cols, drop = FALSE])
  rownames(out) <- NULL

  attr(out, "qc") <- list(
    n_rows = n_in,
    n_invalid_id = n_invalid_id,
    n_invalid_alleles = n_invalid_alleles,
    n_no_direction = sum(!has_direction),
    effect_allele_rule = ea_label
  )
  out
}

read_gwas_rows_at_positions <- function(filepath, positions, rename_cols = NULL) {
  # Read only the rows of a (possibly gzipped) summary-statistics file whose
  # chromosome:position is in `positions` ("1:12345"). Matching is exact and
  # does not depend on allele order. All columns are returned as character.

  header <- names(fread(filepath, nrows = 0L))
  header_std <- header
  if (!is.null(rename_cols) && length(rename_cols) > 0L) {
    hit <- match(header_std, unname(rename_cols))
    header_std[!is.na(hit)] <- names(rename_cols)[hit[!is.na(hit)]]
  }
  empty <- as.data.frame(
    matrix(character(0), nrow = 0L, ncol = length(header_std),
           dimnames = list(NULL, header_std)),
    stringsAsFactors = FALSE
  )
  positions <- unique(positions[!is.na(positions)])
  if (length(positions) == 0L) return(empty)

  is_gz <- endsWith(filepath, ".gz")
  con <- if (is_gz) gzfile(filepath) else file(filepath)
  first_line <- tryCatch(readLines(con, n = 1L), finally = close(con))
  tab_delimited <- grepl("\t", first_line, fixed = TRUE)

  strip_chr <- 'sub(/^[Cc][Hh][Rr]/, "", %s)'
  if ("VAR_ID" %in% header_std) {
    rule <- sprintf('{ split($%d, a, "_"); if ((a[1] ":" a[2]) in k) print }',
                    match("VAR_ID", header_std))
  } else if (all(c("CHR", "POS") %in% header_std)) {
    rule <- sprintf('{ ch = $%d; %s; if ((ch ":" $%d) in k) print }',
                    match("CHR", header_std), sprintf(strip_chr, "ch"),
                    match("POS", header_std))
  } else if ("SNP" %in% header_std) {
    rule <- sprintf('{ s = $%d; %s; split(s, a, ":"); if ((a[1] ":" a[2]) in k) print }',
                    match("SNP", header_std), sprintf(strip_chr, "s"))
  } else {
    stop("Cannot locate variants in ", filepath, ": it needs VAR_ID, or CHR and ",
         "POS, or a SNP (chr:pos) column. Columns found: ",
         paste(header, collapse = ", "))
  }

  key_file <- tempfile("positions_", fileext = ".tmp")
  on.exit(unlink(key_file), add = TRUE)
  writeLines(positions, key_file)

  awk_program <- paste("NR == FNR { k[$1] = 1; next }", rule)
  awk_call <- paste(
    "awk", if (tab_delimited) "-F'\\t'" else "", shQuote(awk_program), shQuote(key_file)
  )
  cmd <- if (is_gz) {
    sprintf("gzip -cd %s | %s -", shQuote(filepath), awk_call)
  } else {
    sprintf("%s %s", awk_call, shQuote(filepath))
  }

  rows <- tryCatch(
    fread(cmd = cmd, header = FALSE, col.names = header_std,
          sep = if (tab_delimited) "\t" else "auto",
          colClasses = "character", data.table = FALSE, showProgress = FALSE),
    error = function(e) empty
  )
  rows
}

get_sig_snps <- function(gwas, PVCUTOFF = 5e-8, rename_cols = NULL) {
  # Collect every variant with P <= PVCUTOFF from each discovery GWAS.
  #
  # gwas: manifest data frame with columns ID, population, full_path and,
  #       optionally, effect_allele (per-file override of the effect-allele
  #       rule: "ALT", "REF", or a column name; blank = detect).
  #
  # Returns one row per variant per population (the lowest-P record within a
  # population), keyed by VAR_KEY, with the discovery P in PVALUE_discovery and
  # the risk allele taken from the sign of the effect on the effect allele.
  # attr(, "input_qc") is a per-file table of what was read and dropped.

  missing_cols <- setdiff(c("ID", "population", "full_path"), names(gwas))
  if (length(missing_cols) > 0L) {
    stop("The GWAS manifest is missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  has_override_col <- "effect_allele" %in% names(gwas)

  output_cols <- c("VAR_KEY", "VAR_ID", "CHR", "POS", "REF", "ALT",
                   "PVALUE_discovery", "Risk_Allele", "Nonrisk_Allele",
                   "GWAS", "Population")
  vars_sig_list <- list()
  qc_list <- list()

  for (cur_pop in unique(gwas$population)) {
    message(sprintf("Pulling significant SNPs from %s GWAS...", cur_pop))
    gwas_pop <- gwas[gwas$population == cur_pop, , drop = FALSE]

    vars_pop_list <- lapply(seq_len(nrow(gwas_pop)), function(i) {
      cur_id <- gwas_pop$ID[i]
      cur_file <- gwas_pop$full_path[i]
      cur_override <- if (has_override_col) gwas_pop$effect_allele[i] else NULL
      message(sprintf("...Reading %s...", cur_id))

      vars <- fread(cur_file, stringsAsFactors = FALSE, data.table = FALSE) %>%
        dplyr::rename(any_of(rename_cols))
      if (!"P_VALUE" %in% names(vars)) {
        stop(cur_id, ": no P_VALUE column. Map your column with rename_cols, ",
             "e.g. rename_cols = c(P_VALUE = \"P\"). Columns found: ",
             paste(names(vars), collapse = ", "))
      }
      n_read <- nrow(vars)
      vars$P_VALUE <- suppressWarnings(as.numeric(vars$P_VALUE))
      vars <- vars[!is.na(vars$P_VALUE) & vars$P_VALUE <= PVCUTOFF, , drop = FALSE]
      n_pass <- nrow(vars)

      std <- standardize_gwas_variants(vars, effect_allele = cur_override)
      qc <- attr(std, "qc")

      out <- std %>%
        mutate(PVALUE_discovery = if_else(P_VALUE == 0, .Machine$double.xmin, P_VALUE),
               GWAS = cur_id,
               Population = cur_pop) %>%
        select(all_of(output_cols))

      qc_list[[cur_id]] <<- data.frame(
        GWAS = cur_id, Population = cur_pop,
        n_rows_read = n_read,
        n_pass_pvalue = n_pass,
        n_invalid_id = qc$n_invalid_id,
        n_invalid_alleles = qc$n_invalid_alleles,
        n_no_direction = qc$n_no_direction,
        n_kept = nrow(out),
        effect_allele_rule = qc$effect_allele_rule,
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "...%s: %i variants with P <= %.1e (%i dropped as not biallelic autosomal SNPs); effect allele: %s",
        cur_id, nrow(out), PVCUTOFF,
        qc$n_invalid_id + qc$n_invalid_alleles, qc$effect_allele_rule
      ))
      out
    })

    # Within a population keep the lowest-P record per variant
    vars_pop <- bind_rows(vars_pop_list) %>%
      arrange(VAR_KEY, PVALUE_discovery) %>%
      distinct(VAR_KEY, .keep_all = TRUE)

    message(sprintf("Final # SNPs for %s: %i", cur_pop, nrow(vars_pop)))
    vars_sig_list[[cur_pop]] <- vars_pop
  }

  vars_sig <- bind_rows(vars_sig_list)
  message(sprintf("Final total # SNPs: %i", nrow(vars_sig)))

  attr(vars_sig, "input_qc") <- bind_rows(qc_list)
  vars_sig
}


get_biggest_gwas <- function(main_ss_filepath, vars_sig,
                             main_pval_cutoff = 0.05,
                             rename_cols = NULL,
                             effect_allele = NULL,
                             concordance_warn = 0.8,
                             concordance_min_n = 20L) {
  # Cross-reference the discovery variants against the primary (largest) GWAS
  # and apply the primary-GWAS filters, BEFORE any clumping or LD pruning.
  #
  # A variant is kept when all of these hold:
  #   - it is present in the primary GWAS (matched on position + allele pair,
  #     whatever allele order either file uses)
  #   - both GWAS give it a direction (effect is not missing or zero)
  #   - its primary-GWAS P is <= main_pval_cutoff
  #   - the discovery and primary GWAS agree on the risk allele
  #
  # Returns one row per kept variant with two P columns:
  #   PVALUE_discovery : lowest P across the discovery GWAS (sentinel eligibility)
  #   PVALUE_main      : P in the primary GWAS
  # VAR_ID, REF and ALT follow the primary GWAS's allele order, and BETA_main is
  # the primary-GWAS effect of ALT. Risk_Allele / Nonrisk_Allele are the
  # primary's (equal to the discovery GWAS's, since mismatches are removed).
  #
  # Attributes:
  #   "dropped"      one row per removed variant, with the reason
  #   "drop_summary" removed-variant counts by discovery GWAS and reason
  #   "concordance"  per discovery GWAS: agreement with the primary GWAS on the
  #                  risk allele, and how many variants use the opposite allele
  #                  order to the primary

  if (!file.exists(main_ss_filepath)) {
    stop("The primary summary stats file does not exist at: ", main_ss_filepath)
  }
  needed <- c("VAR_KEY", "CHR", "POS", "PVALUE_discovery", "Risk_Allele",
              "GWAS", "Population")
  if (!all(needed %in% names(vars_sig))) {
    stop("vars_sig must be the output of get_sig_snps(); missing column(s): ",
         paste(setdiff(needed, names(vars_sig)), collapse = ", "))
  }
  if (!is.numeric(main_pval_cutoff) || length(main_pval_cutoff) != 1L ||
      is.na(main_pval_cutoff) || main_pval_cutoff <= 0 || main_pval_cutoff > 1) {
    stop("main_pval_cutoff must be one number in (0, 1].")
  }

  # ---- 1. Primary-GWAS rows at the candidate positions ----
  message("Reading candidate positions from the primary GWAS...")
  positions <- unique(paste(vars_sig$CHR, vars_sig$POS, sep = ":"))
  main_raw <- read_gwas_rows_at_positions(main_ss_filepath, positions, rename_cols)
  if (!"P_VALUE" %in% names(main_raw) && "PVALUE" %in% names(main_raw)) {
    names(main_raw)[names(main_raw) == "PVALUE"] <- "P_VALUE"
  }
  if (!"P_VALUE" %in% names(main_raw)) {
    stop("The primary GWAS has no P_VALUE column. Map your column with ",
         "rename_cols, e.g. rename_cols = c(P_VALUE = \"P\"). Columns found: ",
         paste(names(main_raw), collapse = ", "))
  }
  main_std <- standardize_gwas_variants(main_raw, effect_allele = effect_allele)
  message(sprintf("  Primary GWAS effect allele: %s",
                  attr(main_std, "qc")$effect_allele_rule))

  p_main <- suppressWarnings(as.numeric(main_std$P_VALUE))
  p_main[!is.na(p_main) & p_main == 0] <- .Machine$double.xmin
  if ("SE" %in% names(main_std)) {
    se_main <- suppressWarnings(as.numeric(main_std$SE))
  } else {
    message("  Primary GWAS has no SE column; deriving SE from BETA and P.")
    se_main <- abs(main_std$BETA / qnorm(p_main / 2))
  }

  main <- data.frame(
    VAR_KEY = main_std$VAR_KEY,
    VAR_ID = main_std$VAR_ID,
    REF = main_std$REF,
    ALT = main_std$ALT,
    # Reported as the effect of ALT, whichever allele the file calls the effect
    # allele, so BETA_main always pairs with the REF/ALT columns returned here.
    BETA_main = ifelse(main_std$EA == main_std$ALT, main_std$BETA, -main_std$BETA),
    SE_main = se_main,
    PVALUE_main = p_main,
    Risk_Allele_main = main_std$Risk_Allele,
    Nonrisk_Allele_main = main_std$Nonrisk_Allele,
    stringsAsFactors = FALSE
  )
  n_main_dup <- sum(duplicated(main$VAR_KEY))
  if (n_main_dup > 0L) {
    message(sprintf("  %i duplicated variant(s) in the primary GWAS; keeping the lowest-P row of each.",
                    n_main_dup))
    main <- main[order(main$PVALUE_main), , drop = FALSE]
    main <- main[!duplicated(main$VAR_KEY), , drop = FALSE]
  }

  # ---- 2. Join every discovery record to the primary ----
  records <- vars_sig %>%
    as.data.frame() %>%
    dplyr::rename(any_of(c(VAR_ID_discovery = "VAR_ID",
                           REF_discovery = "REF",
                           ALT_discovery = "ALT",
                           Risk_Allele_discovery = "Risk_Allele",
                           Nonrisk_Allele_discovery = "Nonrisk_Allele"))) %>%
    left_join(main, by = "VAR_KEY") %>%
    mutate(
      reason = case_when(
        is.na(VAR_ID)                  ~ "not_in_main_gwas",
        is.na(Risk_Allele_discovery)   ~ "no_discovery_direction",
        is.na(Risk_Allele_main)        ~ "no_main_direction",
        is.na(PVALUE_main)             ~ "main_pvalue_missing",
        PVALUE_main > main_pval_cutoff ~ "main_pvalue_above_cutoff",
        Risk_Allele_discovery != Risk_Allele_main ~ "risk_allele_mismatch",
        TRUE                           ~ NA_character_
      )
    )

  # ---- 3. Per-file agreement with the primary (all records, before dedup) ----
  concordance <- records %>%
    group_by(GWAS, Population) %>%
    summarise(
      n_variants = n(),
      n_in_main = sum(!is.na(VAR_ID)),
      n_flipped_order = if ("REF_discovery" %in% names(records)) {
        sum(!is.na(VAR_ID) & REF_discovery != REF)
      } else {
        NA_integer_
      },
      n_testable = sum(!is.na(Risk_Allele_discovery) & !is.na(Risk_Allele_main) &
                         !is.na(PVALUE_main) & PVALUE_main <= main_pval_cutoff),
      n_concordant = sum(!is.na(Risk_Allele_discovery) & !is.na(Risk_Allele_main) &
                           !is.na(PVALUE_main) & PVALUE_main <= main_pval_cutoff &
                           Risk_Allele_discovery == Risk_Allele_main),
      .groups = "drop"
    ) %>%
    mutate(concordance = if_else(n_testable > 0L, n_concordant / n_testable, NA_real_)) %>%
    as.data.frame()

  message("Risk-allele agreement with the primary GWAS, by discovery GWAS:")
  for (i in seq_len(nrow(concordance))) {
    message(sprintf(
      "  %-30s %5d variants | %5d in primary | %5d listed in the opposite allele order | agreement %s (%d of %d testable)",
      concordance$GWAS[i], concordance$n_variants[i], concordance$n_in_main[i],
      concordance$n_flipped_order[i],
      if (is.na(concordance$concordance[i])) "n/a" else sprintf("%.1f%%", 100 * concordance$concordance[i]),
      concordance$n_concordant[i], concordance$n_testable[i]
    ))
  }
  suspicious <- concordance[!is.na(concordance$concordance) &
                              concordance$n_testable >= concordance_min_n &
                              concordance$concordance < concordance_warn, , drop = FALSE]
  if (nrow(suspicious) > 0L) {
    warning(
      "Low risk-allele agreement with the primary GWAS for: ",
      paste(sprintf("%s (%.0f%%)", suspicious$GWAS, 100 * suspicious$concordance),
            collapse = ", "),
      ". Agreement near 0% usually means that file's effect allele is the ",
      "opposite of what was assumed; agreement near 50% suggests a build or ",
      "allele-labelling problem. Check the effect-allele rule reported for the file.",
      call. = FALSE
    )
  }

  # ---- 4. One record per variant: the strongest discovery record ----
  best <- records %>%
    arrange(PVALUE_discovery) %>%
    distinct(VAR_KEY, .keep_all = TRUE)

  dropped <- best %>%
    filter(!is.na(reason)) %>%
    transmute(VAR_KEY, CHR, POS,
              VAR_ID_discovery = if ("VAR_ID_discovery" %in% names(best)) VAR_ID_discovery else NA_character_,
              GWAS, Population, PVALUE_discovery, PVALUE_main,
              Risk_Allele_discovery, Risk_Allele_main, reason) %>%
    as.data.frame()

  drop_summary <- dropped %>%
    dplyr::count(GWAS, Population, reason, name = "n") %>%
    as.data.frame()

  kept <- best %>%
    filter(is.na(reason)) %>%
    transmute(
      VAR_KEY, VAR_ID,
      CHR = as.integer(CHR), POS = as.integer(POS), REF, ALT,
      ChrPos = paste(CHR, POS, sep = ":"),
      ChrPos_LDlink = paste0("chr", ChrPos),
      PVALUE_discovery, PVALUE_main, BETA_main, SE_main,
      Risk_Allele = Risk_Allele_main,
      Nonrisk_Allele = Nonrisk_Allele_main,
      GWAS, Population
    ) %>%
    as.data.frame()

  message(sprintf(
    "%i of %i discovery variants kept after the primary-GWAS filters (P <= %g and the same risk allele).",
    nrow(kept), nrow(best), main_pval_cutoff
  ))
  if (nrow(dropped) > 0L) {
    reason_counts <- table(dropped$reason)
    message("  Removed: ",
            paste(sprintf("%s = %d", names(reason_counts), as.integer(reason_counts)),
                  collapse = "; "))
  }

  attr(kept, "dropped") <- dropped
  attr(kept, "drop_summary") <- drop_summary
  attr(kept, "concordance") <- concordance
  kept
}


prune <- function(var_df, my_token, population="CEU",
                  r2=0.1, method="curl"){
  pruned_vars <- c()
  if (nrow(var_df) == 1) {
    pruned_vars <- c(pruned_vars, var_df$rsID)
  }
  else if (nrow(var_df) > 1){
    if (method=="curl") {
      print(sprintf("Generating LD matrix using LDlinkRest (%i variants)...", nrow(var_df)))
      my_snps <- paste(var_df$rsID, collapse="%0A")

      ld_mat <- fread(cmd=sprintf("curl -k -X GET 'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?snps=%s&pop=%s&r2d=r2&token=%s'",
                          my_snps,
                          population,
                          my_token)
      )

    } else if (method=="LDlinkR") {
      print("Generating LD matrix using LDlinkR...")
      ld_mat <- LDlinkR::LDmatrix(snps=var_df$rsID,
                                  pop=population,
                                  r2d="r2",
                                  token=my_token)  ## This should be replaced by each user's own token (retrieve at: https://ldlink.nci.nih.gov/?tab=apiaccess)
      if (!is.null(ld_mat$X.)) {
        cat("Error within LDlinkR! Chromosome not pruned!")
      }
    }
    ld_mat <- ld_mat %>% column_to_rownames('RS_number')
    ld_mat <- as.matrix(ld_mat)
    ld_mat <- ld_mat[rowSums(is.na(ld_mat)) != ncol(ld_mat),
                     colSums(is.na(ld_mat)) != nrow(ld_mat)]
    remaining_snps <- var_df$rsID
    while(length(remaining_snps) > 0) {
      if (remaining_snps[1] %in% rownames(ld_mat)) {
        # only add to pruned_vars if found in LD matrix
        pruned_vars <- c(pruned_vars, remaining_snps[1])
        
        # remove current SNP and any other SNPs in LD with it
        remaining_snps <- setdiff(
          remaining_snps,
          rownames(ld_mat)[ld_mat[, remaining_snps[1]] >= r2]
        )
      } 
      else {
        # if current SNP not found in LD matrix, move on to next SNP
        remaining_snps <- setdiff(remaining_snps, remaining_snps[1])
      }
    }
  }
  pruning_output <- list(pruned_vars, ld_mat)
  return(pruning_output)
}

ld_pruning <- function(gwas_variants, rsID_map_file, my_token,
                       method="curl", population="CEU", r2=0.1) {
  
  # Given a data frame of original GWAS variants (VAR_ID) and p-values (PVALUE), 
  # prune to a set of independent variants based on some LD threshold
  # Leverage the LDlinkR package to fetch LD relationships for a set of input SNPs
  
  write(gwas_variants$VAR_ID, "all_gwas_varid.tmp")
  
  print("Grepping for VAR_IDs in rsID map...")
  all_var_df <- fread(cmd=paste0("grep -wFf all_gwas_varid.tmp ", rsID_map_file),
                      header=F, col.names=c("VAR_ID", "rsID"),
                      data.table=F, stringsAsFactors=F) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove=F) %>%
    inner_join(gwas_variants, by="VAR_ID") %>%
    arrange(PVALUE)  # This ordering is important for the pruning steps below!
  
  print(head(all_var_df))
  
  bad_rsIDs <- subset(all_var_df, !startsWith(rsID,"rs"))

  print(paste("Num. SNPs mapped to rsID:",nrow(all_var_df)))
  
  # system("rm all_gwas_varid.tmp")
  
  pruned_vars <- c()
  ld_mats <- c()
  
  for (i in 1:22) {
    start=Sys.time()
    
    var_df_chr <- filter(all_var_df, CHR == i)
    n_snps <- nrow(var_df_chr)
    
    n_split <- ifelse(method=="curl", 250, 1000)
    mat_name <- as.character(i)
    
    if (between(n_snps, 1 ,n_split)) {  # LDlinkR has query size limit of 1000
      
      print(paste0("Pruning chromosome ", i, "..."))
      var_df <- var_df_chr
  
      prune_output_chr <- prune(var_df=var_df,
                         my_token=my_token,
                         population=population,
                         r2=r2,
                         method=method)
      
      # store pruned vars and LD matrix
      pruned_vars <- c(pruned_vars, prune_output_chr[[1]])
      ld_mats[[mat_name]] <- prune_output_chr[[2]]
      
      print(sprintf("Chr%i pruned from %i to %i variants",
                    i, n_snps, length(prune_output_chr[[1]])))
    }
    else if (n_snps > n_split) {
      var_df_list <- split(var_df_chr, (seq(nrow(var_df_chr))-1) %/% n_split)
      pruned_vars_list <- c()
      for (j in 1:length(var_df_list)){
        print(sprintf("Pruning subset %i for chromosome %i...", j, i))
        var_df <- var_df_list[[j]]
        
        prune_output_chr_j <- prune(var_df=var_df,
                                  my_token=my_token,
                                  population=population,
                                  r2=r2,
                                  method=method)
        
        # store temp pruned vars and LD matrix
        pruned_vars_list <- c(pruned_vars_list, prune_output_chr_j[[1]])
        ld_mats[[paste(mat_name, j, sep="_")]] <- prune_output_chr_j[[2]]
        
      }
      print(sprintf("Final pruning for chromosome %i...", i))
      var_df <- filter(var_df_chr, rsID %in% pruned_vars_list)
      
      prune_output_chr <- prune(var_df=var_df,
                               my_token=my_token,
                               population=population,
                               r2=r2,
                               method=method)
      
      # store final pruned vars and LD matrix
      pruned_vars <- c(pruned_vars, prune_output_chr[[1]])
      ld_mats[[mat_name]] <- prune_output_chr[[2]]
      
      print(sprintf("Chr%i pruned from %i to %i variants",
                    i, n_snps, length(prune_output_chr[[1]])))
    }
    else {
      print(sprintf("No SNPs on Chr%i!", i))
    }
    end=Sys.time()
    print(end-start)
  }
  print(paste0(length(pruned_vars), " VARIANTS REMAIN AFTER ALL PRUNING!"))
  
  # final outputs
  df_pruned <- filter(all_var_df, rsID %in% pruned_vars)
  return(list(df_pruned, ld_mats))
}

ld_pruning_topLD_api <-  function(gwas_variants, api_path, rsID_map_file, r2=0.1, population="EUR") {
  write(gwas_variants$VAR_ID, "all_gwas_varid.tmp")

  print("Grepping for VAR_IDs in rsID map...")
  all_var_df <- fread(cmd=paste0("grep -wFf all_gwas_varid.tmp ", rsID_map_file),
                      header=F, col.names=c("VAR_ID", "rsID"),
                      data.table=F, stringsAsFactors=F) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove=F) %>%
    inner_join(gwas_variants, by="VAR_ID") %>%
    arrange(PVALUE)

  system("touch outputLD.txt")
  system("touch outputInfo.txt")
  pruned_vars <- c()
  for (i in 1:22) {
    tmp <- all_var_df %>%
      filter(CHR==i)
    if (nrow(tmp)>1) {
      print(sprintf("Getting LD for CHR %i (%i variants)...", i, nrow(tmp)))
      rsID_pairs <- combn(tmp$rsID, 2, FUN=paste, collapse=',')
      cat(sprintf("%i rsID combinations...\n\n", length(rsID_pairs)))
      write(rsID_pairs, "to_prune_rsIDs.tmp")
      
      system(sprintf("%s -thres 0.0 -pop %s -maf 0.01 -inFile to_prune_rsIDs.tmp -outputLD outputLD_temp.txt -outputInfo outputInfo_temp.txt", api_path, population))
      system("awk 'FNR>1' outputLD_temp.txt >> outputLD.txt")
      system("awk 'FNR>1' outputInfo_temp.txt >> outputInfo.txt")

      ld_mat <- fread("outputLD.txt",
                      stringsAsFactors = F, data.table = F) %>%
        subset(rsID1 %like% "rs" & rsID1 != "rsID1")
      
      info_mat <- fread("outputInfo.txt",
                        stringsAsFactors = F, data.table = F) %>%
        subset(rsID %in% all_var_df$rsID) %>%
        subset(!duplicated(rsID))
      
      remaining_snps <- tmp %>%
        filter(rsID %in% info_mat$rsID) %>%
        arrange(PVALUE) %>%
        pull(rsID)
    
    pruned_vars_chr <- c()
    all_snps <- unique(c(ld_mat$rsID2, ld_mat$rsID2))
    while(length(remaining_snps) > 0) {
      if (remaining_snps[1] %in% all_snps) {
        # add to pruned_vars list
        pruned_vars_chr <- c(pruned_vars_chr, remaining_snps[1])
        # find rows where current SNP is present and R2>threshold
        tmp_mat <- ld_mat %>%
          subset(rsID1 %in% remaining_snps[1] | rsID2 %in% remaining_snps[1]) %>%
          subset(R2>r2)
        # remove any SNPs found in those rows from remaining_snps (this should also remove current SNP)
        remaining_snps <- setdiff(remaining_snps, unique(c(tmp_mat$rsID1, tmp_mat$rsID2, remaining_snps[1])))
      } else{
        # do  keep SNPs that are not found
        remaining_snps <- setdiff(remaining_snps, remaining_snps[1])
        }
      }
  } else {
    pruned_vars_chr <- tmp$rsID
  }
  print(sprintf("CHR %i pruned from %i to %i SNPs...",i, nrow(tmp), length(pruned_vars_chr)))
  pruned_vars <- c(pruned_vars, pruned_vars_chr)
  print(length(pruned_vars))
  }
    
  pruned_variants <- all_var_df %>%
    filter(rsID %in% pruned_vars)
  return(pruned_variants)
}

check_order_by <- function(df, order_by, caller) {
  # The variant kept at each locus is the first one under `order_by`, so the
  # caller must say which P-value column(s) decide it. There is no default.
  if (missing(order_by) || is.null(order_by) || length(order_by) == 0L) {
    stop(caller, "() needs order_by: the column(s) that decide which variant ",
         "leads each locus, most important first, e.g. ",
         "order_by = c(\"PVALUE_main\", \"PVALUE_discovery\").", call. = FALSE)
  }
  missing_cols <- setdiff(order_by, names(df))
  if (length(missing_cols) > 0L) {
    stop(caller, "(): order_by column(s) not found in the input: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(order_by)
}

snp_clump <- function(df_snps,
                      id="rsID",
                      window=500000,
                      chr=c(1:22),
                      pos_range=c(1,Inf),
                      order_by) {

  check_order_by(df_snps, if (missing(order_by)) NULL else order_by, "snp_clump")

  clumped_snps <- c()
  uniq_chr <- sort(as.integer(unique(df_snps$CHR)))

  for (i in uniq_chr) {
    tmp <- df_snps %>%
      filter(CHR==i) %>%
      arrange(across(all_of(order_by))) %>%
      mutate(POS = as.integer(POS)) %>%
      data.frame()
    
    if (i %in% chr) {
      print(sprintf("Clumping Chr. %i...",i))
      
      do_clump <- tmp %>%
        filter(between(POS, pos_range[1], pos_range[2]))
      dont_clump <- tmp %>%
        filter(!between(POS, pos_range[1], pos_range[2]))
      clumped_snps <- c(clumped_snps, dont_clump[,id])
      print(sprintf("Clumping %i variants, not clumping %i variants...",
                    nrow(do_clump), nrow(dont_clump)))
      
      remaining_snps <- do_clump[,id]
      
      j=0
      while (length(remaining_snps)>0){
        clumped_snps <- c(clumped_snps, remaining_snps[1])
        
        cur_pos <- do_clump$POS[do_clump[id]==remaining_snps[1]]
        
        close_snps <- do_clump[abs(do_clump$POS-cur_pos)<=window, id]

        remaining_snps <- setdiff(remaining_snps, close_snps)

        j=j+1
      }
    num_clumped <- j + nrow(dont_clump)
    } else {
      print(sprintf("No clumping for Chr %i!",i))
      num_clumped=nrow(tmp)
      clumped_snps <- c(clumped_snps, tmp[,id])
    }
    cat(sprintf("Chr %i clumped from %i to %i SNPs\n\n",i, nrow(tmp), num_clumped))
  }
  print(sprintf("No. SNPs after clumping: %i",length(clumped_snps)))
  return(clumped_snps)
}

ld_pruning_SNP.clip <- function(df_snps,
                                pop,
                                r2 = 0.1,
                                maf = 0.01,
                                chr = 1:22,
                                output_dir = "./",
                                token,
                                parallel = FALSE,
                                order_by) {

  # SNPclip keeps the first variant of each LD block in input order, so the
  # caller must name the column(s) that decide the order. Checked before any
  # request is sent to LDlink.
  check_order_by(df_snps, if (missing(order_by)) NULL else order_by,
                 "ld_pruning_SNP.clip")

  if (parallel) library(future.apply)

  #— 0) PREP: ensure VAR_ID, split into fields, add ChrPos & order
  if (!"VAR_ID" %in% colnames(df_snps)) {
    stop("df_snps must contain a VAR_ID column.")
  }

  var_id_parts <- data.table::tstrsplit(df_snps$VAR_ID, "_", fixed = TRUE)
  snp_clip_input <- df_snps %>%
    mutate(
      CHR = as.integer(var_id_parts[[1]]),
      POS = as.integer(var_id_parts[[2]]),
      ChrPos = paste0("chr", CHR, ":", POS)
    ) %>%
    arrange(across(all_of(order_by)))
  
  #— 1) Find which chr‐files are already done
  pattern <- sprintf("^snpClip_results_%s_chr([0-9]+)\\.txt$", pop)
  done_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  done_chr <- if (length(done_files) > 0) {
    as.integer(sub(pattern, "\\1", basename(done_files)))
  } else integer(0)
  
  #— 2) Subset to only pending chromosomes
  pending_chr <- setdiff(chr, done_chr)
  if (length(pending_chr) == 0) {
    message("All population–chromosome pairs already completed. Reading results from disk...")
    out_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
    return(bind_rows(lapply(out_files, data.table::fread)))
  }
  message("Will process chromosomes: ", paste(pending_chr, collapse = ", "))
  
  #— 3) Single‐chromosome handler
  process_chromosome <- function(i) {
    st <- Sys.time()
    out_file <- file.path(output_dir, sprintf("snpClip_results_%s_chr%i.txt", pop, i))
    if (file.exists(out_file)) {
      message(" Chr ", i, " already done; skipping.")
      return(TRUE)
    }
    message(" Processing Chr ", i, " …")
    cur_chr <- snp_clip_input %>% filter(CHR == i)
    if (nrow(cur_chr) == 0) {
      message("  No SNPs on Chr ", i, "; skipping.")
      return(TRUE)
    }
    
    clipped_res <- NULL
    
    if (nrow(cur_chr) == 1) {
      # Single SNP: use LDhap
      clipped_res <- tryCatch({
        LDhap(
          snps = cur_chr$ChrPos,
          pop = pop,
          token = token,
          genome_build = "grch37",
          table_type = "variant"
        ) %>%
          rename(Alleles = Allele_Frequency) %>%
          mutate(Details = "Variant kept.")
      }, error = function(e) {
        message("LDhap error: ", e$message)
        NULL
      })
      
    } else if (nrow(cur_chr) <= 5000) {
      # Moderate size: prune directly
      clipped_res <- tryCatch({
        LDlinkR::SNPclip(
          snps = cur_chr$ChrPos,
          pop = pop,
          r2_threshold = r2,
          maf_threshold = maf,
          token = token,
          file = FALSE,
          genome_build = "grch37"
        )
      }, error = function(e) {
        message("SNPclip error: ", e$message)
        NULL
      })
      
    } else {
      # Large: split into chunks
      message(" Chromosome ", i, " has >5000 SNPs; breaking into sections...")
      var_df_list <- split(cur_chr, (seq(nrow(cur_chr)) - 1) %/% 5000)
      kept_snps <- c()
      chunk_failed <- FALSE
      
      for (j in seq_along(var_df_list)) {
        message(sprintf(" Pruning subset %i for chromosome %i...", j, i))
        var_df <- var_df_list[[j]]
        cur_snps_j <- var_df %>% pull(ChrPos)
        
        clipped_res_split <- tryCatch({
          LDlinkR::SNPclip(
            snps = cur_snps_j,
            pop = pop,
            r2_threshold = r2,
            maf_threshold = maf,
            token = token,
            file = FALSE,
            genome_build = "grch37"
          )
        }, error = function(e) {
          message("Subset SNPclip error: ", e$message)
          NULL
        })
        
        if (!is.null(clipped_res_split)) {
          kept_snps_subset <- clipped_res_split %>%
            filter(Details == "Variant kept.") %>%
            pull(RS_Number)
          kept_snps <- c(kept_snps, kept_snps_subset)
          message(sprintf(" Subset %i pruned to %i SNPs...", j, length(kept_snps_subset)))
        } else {
          chunk_failed <- TRUE
        }
      }
      
      if (chunk_failed || length(kept_snps) == 0L) {
        message(" At least one SNPclip subset failed; discarding the partial chromosome result.")
        clipped_res <- NULL
      } else {
        # Final pruning pass over combined kept SNPs
        message(sprintf(" Performing final chromosomal pruning for %i SNPs...", length(kept_snps)))
        clipped_res <- tryCatch({
          LDlinkR::SNPclip(
            snps = kept_snps,
            pop = pop,
            r2_threshold = r2,
            maf_threshold = maf,
            token = token,
            file = FALSE,
            genome_build = "grch37"
          )
        }, error = function(e) {
          message("Final SNPclip error: ", e$message)
          NULL
        })
      }
    }
    
    #— Validate the API response before trusting it. LDlinkR can return a
    # malformed object (rate limit, transient API error, empty body) without
    # actually throwing an R error, which would otherwise crash the filter()/
    # fwrite() below with an obscure "object 'Details' not found" error and take
    # the whole multi-chromosome loop down. Treat anything unexpected as a
    # failed chromosome instead -- it's skipped (no output file written), so
    # the file.exists() cache-skip above will retry just this chromosome on
    # the next run.
    if (!is.null(clipped_res) &&
        (!is.data.frame(clipped_res) || !all(c("Details", "RS_Number") %in% names(clipped_res)))) {
      message("  Chr ", i, ": unexpected API response (missing Details/RS_Number column) -- treating as failed.")
      message("  Raw response:\n  ", paste(utils::capture.output(str(clipped_res)), collapse = "\n  "))
      clipped_res <- NULL
    }

    #— Write out if successful
    if (!is.null(clipped_res)) {
      data.table::fwrite(clipped_res, file = out_file, sep = "\t", quote = FALSE)
      msg <- sprintf("  Chr %i: %i → %i SNPs",
                     i, nrow(cur_chr), nrow(clipped_res %>% filter(Details == "Variant kept.")))
      message(msg)
    } else {
      message("  No results for Chr ", i)
    }
    message("  elapsed: ", round(difftime(Sys.time(), st, units = "secs"), 1), "s")
    return(!is.null(clipped_res))
  }
  
  #— 4) Run chromosomes (parallel or not)
  if (parallel) {
    completed_ok <- future.apply::future_lapply(pending_chr, process_chromosome)
  } else {
    completed_ok <- lapply(pending_chr, process_chromosome)
  }
  failed_chr <- pending_chr[!unlist(completed_ok, use.names = FALSE)]
  if (length(failed_chr) > 0L) {
    stop("LD pruning failed for chromosome(s) ", paste(failed_chr, collapse = ", "),
         ". No partial result will be returned; rerun to retry them.")
  }
  
  #— 5) Read all results
  final_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  all_res <- bind_rows(lapply(final_files, data.table::fread))
  message("LD pruning complete across populations & chromosomes.")
  return(all_res)
}
count_traits_per_variant <- function(gwas_variants, ss_files) {
  # DEPRECATED: Use count_traits_per_variant_2025() for parallel execution.
  #
  # Given a vector of variants and a named vector of summary statistics files
  # for traits to be clustered, output a data frame of sample sizes per variant-trait.
  
  print("Assessing variant missingness across traits...")
  write(gwas_variants, "all_snps_varids.tmp")
  
  rename_cols <- c(N_PH="N", N_PH="Neff")
  
  variant_df_list <- lapply(1:length(ss_files), function(i) {
    print(sprintf("...Reading %s...", names(ss_files)[i]))
    
    headers <- as.character(fread(ss_files[i], nrows=1,
                                  data.table=F, stringsAsFactors=F, header=F))
    if (endsWith(ss_files[i],".gz")) {
      df <- fread(cmd=sprintf("gzip -cd %s | fgrep -wf all_snps_varids.tmp ",ss_files[i]),
                  header=F,
                  col.names=headers,
                  data.table=F,
                  stringsAsFactors=F) %>%
        rename(any_of(rename_cols))
      
      
    } else {
      df <- fread(cmd=sprintf("fgrep -wf all_snps_varids.tmp %s ",ss_files[i]),
                  header=F,
                  col.names=headers,
                  data.table=F,
                  stringsAsFactors=F) %>%
        rename(any_of(rename_cols))
      
    }
    if (nrow(df)==0) {
      df <- data.frame(VAR_ID="No matches", N_PH=0)
    }
    df <- df %>%
      filter(VAR_ID %in% gwas_variants) %>%
      select(VAR_ID, N_PH)
    print(nrow(df))
    return(df)
  })
  
  # make dataframe of Ns
  df_N <- variant_df_list %>%
    setNames(names(ss_files)) %>%
    bind_rows(.id="trait") %>%
    select(trait, VAR_ID, N_PH) %>%
    pivot_wider(names_from="trait", values_from="N_PH") %>%
    data.frame()
}


count_traits_per_variant_2025 <- function(gwas_variants, ss_files) {
  # Set up a parallel plan: use the minimum of the number of files and 8 workers
  plan(multisession, workers = 4)
  
  message("Assessing variant missingness across traits...")
  
  # Write variant IDs to a temporary file (used by the external grep command)
  tmp_file <- "all_snps_varids.tmp"
  writeLines(gwas_variants, con = tmp_file)
  
  # Use future_lapply to process each summary stats file in parallel
  variant_df_list <- future_lapply(seq_along(ss_files), function(i) {
    trait_name <- names(ss_files)[i]
    # Read header from the file (assume header is on the first line)
    headers <- as.character(fread(ss_files[i], nrows = 1,
                                  data.table = FALSE, header = FALSE, stringsAsFactors = FALSE))
    
    # Build the system command:
    cmd <- if (endsWith(ss_files[i], ".gz")) {
      sprintf("gzip -cd %s | fgrep -wf %s", ss_files[i], tmp_file)
    } else {
      sprintf("fgrep -wf %s %s", tmp_file, ss_files[i])
    }
    
    # Run the command and capture output using fread
    df <- fread(cmd = cmd,
                header = FALSE,
                col.names = headers,
                data.table = FALSE,
                stringsAsFactors = FALSE)
    
    # If no matches are found, create a placeholder data frame
    if (nrow(df) == 0) {
      df <- data.frame(VAR_ID = "No matches", N_PH = 0, stringsAsFactors = FALSE)
    }
    
    # Filter to keep only variants in our gwas_variants vector
    df <- df %>%
      filter(VAR_ID %in% gwas_variants) %>%
      select(VAR_ID, N_PH)
    
    return(df)
  })
  
  # Combine results into one data frame
  df_N <- variant_df_list %>%
    setNames(names(ss_files)) %>%
    bind_rows(.id = "trait") %>%
    select(trait, VAR_ID, N_PH) %>%
    pivot_wider(names_from = "trait", values_from = "N_PH") %>%
    as.data.frame()
  
  file.remove(tmp_file)
  return(df_N)
}# Example usage:


# -----------------------------------------------------------------------------
#' Filter a variant pool to within a genomic window of sentinel variants
#'
#' Used to create a tractable proxy candidate pool for fetch_summary_stats():
#' instead of fetching all P < PVCUTOFF_PROXY variants genome-wide (which can
#' be 100k+ for polygenic traits), restrict to variants near LD-pruned sentinels.
#' Any proxy returned by LDlinkR will by definition be in this neighbourhood.
#'
#' @param candidates  data.frame with CHR (int/char) and POS (int) columns
#' @param sentinels   data.frame with CHR and POS columns (LD-pruned sentinel variants)
#' @param window_kb   numeric; half-window size in kilobases (default 500)
#' @return Subset of candidates (same class/columns) within window_kb of at
#'   least one sentinel. Rows are unique; order is preserved.
window_to_sentinels <- function(candidates, sentinels, window_kb = 500) {
  window_bp <- as.integer(round(window_kb * 1000))

  cand_dt <- copy(as.data.table(candidates))
  sent_dt <- as.data.table(sentinels)

  # Point intervals for candidates
  cand_dt[, `:=`(chr_key = as.character(CHR),
                 start   = as.integer(POS),
                 end     = as.integer(POS))]

  # Sentinel windows
  sent_windows <- sent_dt[, .(
    chr_key = as.character(CHR),
    start   = pmax(0L, as.integer(POS) - window_bp),
    end     = as.integer(POS) + window_bp
  )]

  setkey(sent_windows, chr_key, start, end)
  setkey(cand_dt,      chr_key, start, end)

  hits <- foverlaps(cand_dt, sent_windows, type = "any", nomatch = NULL)

  # Return unique rows from candidates (VAR_ID preferred; fall back to row index)
  if ("VAR_ID" %in% names(candidates)) {
    out <- candidates[candidates$VAR_ID %in% unique(hits$VAR_ID), ]
  } else {
    out <- candidates[candidates$ChrPos %in% unique(hits$ChrPos), ]
  }

  message(sprintf(
    "  window_to_sentinels: %d / %d candidates retained within %g kb of a sentinel",
    nrow(out), nrow(candidates), window_kb
  ))
  out
}


# -----------------------------------------------------------------------------
find_variants_needing_proxies <- function(gwas_variant_df,
                                          var_nonmissingness,
                                          missing_cutoff = 0.8) {

  print("Choosing variants in need of proxies...")

  # Track reasons for needing proxies
  proxy_reasons_list <- list()

  # Find strand-ambiguous variants
  strand_ambig <- with(gwas_variant_df, VAR_ID[paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")])
  print(paste0("...", length(strand_ambig), " strand-ambiguous variants"))
  if (length(strand_ambig) > 0) {
    proxy_reasons_list[["ambiguous"]] <- data.frame(VAR_ID = strand_ambig, reason = "ambiguous", stringsAsFactors = FALSE)
  }

  # Find multi-allelic variants
  multi_allelic <- grep("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", gwas_variant_df$VAR_ID, value = TRUE)
  print(paste0("...", length(multi_allelic), " multi-allelic variants"))
  if (length(multi_allelic) > 0) {
    proxy_reasons_list[["multiallelic"]] <- data.frame(VAR_ID = multi_allelic, reason = "multiallelic", stringsAsFactors = FALSE)
  }

  # Every sentinel must have a non-missingness value. fetch_summary_stats()
  # returns a z-matrix row for every variant it was given (all NA when no trait
  # has data), so a sentinel with no value here was lost upstream. Treating that
  # as "high missingness" would send a rejected variant to the proxy search.
  no_value <- setdiff(gwas_variant_df$VAR_ID, names(var_nonmissingness))
  no_value <- union(no_value,
                    names(var_nonmissingness)[is.na(var_nonmissingness) &
                                                names(var_nonmissingness) %in% gwas_variant_df$VAR_ID])
  if (length(no_value) > 0L) {
    stop(sprintf(
      paste0("%d sentinel variant(s) have no non-missingness value (first: %s). ",
             "Every sentinel should have a row in the z-matrix; a missing one was ",
             "dropped between variant selection and the trait fetch and must be ",
             "investigated, not sent to the proxy search."),
      length(no_value), paste(head(no_value, 5), collapse = ", ")
    ), call. = FALSE)
  }

  # Find variants with excessive missingness (non-missing fraction below cutoff)
  low_cnt <- with(gwas_variant_df, VAR_ID[var_nonmissingness[VAR_ID] < missing_cutoff])
  print(paste0("...", length(low_cnt), " variants with excessive missingness"))
  if (length(low_cnt) > 0) {
    proxy_reasons_list[["high_missingness"]] <- data.frame(VAR_ID = low_cnt, reason = "high_missingness", stringsAsFactors = FALSE)
  }

  # Combine all reasons and handle overlaps
  if (length(proxy_reasons_list) > 0) {
    all_reasons <- bind_rows(proxy_reasons_list)
    
    # For variants with multiple reasons, combine them
    reasons_summary <- all_reasons %>%
      group_by(VAR_ID) %>%
      summarise(reason = paste(reason, collapse = ","), .groups = 'drop')
    
    need_proxies_varid <- unique(all_reasons$VAR_ID)
  } else {
    need_proxies_varid <- character(0)
    reasons_summary <- data.frame(VAR_ID = character(0), reason = character(0))
  }
  
  print(paste0("...", length(need_proxies_varid), " unique variants in total"))
  
  # Simply return the variants with their reasons and existing rsIDs
  result <- gwas_variant_df %>%
    filter(VAR_ID %in% need_proxies_varid) %>%
    left_join(reasons_summary, by = "VAR_ID") %>%
    select(VAR_ID, rsID = RS_Number, reason,
           any_of(c("PVALUE_main", "PVALUE_discovery", "PVALUE")))
  
  return(result)
}

choose_proxies <- function(need_proxies,
                           rsid_map_dir = "rsid_maps_by_chr",
                           pruned_variants,
                           zmat_fullset,           # matrix/data frame with rownames in "CHR:POS" format
                           token,                  # LDlink API token (Sys.getenv("LDLINK_TOKEN"))
                           population = "EUR",
                           frac_nonmissing_num = 0.8,
                           r2_num = 0.8,
                           variant_metadata = NULL,
                           output_dir = ".") {

  message(sprintf("Number of rows in need_proxies: %d", nrow(need_proxies)))

  # Initialize candidate storage
  candidates_with_data <- NULL

  # ----- Proxy search via LDlinkR::LDproxy_batch -----
  message(sprintf("Using LDlinkR to find proxies for %s population...", population))

  need_proxies <- need_proxies %>%
    separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
    mutate(query_snp = paste0("chr", CHR, ":", POS)) %>%
    select(-c(CHR, POS, REF, ALT))
  need_proxies_snps <- need_proxies$query_snp

  # Keep LDlink batch output and the failed-proxy log inside the run-specific
  # results directory so one analysis cannot accidentally reuse another's cache.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  proxy_out_file <- file.path(output_dir, "combined_query_snp_list_grch37.txt")
  proxy_meta_file <- file.path(output_dir, "combined_query_snp_list_grch37.meta.rds")
  expected_proxy_meta <- list(
    query_snps = sort(unique(need_proxies_snps)),
    population = population,
    genome_build = "grch37"
  )
  if (!file.exists(proxy_out_file)) {
    calling_dir <- getwd()
    on.exit(setwd(calling_dir), add = TRUE)
    setwd(output_dir)
    LDlinkR::LDproxy_batch(need_proxies_snps,
                           pop            = population,
                           r2d            = "r2",
                           token          = token,
                           append         = TRUE,
                           genome_build   = "grch37")
    if (file.exists(proxy_out_file)) saveRDS(expected_proxy_meta, proxy_meta_file)
  } else {
    if (!file.exists(proxy_meta_file) ||
        !identical(readRDS(proxy_meta_file), expected_proxy_meta)) {
      stop("Existing LDproxy cache does not match this query/population: ",
           proxy_out_file, ". Move or remove that cache before rerunning.")
    }
    message("Skipping LDproxy_batch — validated matching cache already exists.")
  }
  if (file.exists(proxy_out_file)) {
    proxy_df <- read.table(proxy_out_file, sep = "\t", header = TRUE, row.names = NULL) %>%
      filter(R2 > 0.1) %>%
      filter(!Coord %in% need_proxies_snps) %>%
      inner_join(need_proxies, by = "query_snp") %>%   # join on query_snp col present in LDproxy batch output
      distinct(rsID, RS_Number, .keep_all = TRUE) %>%
      select(rsID, proxy_rsID = RS_Number, r2 = R2)
  } else {
    proxy_df <- data.frame(rsID = character(), proxy_rsID = character(), r2 = numeric())
  }

  message(sprintf("Found %i potential proxies for the %i queried SNPs.", nrow(proxy_df), nrow(need_proxies)))
  need_proxies <- need_proxies %>% select(-query_snp)
  
  message(sprintf("Number of possible proxies found: %d", nrow(proxy_df)))
  
  # ----- Use chromosome-specific rsID maps to find proxy information -----
  
  if (nrow(proxy_df) > 0) {
    message("Creating proxy rsID map using chromosome-specific files...")
    
    proxy_rsids <- unique(proxy_df$proxy_rsID)
    potential_proxies_list <- list()
    chr_files <- list.files(rsid_map_dir, pattern = "chr.*\\.txt$", full.names = TRUE)
    
    if (length(chr_files) == 0) stop(sprintf("No chromosome files found in %s", rsid_map_dir))
    
    for (chr_file in chr_files) {
      chr_name <- basename(chr_file)
      # message(sprintf("Searching %s...", chr_name))
      chr_rsid_map <- fread(chr_file,
                            header = FALSE,
                            col.names = c("hg19_posID", "proxy_rsID", "ref_allele", "alt_allele"),
                            key = "proxy_rsID")
      chr_matches <- chr_rsid_map[proxy_rsID %in% proxy_rsids]
      if (nrow(chr_matches) > 0) potential_proxies_list[[chr_name]] <- chr_matches
    }
    
    # Combine results from all chromosomes
    if (length(potential_proxies_list) > 0) {
      potential_proxies_map <- rbindlist(potential_proxies_list)
      potential_proxies_map <- as.data.frame(potential_proxies_map)

      # Process Map: Extract Coords and Fix Orientation
      potential_proxies_map <- potential_proxies_map %>%
        separate(hg19_posID, into = c("CHR", "POS"), sep = ":", remove = FALSE) %>%
        mutate(CHR = gsub("chr", "", CHR),
               proxy_ChrPos = paste(CHR, POS, sep = ":")) %>%
        filter(proxy_ChrPos %in% rownames(zmat_fullset))

      if (!is.null(variant_metadata)) {
        metadata <- as.data.frame(variant_metadata)
        if (!"ChrPos" %in% names(metadata) && "SNP" %in% names(metadata)) {
          metadata <- metadata %>% rename(ChrPos = SNP)
        }
        required_metadata <- c("ChrPos", "VAR_ID", "REF", "ALT")
        if (!all(required_metadata %in% names(metadata))) {
          stop("variant_metadata must contain ChrPos (or SNP), VAR_ID, REF, and ALT.")
        }
        metadata <- metadata %>%
          select(proxy_ChrPos = ChrPos, proxy_VAR_ID = VAR_ID, REF, ALT) %>%
          distinct(proxy_ChrPos, .keep_all = TRUE)
        potential_proxies_map <- potential_proxies_map %>%
          select(-ref_allele, -alt_allele) %>%
          inner_join(metadata, by = "proxy_ChrPos")
      } else {
        warning("variant_metadata was not supplied; proxy VAR_ID allele orientation is taken from the rsID map.")
        potential_proxies_map <- potential_proxies_map %>%
          mutate(
            proxy_VAR_ID = paste(CHR, POS, ref_allele, alt_allele, sep = "_"),
            REF = ref_allele,
            ALT = alt_allele
          ) %>%
          select(-ref_allele, -alt_allele)
      }

      print(sprintf("%i of %i potential proxies in the full z-matrix...", nrow(potential_proxies_map), nrow(proxy_df)))

    } else {
      potential_proxies_map <- data.frame()
    }

    if (nrow(potential_proxies_map) > 0) {
      # ----- Assess missingness -----
      z_mat <- as.matrix(zmat_fullset)
      proxy_missingness <- rowSums(!is.na(z_mat)) / ncol(z_mat)
      proxy_missingness_df <- data.frame(
        proxy_ChrPos = names(proxy_missingness),
        frac_nonmissing = proxy_missingness,
        stringsAsFactors = FALSE
      )
      
      # ----- SAVE CANDIDATES WITH DATA (Before Filtering) -----
      candidates_with_data <- proxy_df %>%
        inner_join(potential_proxies_map, by = "proxy_rsID") %>%
        inner_join(proxy_missingness_df, by = "proxy_ChrPos")
      
      # ----- Apply Strict Filters for "Success" -----
      final_proxy_df <- candidates_with_data %>%
        filter(r2 >= r2_num) %>%
        filter(!(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")) &
                 !grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", proxy_VAR_ID)) %>%
        filter(frac_nonmissing >= frac_nonmissing_num)
      
      # Select best proxy
      final_proxy_df <- final_proxy_df %>% 
        group_by(rsID) %>%
        arrange(desc(frac_nonmissing), desc(r2), CHR) %>% 
        dplyr::slice(1) %>%
        ungroup() %>%
        inner_join(need_proxies, by = "rsID") 
      
    } else {
      final_proxy_df <- NULL
    }
  } else {
    final_proxy_df <- NULL
  }
  
  # Calculate Stats
  proxies_found <- if (!is.null(final_proxy_df)) final_proxy_df$rsID else character(0)
  no_proxies_found <- setdiff(need_proxies$rsID, proxies_found)
  
  message(sprintf("No proxies needed for %d variants.", length(setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID))))
  message(sprintf("Proxies found for %d variants.", length(proxies_found)))
  message(sprintf("No adequate proxies found for %d variants.", length(no_proxies_found)))
  
  if (length(no_proxies_found) > 0) {
    no_proxy_file <- file.path(output_dir, "no_proxies_found.txt")
    write(no_proxies_found, no_proxy_file)
    message("See ", no_proxy_file, " for a list of these variants.")
  }

  # ==========================================================
  # FINAL SELECTION LOGIC: FALLBACK WITH GUARDRAILS (FIXED)
  # ==========================================================
  
  # 1. Identify who got a proxy
  if (!is.null(final_proxy_df) && nrow(final_proxy_df) > 0) {
    replaced_originals <- unique(final_proxy_df$VAR_ID)
  } else {
    replaced_originals <- character(0)
  }
  
  # 2. Identify who FAILED to get a proxy
  failed_search_vars <- setdiff(need_proxies$VAR_ID, replaced_originals)
  
  # 3. Assess "Rescue" Potential for Failed Variants
  rescued_vars <- character(0)
  
  if (length(failed_search_vars) > 0) {
    message("Assessing failed variants for rescue (Strict QC)...")
    
    check_df <- data.frame(VAR_ID = failed_search_vars, stringsAsFactors = FALSE) %>%
      separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
      mutate(
        ChrPos_prefix = paste0("chr", CHR, ":", POS),  # Format: chr1:100
        ChrPos_noprefix = paste0(CHR, ":", POS)        # Format: 1:100
      )
    
    # Calculate missingness using BOTH formats to be safe
    if (!is.null(zmat_fullset)) {
      z_mat <- as.matrix(zmat_fullset)
      
      # Check which format matches the Z-matrix keys
      keys_prefix <- intersect(check_df$ChrPos_prefix, rownames(z_mat))
      keys_noprefix <- intersect(check_df$ChrPos_noprefix, rownames(z_mat))
      
      # Decide which key column to use
      if (length(keys_prefix) >= length(keys_noprefix)) {
        use_col <- "ChrPos_prefix"
        valid_keys <- keys_prefix
      } else {
        use_col <- "ChrPos_noprefix"
        valid_keys <- keys_noprefix
      }
      # message(sprintf("Matching Z-matrix using format: %s (Matches: %d)", use_col, length(valid_keys)))
      
      # Calculate fractions
      miss_stats <- rowSums(!is.na(z_mat[valid_keys, , drop=FALSE])) / ncol(z_mat)
      
      check_df$frac_nonmissing <- 0 
      match_idx <- match(check_df[[use_col]], names(miss_stats))
      check_df$frac_nonmissing[!is.na(match_idx)] <- miss_stats[na.omit(match_idx)]
      
    } else {
      warning("zmat_fullset missing; cannot calculate missingness for rescue. Defaulting to 0.")
      check_df$frac_nonmissing <- 0
    }
    
    # APPLY THE FILTERS
    rescue_candidates <- check_df %>%
      mutate(
        is_ambiguous = paste0(REF, ALT) %in% c("AT", "TA", "GC", "CG"),
        is_multiallelic = grepl(",", ALT),
        is_missing_too_high = frac_nonmissing < 0.5
      ) 
    
    # Corrected Logging: Use 'rescue_candidates' (the mutated df), not 'check_df'
    message(sprintf("  - Dropped %d ambiguous variants", sum(rescue_candidates$is_ambiguous)))
    message(sprintf("  - Dropped %d multi-allelic variants", sum(rescue_candidates$is_multiallelic)))
    message(sprintf("  - Dropped %d variants with <50%% data", sum(rescue_candidates$is_missing_too_high)))
    
    # Final Selection
    rescued_final <- rescue_candidates %>%
      filter(!is_ambiguous & !is_multiallelic & !is_missing_too_high)
    
    rescued_vars <- rescued_final$VAR_ID
    message(sprintf("  - RESCUED %d variants (Clean & >50%% data)", length(rescued_vars)))
  }
  
  # 4. Construct Final Sets
  vars_never_needed_proxy <- setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID)
  final_originals_to_keep <- c(vars_never_needed_proxy, rescued_vars)
  
  final_variant_set <- list(
    final_originals_to_keep,  # 1. Originals
    final_proxy_df           # 2. Proxies
  )
  
  return(final_variant_set)
}

suzuki_pruning <- function(assoc_data_filtered, my_LDlinkR_token, populations = c("EUR"), chunk_size = 2500, output_file = "ld_pruning_results.rds") {
  
  # Check if a previous run exists
  if (file.exists(output_file)) {
    message(sprintf("Resuming from previous run: %s", output_file))
    results <- readRDS(output_file)
    
    chr_NAs <- results %>% group_by(CHR) %>%
      summarize(n = sum(!is.na(clump))) %>%
      arrange(desc(n)) %>%
      dplyr::rename(non_NA = n)
    chr_counts <- results %>% dplyr::count(CHR) %>% dplyr::rename(snps=n)
    chr_checks <- chr_NAs %>% inner_join(chr_counts, by='CHR')
    processed_chr <- chr_checks %>%
      filter(snps==non_NA) %>%
      pull(CHR)
    
    assoc_data_sorted <- results
  } else {
    message("Starting new LD pruning run.")
    assoc_data_sorted <- assoc_data_filtered %>%
      arrange(P, CHR, BP) %>%
      mutate(clump = NA_integer_, index_snp = FALSE)
    processed_chr <- character(0) # Initialize empty character vector
  }
  
  clump_index <- max(assoc_data_sorted$clump, na.rm = TRUE) + 1
  if (!is.finite(clump_index)) clump_index <- 1 #set clump index to 1 if no clump has been assigned yet.
  unique_chr <- unique(assoc_data_sorted$CHR)
  
  for (chr in unique_chr) {
    if (chr %in% processed_chr) {
      message(sprintf("Chromosome %s already processed. Skipping.", chr))
      next
    }
    
    message(sprintf("Processing chromosome %s", chr))
    
    assoc_chr <- assoc_data_sorted %>% filter(CHR == chr)
    
    # Important: We need to reset these vectors based on the current state of assoc_data_sorted
    # This ensures we're working with the latest data
    snp_list <- assoc_chr$rsID
    bp_list <- assoc_chr$BP
    p_list <- assoc_chr$P
    
    if (length(snp_list) < 2) {
      message(sprintf("  Skipping chromosome %s due to less than 2 SNPs.", chr))
      processed_chr <- c(processed_chr, chr)
      saveRDS(assoc_data_sorted, file = output_file) # Save progress
      next
    }
    
    snp_chunks <- split(snp_list, ceiling(seq_along(snp_list) / chunk_size))
    bp_chunks <- split(bp_list, ceiling(seq_along(bp_list) / chunk_size))
    
    combined_ld_chunks <- list()
    
    for (chunk_idx in seq_along(snp_chunks)) {
      chunk_snps <- snp_chunks[[chunk_idx]]
      message(sprintf("  Processing chunk %d of %d", chunk_idx, length(snp_chunks)))
      
      ld_list <- list()
      for (pop in populations) {
        message(sprintf("    Calling LDmatrix for population %s", pop))
        ld_mat <- tryCatch({
          LDlinkR::LDmatrix(snps = chunk_snps,
                            pop = pop,
                            token = my_LDlinkR_token,
                            genome_build = "grch37")
        }, error = function(e) {
          message(sprintf("      Error for pop %s: %s", pop, e$message))
          return(NULL)
        })
        if (!is.null(ld_mat)) {
          mat <- as.matrix(ld_mat[, -1])
          rownames(mat) <- ld_mat$RS_number
          ld_list[[pop]] <- mat
        }
      }
      
      combined_ld_chunk <- Reduce(function(x, y) pmax(x, y, na.rm = TRUE), ld_list)
      combined_ld_chunks[[chunk_idx]] <- combined_ld_chunk
    }
    
    if (length(combined_ld_chunks) > 1) {
      all_snps_in_chr = unique(unlist(lapply(combined_ld_chunks, rownames)))
      combined_ld = matrix(NA, nrow=length(all_snps_in_chr), ncol=length(all_snps_in_chr), dimnames = list(all_snps_in_chr, all_snps_in_chr))
      
      for (ld_chunk in combined_ld_chunks){
        combined_ld[rownames(ld_chunk), colnames(ld_chunk)] = ld_chunk
      }
    } else if (length(combined_ld_chunks) == 1) {
      combined_ld = combined_ld_chunks[[1]]
    } else {
      # No LD chunks were created - likely due to errors
      message(sprintf("  No valid LD data retrieved for chromosome %s. Skipping.", chr))
      processed_chr <- c(processed_chr, chr)
      saveRDS(assoc_data_sorted, file = output_file) # Save progress
      next
    }
    
    # Process SNPs by p-value order (most significant first)
    chr_indices <- which(assoc_data_sorted$CHR == chr & is.na(assoc_data_sorted$clump))
    
    # Sort by P-value - crucial for getting the right index SNPs
    sorted_indices <- chr_indices[order(assoc_data_sorted$P[chr_indices])]
    
    for (idx in sorted_indices) {
      # Skip if this SNP has already been assigned to a clump
      if (!is.na(assoc_data_sorted$clump[idx])) next
      
      snp_i <- assoc_data_sorted$rsID[idx]
      current_bp <- assoc_data_sorted$BP[idx]
      
      # Mark this SNP as an index SNP and assign it to a new clump
      assoc_data_sorted$clump[idx] <- clump_index
      assoc_data_sorted$index_snp[idx] <- TRUE
      
      message(sprintf("    Created new clump %d with index SNP %s", clump_index, snp_i))
      
      # Find SNPs in LD with this index SNP
      # Look for unassigned SNPs within 5Mb
      candidate_indices <- which(
        assoc_data_sorted$CHR == chr &
          is.na(assoc_data_sorted$clump) &
          abs(assoc_data_sorted$BP - current_bp) <= 5e6
      )
      
      clump_size <- 1  # Start with 1 (the index SNP)
      
      for (j in candidate_indices) {
        snp_j <- assoc_data_sorted$rsID[j]
        
        # Check if SNPs are in LD
        if (!is.na(snp_i) && !is.na(snp_j) && 
            snp_i %in% rownames(combined_ld) && 
            snp_j %in% colnames(combined_ld)) {
          
          ld_value <- combined_ld[snp_i, snp_j]
          
          if (!is.na(ld_value) && ld_value > 0.05) {
            # Assign to the current clump but mark as NOT an index SNP
            assoc_data_sorted$clump[j] <- clump_index
            assoc_data_sorted$index_snp[j] <- FALSE
            clump_size <- clump_size + 1
          }
        }
      }
      
      message(sprintf("    Clump %d has %d SNPs in total", clump_index, clump_size))
      
      # Move to the next clump
      clump_index <- clump_index + 1
    }
    
    processed_chr <- c(processed_chr, chr)
    saveRDS(assoc_data_sorted, file = output_file) # Save progress
  }
  
  # Final summary
  clump_summary <- assoc_data_sorted %>%
    filter(!is.na(clump)) %>%
    group_by(clump) %>%
    summarize(
      n_snps = n(),
      has_index = any(index_snp),
      chr = CHR[1],
      index_snp = ifelse(any(index_snp), rsID[which(index_snp)[1]], NA)
    )
  
  n_clumps <- nrow(clump_summary)
  n_index_snps <- sum(clump_summary$has_index)
  
  cat(sprintf("Clumping complete. Processed %i chromosomes.\n", length(processed_chr)))
  cat(sprintf("Found %i clumps with %i index SNPs.\n", n_clumps, n_index_snps))
  
  return(assoc_data_sorted)
}
