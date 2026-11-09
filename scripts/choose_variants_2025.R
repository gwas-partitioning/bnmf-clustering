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
#   read_primary_gwas_at_positions() — Read the primary GWAS at given positions in the
#                                  standard form (shared by sentinels and proxies)
#   snp_clump()                  — Position-based clumping (no external API required)
#   ld_pruning_SNP.clip()        — LD pruning via LDlink SNPclip (multi-population)
#   find_variants_needing_proxies() — Flag strand-ambiguous or high-missingness sentinels
#   find_proxy_candidates()      — Proxy stage 1: LD partners of flagged sentinels from LDlink
#   qualify_proxy_candidates()   — Proxy stage 2: judge each partner against the primary GWAS
#   select_proxies()             — Proxy stage 3: choose one proxy per sentinel (no rescue)
#
# Legacy, kept for older run scripts and not used by main_script_example.R:
#   prune(), ld_pruning(), ld_pruning_topLD_api()   — earlier LD pruning routines
#   count_traits_per_variant(), count_traits_per_variant_2025() — trait counts per variant
#   window_to_sentinels()        — Restrict a pre-fetched proxy pool to windows around sentinels
#   choose_proxies()             — Single-call proxy search over a pre-fetched pool
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


read_primary_gwas_at_positions <- function(main_ss_filepath, positions,
                                           rename_cols = NULL, effect_allele = NULL) {
  # Read the primary GWAS at the given positions ("1:12345") and return one row
  # per variant in the pipeline's standard form:
  #   VAR_KEY, VAR_ID, CHR, POS, REF, ALT (the primary's allele order),
  #   BETA_main (effect of ALT), SE_main, PVALUE_main,
  #   Risk_Allele_main, Nonrisk_Allele_main
  # Used by get_biggest_gwas() and by the proxy stage, so sentinels and proxy
  # candidates are read from the primary GWAS in exactly the same way.

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
    CHR = main_std$CHR,
    POS = main_std$POS,
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
  main
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
  main <- read_primary_gwas_at_positions(main_ss_filepath, positions,
                                         rename_cols = rename_cols,
                                         effect_allele = effect_allele)
  # CHR and POS are taken from the discovery records in the join below.
  main$CHR <- NULL
  main$POS <- NULL

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
      # Kept so the QC report can re-check that the two GWAS agree.
      Risk_Allele_discovery,
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

# -----------------------------------------------------------------------------
# Proxy stage
# -----------------------------------------------------------------------------
# A sentinel flagged by find_variants_needing_proxies() is replaced by a variant
# in strong LD with it, or dropped if none qualifies. There is no rescue: a
# flagged sentinel is never kept as itself.
#
# The stage is three steps, called in turn by the run script:
#   1. find_proxy_candidates()     ask LDlink for LD partners of the flagged
#                                  sentinels (the only step that needs a network)
#   2. qualify_proxy_candidates()  judge each candidate against the primary
#                                  GWAS with the same rules a sentinel must pass
#   3. select_proxies()            given trait data for the qualified candidates
#                                  (a second fetch_summary_stats() call), choose
#                                  one proxy per sentinel
# Steps 2 and 3 need no network and are covered by the unit tests.

query_ldproxy_batch <- function(snps, population, token, out_file,
                                genome_build = "grch37") {
  # Ask LDlink for the LD partners of each variant in `snps` ("chr1:12345") and
  # append the answers to `out_file`. LDlinkR::LDproxy_batch() writes to a fixed
  # file name in the working directory and silently skips any variant LDlink
  # returns an error for; find_proxy_candidates() detects and reports those.
  expected_name <- sprintf("combined_query_snp_list_%s.txt", genome_build)
  if (basename(out_file) != expected_name) {
    stop("out_file must be named ", expected_name)
  }
  calling_dir <- getwd()
  on.exit(setwd(calling_dir), add = TRUE)
  setwd(dirname(out_file))
  tryCatch(
    LDlinkR::LDproxy_batch(snps, pop = population, r2d = "r2", token = token,
                           append = TRUE, genome_build = genome_build),
    error = function(e) {
      message("LDlink request stopped early: ", conditionMessage(e),
              " (answers received so far are kept).")
    }
  )
  invisible(NULL)
}

find_proxy_candidates <- function(need_proxies, token, population = "EUR",
                                  output_dir = ".", r2_min = 0.8,
                                  genome_build = "grch37", query_fun = NULL) {
  # Step 1 of the proxy stage.
  #
  # need_proxies: data frame with VAR_ID (CHR_POS_REF_ALT) for the flagged sentinels.
  # output_dir  : where LDlink's answers are kept. A rerun reuses the answers on
  #               disk and asks again only for sentinels with no answer, so a
  #               failed or rate-limited request is retried, never recorded as
  #               "no proxies".
  # query_fun   : function(snps, population, token, out_file, genome_build) that
  #               appends LDlink answers to out_file; defaults to
  #               query_ldproxy_batch(). Replaceable for testing.
  #
  # Returns one row per (sentinel, candidate) with r2 >= r2_min, biallelic SNPs
  # only, excluding the sentinel itself:
  #   VAR_ID, query_snp, proxy_rsID, proxy_ChrPos, proxy_A1, proxy_A2, r2,
  #   Correlated_Alleles (LDlink's map from each sentinel allele to the proxy
  #   allele it travels with, e.g. "A=C,G=T")
  # attr(, "unanswered") lists the sentinels LDlink gave no answer for.

  if (!"VAR_ID" %in% names(need_proxies)) stop("need_proxies must have a VAR_ID column.")
  if (is.null(query_fun)) query_fun <- query_ldproxy_batch

  empty <- data.frame(
    VAR_ID = character(0), query_snp = character(0), proxy_rsID = character(0),
    proxy_ChrPos = character(0), proxy_A1 = character(0), proxy_A2 = character(0),
    r2 = numeric(0), Correlated_Alleles = character(0), stringsAsFactors = FALSE
  )
  attr(empty, "unanswered") <- character(0)
  if (nrow(need_proxies) == 0L) return(empty)

  id_parts <- data.table::tstrsplit(need_proxies$VAR_ID, "_", fixed = TRUE)
  queries <- data.frame(
    VAR_ID = need_proxies$VAR_ID,
    query_snp = paste0("chr", id_parts[[1]], ":", id_parts[[2]]),
    stringsAsFactors = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  out_file <- file.path(output_dir, sprintf("combined_query_snp_list_%s.txt", genome_build))
  settings_file <- file.path(output_dir, "ldproxy_settings.rds")
  settings <- list(population = population, genome_build = genome_build, r2d = "r2")
  if (file.exists(out_file) && file.size(out_file) > 0) {
    if (!file.exists(settings_file) || !identical(readRDS(settings_file), settings)) {
      stop("LDlink answers already in ", output_dir, " were made with different settings ",
           "(population or genome build), or their settings are unknown. Use a different ",
           "output_dir for this search.")
    }
  } else {
    saveRDS(settings, settings_file)
  }

  read_answers <- function() {
    if (!file.exists(out_file) || file.size(out_file) == 0) return(NULL)
    # LDproxy_batch writes row numbers without a header for them, so the header
    # has one field fewer than the rows; row.names = NULL keeps the columns aligned.
    read.table(out_file, sep = "\t", header = TRUE, row.names = NULL,
               stringsAsFactors = FALSE, quote = "", comment.char = "")
  }
  answered_snps <- function(px) if (is.null(px)) character(0) else unique(px$query_snp)

  all_snps <- unique(queries$query_snp)
  for (attempt in 1:2) {
    pending <- setdiff(all_snps, answered_snps(read_answers()))
    if (length(pending) == 0L) break
    message(sprintf(
      "Asking LDlink (%s) for LD partners of %d flagged sentinel(s)%s...",
      population, length(pending), if (attempt == 2L) ", second attempt" else ""
    ))
    query_fun(pending, population, token, out_file, genome_build)
  }

  px <- read_answers()
  unanswered_snps <- setdiff(all_snps, answered_snps(px))
  unanswered <- queries$VAR_ID[queries$query_snp %in% unanswered_snps]

  candidates <- empty
  if (!is.null(px)) {
    px <- px[px$query_snp %in% all_snps, , drop = FALSE]
    px$R2 <- suppressWarnings(as.numeric(px$R2))
    is_snp <- grepl("^\\([ACGT]/[ACGT]\\)$", px$Alleles)
    keep <- px$Coord != px$query_snp & !is.na(px$R2) & px$R2 >= r2_min & is_snp
    px <- px[keep, , drop = FALSE]
    if (nrow(px) > 0L) {
      found <- data.frame(
        query_snp = px$query_snp,
        proxy_rsID = ifelse(is.na(px$RS_Number) | px$RS_Number %in% c(".", ""),
                            NA_character_, as.character(px$RS_Number)),
        proxy_ChrPos = sub("^chr", "", px$Coord),
        proxy_A1 = substr(px$Alleles, 2, 2),
        proxy_A2 = substr(px$Alleles, 4, 4),
        r2 = px$R2,
        Correlated_Alleles = as.character(px$Correlated_Alleles),
        stringsAsFactors = FALSE
      )
      candidates <- merge(queries, found, by = "query_snp", sort = FALSE)
      candidates <- candidates[!duplicated(candidates[, c("VAR_ID", "proxy_ChrPos",
                                                          "proxy_A1", "proxy_A2")]), ]
      candidates <- candidates[order(match(candidates$VAR_ID, queries$VAR_ID),
                                     -candidates$r2), names(empty)]
      rownames(candidates) <- NULL
    }
  }

  message(sprintf(
    "LDlink: %d of %d flagged sentinels answered; %d candidate(s) with r2 >= %g for %d of them.",
    length(all_snps) - length(unanswered_snps), length(all_snps),
    nrow(candidates), r2_min, length(unique(candidates$VAR_ID))
  ))
  if (length(unanswered) > 0L) {
    warning(sprintf(
      "LDlink returned no answer for %d flagged sentinel(s) (first: %s). They get no proxy in this run; rerun to ask again.",
      length(unanswered), paste(head(unanswered, 5), collapse = ", ")
    ), call. = FALSE)
  }
  attr(candidates, "unanswered") <- unanswered
  candidates
}

qualify_proxy_candidates <- function(candidates, main_ss_filepath, sentinels,
                                     main_pval_cutoff = 0.05,
                                     rename_cols = NULL, effect_allele = NULL,
                                     require_phase_consistency = TRUE) {
  # Step 2 of the proxy stage: judge every candidate against the primary GWAS.
  #
  # sentinels: data frame of ALL sentinels (flagged or not) with VAR_ID and
  #            Risk_Allele, as returned by get_biggest_gwas().
  #
  # A candidate is "qualified" when it
  #   - is not itself a sentinel (kept or flagged)
  #   - is not strand-ambiguous (A/T or C/G)
  #   - is in the primary GWAS (matched on position + allele pair)
  #   - has a direction of effect and a P <= main_pval_cutoff there
  #   - has a risk allele that travels with the sentinel's risk allele
  #     (LDlink's Correlated_Alleles), so the proxy represents the same signal
  #     in the same direction. Skipped when LDlink gives no allele map for the
  #     sentinel's risk allele, or when require_phase_consistency = FALSE.
  # Every other candidate is returned too, with the reason in `status`.
  #
  # Returns the candidates with: proxy_VAR_ID, CHR, POS, REF, ALT (the primary's
  # allele order), Risk_Allele, Nonrisk_Allele, PVALUE_main, BETA_main, SE_main,
  # status. attr(, "unanswered") is carried over from find_proxy_candidates().

  needed <- c("VAR_ID", "proxy_ChrPos", "proxy_A1", "proxy_A2", "r2")
  if (!all(needed %in% names(candidates))) {
    stop("candidates must be the output of find_proxy_candidates(); missing: ",
         paste(setdiff(needed, names(candidates)), collapse = ", "))
  }
  if (!all(c("VAR_ID", "Risk_Allele") %in% names(sentinels))) {
    stop("sentinels must have VAR_ID and Risk_Allele columns.")
  }
  unanswered <- attr(candidates, "unanswered")
  out <- as.data.frame(candidates, stringsAsFactors = FALSE)

  pos_parts <- data.table::tstrsplit(out$proxy_ChrPos, ":", fixed = TRUE)
  if (length(pos_parts) < 2L) pos_parts <- list(character(0), character(0))
  out$CHR <- suppressWarnings(as.integer(pos_parts[[1]]))
  out$POS <- suppressWarnings(as.integer(pos_parts[[2]]))
  out$VAR_KEY <- canonical_variant_key(out$CHR, out$POS, out$proxy_A1, out$proxy_A2)

  if (nrow(out) > 0L) {
    main <- read_primary_gwas_at_positions(main_ss_filepath, unique(out$proxy_ChrPos),
                                           rename_cols = rename_cols,
                                           effect_allele = effect_allele)
  } else {
    main <- data.frame(VAR_KEY = character(0), VAR_ID = character(0), REF = character(0),
                       ALT = character(0), BETA_main = numeric(0), SE_main = numeric(0),
                       PVALUE_main = numeric(0), Risk_Allele_main = character(0),
                       Nonrisk_Allele_main = character(0), stringsAsFactors = FALSE)
  }
  hit <- match(out$VAR_KEY, main$VAR_KEY)
  out$proxy_VAR_ID   <- main$VAR_ID[hit]
  out$REF            <- main$REF[hit]
  out$ALT            <- main$ALT[hit]
  out$Risk_Allele    <- main$Risk_Allele_main[hit]
  out$Nonrisk_Allele <- main$Nonrisk_Allele_main[hit]
  out$PVALUE_main    <- main$PVALUE_main[hit]
  out$BETA_main      <- main$BETA_main[hit]
  out$SE_main        <- main$SE_main[hit]

  # The proxy allele that travels with the sentinel's risk allele, per LDlink.
  sentinel_risk <- sentinels$Risk_Allele[match(out$VAR_ID, sentinels$VAR_ID)]
  corr <- if ("Correlated_Alleles" %in% names(out)) out$Correlated_Alleles else rep(NA_character_, nrow(out))
  expected_risk <- vapply(seq_len(nrow(out)), function(i) {
    if (is.na(corr[i]) || is.na(sentinel_risk[i])) return(NA_character_)
    pairs <- strsplit(strsplit(corr[i], ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
    pairs <- pairs[vapply(pairs, length, integer(1)) == 2L]
    if (length(pairs) == 0L) return(NA_character_)
    from <- vapply(pairs, `[`, character(1), 1L)
    to   <- vapply(pairs, `[`, character(1), 2L)
    matched <- to[match(sentinel_risk[i], from)]
    if (is.na(matched) || !grepl("^[ACGT]$", matched)) NA_character_ else matched
  }, character(1))
  opposes <- require_phase_consistency & !is.na(expected_risk) &
    !is.na(out$Risk_Allele) & expected_risk != out$Risk_Allele

  sentinel_parts <- data.table::tstrsplit(sentinels$VAR_ID, "_", fixed = TRUE)
  sentinel_chrpos <- paste(sentinel_parts[[1]], sentinel_parts[[2]], sep = ":")
  palindromic <- paste0(out$proxy_A1, out$proxy_A2) %in% c("AT", "TA", "CG", "GC")

  out$status <- dplyr::case_when(
    out$proxy_ChrPos %in% sentinel_chrpos                      ~ "is_a_sentinel",
    palindromic                                                ~ "strand_ambiguous",
    is.na(out$proxy_VAR_ID)                                    ~ "not_in_main_gwas",
    is.na(out$Risk_Allele)                                     ~ "no_main_direction",
    is.na(out$PVALUE_main)                                     ~ "main_pvalue_missing",
    out$PVALUE_main > main_pval_cutoff                         ~ "main_pvalue_above_cutoff",
    opposes                                                    ~ "risk_allele_opposes_sentinel",
    TRUE                                                       ~ "qualified"
  )

  out <- out[, c("VAR_ID", intersect("query_snp", names(out)), "proxy_rsID", "proxy_ChrPos",
                 "proxy_VAR_ID", "CHR", "POS", "REF", "ALT", "r2", "Risk_Allele",
                 "Nonrisk_Allele", "PVALUE_main", "BETA_main", "SE_main", "status"),
             drop = FALSE]
  if (nrow(out) > 0L) {
    counts <- table(out$status)
    message("Proxy candidates against the primary GWAS: ",
            paste(sprintf("%s = %d", names(counts), as.integer(counts)), collapse = "; "))
  }
  attr(out, "unanswered") <- if (is.null(unanswered)) character(0) else unanswered
  out
}

select_proxies <- function(need_proxies, qualified, z_candidates,
                           min_nonmissing = 0.8,
                           unanswered = attr(qualified, "unanswered")) {
  # Step 3 of the proxy stage: choose one proxy per flagged sentinel.
  #
  # need_proxies : the flagged sentinels (VAR_ID; PVALUE_main if available)
  # qualified    : output of qualify_proxy_candidates()
  # z_candidates : z-scores for the qualified candidates over the retained
  #                traits (rows named "Chr:Pos"), from a second
  #                fetch_summary_stats() call
  #
  # A candidate must be qualified and have data for at least `min_nonmissing`
  # of the traits. Among those, the one in strongest LD with the sentinel is
  # chosen (then most trait data, then lowest primary-GWAS P). A variant is the
  # proxy for at most one sentinel; sentinels choose in order of primary-GWAS P.
  # A flagged sentinel with no proxy is dropped: there is no rescue.
  #
  # Returns a list:
  #   chosen     one row per sentinel that got a proxy
  #   no_proxy   flagged sentinels with no proxy, with the reason
  #   candidates every candidate with its final status

  if (!"VAR_ID" %in% names(need_proxies)) stop("need_proxies must have a VAR_ID column.")
  if (is.null(unanswered)) unanswered <- character(0)
  audit <- as.data.frame(qualified, stringsAsFactors = FALSE)
  if (!"status" %in% names(audit)) stop("qualified must be the output of qualify_proxy_candidates().")

  coverage <- if (is.null(z_candidates) || nrow(z_candidates) == 0L) {
    numeric(0)
  } else {
    rowMeans(!is.na(as.matrix(z_candidates)))
  }
  audit$coverage <- unname(coverage[audit$proxy_ChrPos])
  too_sparse <- audit$status == "qualified" &
    (is.na(audit$coverage) | audit$coverage < min_nonmissing)
  audit$status[too_sparse] <- "insufficient_trait_data"

  p_main <- if ("PVALUE_main" %in% names(audit)) audit$PVALUE_main else rep(0, nrow(audit))
  sentinel_order <- if ("PVALUE_main" %in% names(need_proxies)) {
    need_proxies$VAR_ID[order(need_proxies$PVALUE_main)]
  } else {
    need_proxies$VAR_ID
  }

  used <- character(0)
  blocked <- character(0)
  chosen_rows <- integer(0)
  for (sentinel in sentinel_order) {
    idx <- which(audit$VAR_ID == sentinel & audit$status == "qualified")
    if (length(idx) == 0L) next
    idx <- idx[order(-audit$r2[idx], -audit$coverage[idx], p_main[idx], audit$proxy_ChrPos[idx])]
    free <- idx[!audit$proxy_ChrPos[idx] %in% used]
    if (length(free) == 0L) {
      blocked <- c(blocked, sentinel)
      next
    }
    chosen_rows <- c(chosen_rows, free[1])
    used <- c(used, audit$proxy_ChrPos[free[1]])
  }
  audit$status[chosen_rows] <- "chosen"
  audit$status[audit$status == "qualified"] <- "qualified_not_chosen"

  chosen <- audit[chosen_rows, , drop = FALSE]
  names(chosen)[names(chosen) == "VAR_ID"] <- "original_SNP"
  chosen <- chosen[, intersect(c("original_SNP", "proxy_VAR_ID", "proxy_ChrPos", "proxy_rsID",
                                 "r2", "coverage", "CHR", "POS", "REF", "ALT", "Risk_Allele",
                                 "Nonrisk_Allele", "PVALUE_main", "BETA_main", "SE_main"),
                               names(chosen)), drop = FALSE]
  rownames(chosen) <- NULL

  without <- setdiff(need_proxies$VAR_ID, chosen$original_SNP)
  no_proxy <- data.frame(
    VAR_ID = without,
    reason = vapply(without, function(sentinel) {
      rows <- audit[audit$VAR_ID == sentinel, , drop = FALSE]
      if (sentinel %in% unanswered) {
        "ldlink_no_answer"
      } else if (nrow(rows) == 0L) {
        "no_ld_candidate"
      } else if (sentinel %in% blocked) {
        "candidate_used_by_another_sentinel"
      } else if (any(rows$status == "insufficient_trait_data")) {
        "no_candidate_with_enough_trait_data"
      } else {
        "no_candidate_passed_primary_filters"
      }
    }, character(1)),
    stringsAsFactors = FALSE
  )
  rownames(no_proxy) <- NULL

  message(sprintf("Proxies: %d flagged sentinel(s); proxy chosen for %d; %d dropped with no proxy.",
                  nrow(need_proxies), nrow(chosen), nrow(no_proxy)))
  if (nrow(no_proxy) > 0L) {
    counts <- table(no_proxy$reason)
    message("  No proxy because: ",
            paste(sprintf("%s = %d", names(counts), as.integer(counts)), collapse = "; "))
  }

  list(chosen = chosen, no_proxy = no_proxy, candidates = audit)
}

# Legacy single-call proxy search, kept for run scripts that pre-fetch a pool of
# proxy candidates (zmat_fullset must already hold every candidate). New run
# scripts should use the three functions above. Like them, it does not rescue:
# a flagged sentinel with no proxy is dropped.
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
  # FINAL SETS
  # ==========================================================
  # Originals that never needed a proxy are kept. A flagged sentinel is kept
  # only through a proxy; one with no acceptable proxy is dropped (no rescue).
  final_originals_to_keep <- setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID)
  dropped_without_proxy <- setdiff(
    need_proxies$VAR_ID,
    if (!is.null(final_proxy_df) && nrow(final_proxy_df) > 0) final_proxy_df$VAR_ID else character(0)
  )
  message(sprintf("Dropped %d flagged variant(s) with no acceptable proxy.",
                  length(dropped_without_proxy)))

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
