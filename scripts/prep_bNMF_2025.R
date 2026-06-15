# =============================================================================
# prep_bNMF_2025.R
# =============================================================================
# Helper functions for fetching GWAS summary statistics and preparing the
# z-score matrix used as input to the bNMF clustering algorithm.
#
# Functions:
#   peek_headers()          — Read column names from the first line of a file (gz-aware)
#   read_single_trait_dt()  — Fetch and align summary stats for one trait via grep/tabix
#   fetch_summary_stats()   — Fetch and align summary stats across all traits; returns z- and N-matrices
#   prep_z_matrix()         — Filter traits, prune by correlation, scale z-scores, expand to non-negative form
#   fill_missing_zscores()  — (legacy) Fill missing z-scores using proxy variants
#
# Assumptions:
#   - Input variants use "Chr:Pos" (colon-separated) as the standard SNP key
#   - Trait summary stat files contain at minimum: a SNP ID column, BETA, SE
#   - Genome build: hg19/GRCh37
# =============================================================================

library(tidyverse)
library(data.table)
library(future.apply)
library(stringr)

# Helper: safe read of column headers from a file (supports .gz)
peek_headers <- function(filepath) {
  if (endsWith(filepath, ".gz")) {
    hdr <- readLines(gzfile(filepath), n = 1)
  } else {
    hdr <- readLines(filepath, n = 1)
  }
  strsplit(trimws(hdr), "\\s+")[[1]]
}

# Rewritten, robust read_single_trait_dt
read_single_trait_dt <- function(trait,
                                 read_trait_input,
                                 trait_ss_files,
                                 trait_ss_size = NULL,
                                 tmp_var_file = "variants_to_query.tmp") {
  
  message(sprintf("Processing %s...", trait))
  
  # ---- 0. Basic checks ----
  if (!trait %in% names(trait_ss_files)) {
    stop("trait not found in trait_ss_files")
  }
  trait_file <- trait_ss_files[[trait]]
  
  # ---- 1. Inspect trait file header to detect ID column ----
  headers <- peek_headers(trait_file)
  id_col_in_file <- if ("SNP" %in% headers) "SNP" else if ("VAR_ID" %in% headers) "VAR_ID" else NULL
  if (is.null(id_col_in_file)) {
    message(sprintf("  No SNP or VAR_ID column found in %s", trait))
    return(NULL)
  }
  
  # ---- 2. Input format detection (from read_trait_input$SNP) ----
  # Ensure read_trait_input is a data.table
  if (!is.data.table(read_trait_input)) read_trait_input <- as.data.table(read_trait_input)
  
  input_snps <- read_trait_input$SNP
  if (length(input_snps) == 0) {
    message("  No input SNPs provided")
    return(NULL)
  }
  sample_input <- input_snps[1]
  
  input_has_chr_prefix        <- grepl("^chr", sample_input, ignore.case = TRUE)
  input_is_chrpos_colon       <- grepl("^(chr)?\\d+:\\d+$", sample_input, ignore.case = TRUE)
  input_is_chrpos_underscore  <- grepl("^\\d+_\\d+$", sample_input)
  input_is_varid              <- grepl("^\\d+_\\d+_[ACGT]+_[ACGT]+$", sample_input, ignore.case = FALSE)
  
  message(sprintf("  Input format: chr_prefix=%s, colon=%s, underscore=%s, varid=%s",
                  input_has_chr_prefix, input_is_chrpos_colon, input_is_chrpos_underscore, input_is_varid))
  
  # ---- 3. Build search_snps to match trait file format ----
  # We'll create a vector search_snps that maps 1:1 to input_snps (original -> search)
  search_snps <- input_snps
  
  # If trait uses VAR_ID, we want "CHR_POS" (no alleles)
  # If trait uses colon format, we want "chrN:pos" or "N:pos" according to trait prefix
  # If trait uses underscore format, we want "N_pos"
  # We'll detect trait sample to know what's needed
  # Read a tiny sample from trait file to detect its ID style
  sample_trait_id <- fread(cmd = if (endsWith(trait_file, ".gz")) sprintf("gzip -cd %s | head -n 5", trait_file) else sprintf("head -n 5 %s", trait_file),
                           select = id_col_in_file, data.table = FALSE, showProgress = FALSE)[1,1]
  
  trait_has_chr_prefix       <- grepl("^chr", sample_trait_id, ignore.case = TRUE)
  trait_is_chrpos_colon      <- grepl("^(chr)?\\d+:\\d+$", sample_trait_id, ignore.case = TRUE)
  trait_is_chrpos_underscore <- grepl("^\\d+_\\d+$", sample_trait_id)
  trait_is_varid             <- grepl("^\\d+_\\d+_[ACGT]+_[ACGT]+$", sample_trait_id)
  
  message(sprintf("  Trait file format (%s): chr_prefix=%s, colon=%s, underscore=%s, varid=%s",
                  id_col_in_file, trait_has_chr_prefix, trait_is_chrpos_colon, trait_is_chrpos_underscore, trait_is_varid))
  message(sprintf("  Sample from trait file: %s", sample_trait_id))
  
  # Convert input -> search format
  # Start from canonical forms derived from input:
  # - make chr_pos_underscore: "1_12345"
  # - make chr_pos_colon: "1:12345" (optionally "chr1:12345")
  canonical_underscore <- if (input_is_chrpos_colon) gsub(":", "_", gsub("^chr", "", input_snps, ignore.case = TRUE)) else if (input_is_chrpos_underscore) input_snps else if (input_is_varid) sub("_[ACGT]+_[ACGT]+$", "", input_snps) else input_snps
  canonical_colon      <- gsub("_", ":", canonical_underscore)
  canonical_chr_colon  <- ifelse(grepl("^chr", input_snps, ignore.case = TRUE), gsub("_", ":", input_snps), paste0("chr", canonical_colon))
  
  if (trait_is_varid) {
    # trait uses VAR_ID: search should be chr_pos_underscore without alleles (e.g., "1_2189477")
    search_snps <- canonical_underscore
  } else if (trait_is_chrpos_colon) {
    # trait uses colon. Respect trait prefix preference
    if (trait_has_chr_prefix) {
      # needs "chr1:12345"
      # If input lacks "chr", add it
      search_snps <- ifelse(grepl("^chr", canonical_colon, ignore.case = TRUE), canonical_colon, paste0("chr", canonical_colon))
    } else {
      # needs "1:12345" (no chr)
      search_snps <- canonical_colon
    }
  } else if (trait_is_chrpos_underscore) {
    # trait uses underscore format "1_12345"
    search_snps <- canonical_underscore
  } else {
    # fallback: try matching exact input if nothing detected
    search_snps <- input_snps
  }
  
  # Map original -> search
  snp_mapping <- data.table(original_snp = input_snps, search_snp = search_snps)
  
  message(sprintf("  Converted search terms. Example: %s -> %s",
                  snp_mapping$original_snp[1], snp_mapping$search_snp[1]))
  
  # Write mapping to temp file used by grep (overwriting tmp_var_file)
  # Caller may already have written one, but make sure it's the correct set for this trait.
  writeLines(unique(snp_mapping$search_snp), tmp_var_file)
  
  # ---- 4. Grep/read trait file for matches ----
  cmd_str <- if (endsWith(trait_file, ".gz")) {
    sprintf("gzip -cd %s | grep -Ff %s", trait_file, tmp_var_file)
  } else {
    sprintf("grep -Ff %s %s", tmp_var_file, trait_file)
  }
  
  dt <- tryCatch({
    # read with names from header (we already peeked headers)
    fread(cmd = cmd_str, header = FALSE, col.names = headers, data.table = TRUE, showProgress = FALSE)
  }, error = function(e) {
    message("  Grep returned no results or fread failed: ", e$message)
    return(data.table())
  })
  
  if (nrow(dt) == 0) {
    message(sprintf("  No SNP matches found for %s", trait))
    return(NULL)
  }
  message(sprintf("  Found %d rows after grep", nrow(dt)))
  
  # ---- 5. Harmonize column names EARLY ----
  # Standard allele column canonicalization: use EA and NEA internally
  if ("Effect_Allele_PH" %in% names(dt)) setnames(dt, "Effect_Allele_PH", "EA")
  if ("Effect_Allele"    %in% names(dt)) setnames(dt, "Effect_Allele", "EA")
  if ("effect_allele"    %in% names(dt)) setnames(dt, "effect_allele", "EA")
  if ("A1" %in% names(dt) && !("EA" %in% names(dt))) setnames(dt, "A1", "EA")
  
  if ("NEA" %in% names(dt)) {
    # ok
  } else if ("A2" %in% names(dt)) {
    setnames(dt, "A2", "NEA")
  }
  
  # Standardize P/BETA/SE names
  if ("p" %in% names(dt)) setnames(dt, "p", "P_VALUE")
  if ("pval" %in% names(dt)) setnames(dt, "pval", "P_VALUE")
  if ("N" %in% names(dt) && !("N_PH" %in% names(dt))) setnames(dt, "N", "N_PH")
  if ("Effect" %in% names(dt) && !("BETA" %in% names(dt))) setnames(dt, "Effect", "BETA")
  
  # ---- 6. Extract NEA from VAR_ID if missing ----
  if (!"NEA" %in% names(dt)) {
    if ("VAR_ID" %in% names(dt)) {
      
      # VAR_ID format: CHR_POS_REF_ALT
      tmp_cols <- tstrsplit(dt$VAR_ID, "_", fixed = TRUE)
      
      if (length(tmp_cols) < 4) {
        message("  ERROR: VAR_ID does not have 4 fields for ", trait)
        return(NULL)
      }
      
      dt[, chr_from_varid := tmp_cols[[1]]]
      dt[, pos_from_varid := tmp_cols[[2]]]
      dt[, ref_from_varid := tmp_cols[[3]]]
      dt[, alt_from_varid := tmp_cols[[4]]]
      
      # assign NEA based on EA vs ref/alt
      if ("EA" %in% names(dt)) {
        dt[, NEA := fifelse(toupper(EA) == toupper(ref_from_varid), alt_from_varid,
                            fifelse(toupper(EA) == toupper(alt_from_varid), ref_from_varid, NA_character_))]
      } else {
        dt[, NEA := alt_from_varid]
      }
      
      # chr_pos key
      dt[, varid_chrpos_us := paste0(chr_from_varid, "_", pos_from_varid)]
    } else {
      message("  ERROR: No NEA column and no VAR_ID to extract from")
      return(NULL)
    }
  }  
  # ---- 7. Build match_key on trait dt so we can join back to original search names ----
  # If trait file uses VAR_ID, match_key should be chr_pos underscore
  if (id_col_in_file == "VAR_ID" || trait_is_varid) {
    # prefer varid_chrpos_us if present, else try to derive from CHR+POS columns
    if ("varid_chrpos_us" %in% names(dt)) {
      dt[, match_key := varid_chrpos_us]
    } else if (all(c("CHR", "POS") %in% names(dt))) {
      dt[, match_key := paste0(CHR, "_", POS)]
    } else {
      # last resort: derive from VAR_ID
      if ("VAR_ID" %in% names(dt)) {
        dt[, match_key := sub("_[^_]+_[^_]+$", "", VAR_ID)]  # drop _REF_ALT
      } else {
        dt[, match_key := NA_character_]
      }
    }
  } else {
    # trait uses explicit SNP column (could be chr:pos or chr_pos or rsID)
    dt[, match_key := get(id_col_in_file)]
  }
  
  # ---- 8. Join dt to mapping table to restore original input SNP IDs ----
  setDT(snp_mapping)
  setkey(snp_mapping, search_snp)
  dt_merged <- dt[snp_mapping, on = .(match_key = search_snp), nomatch = 0]
  dt <- copy(dt_merged)
  setnames(dt, "original_snp", "SNP")
  
  # ───── SPECIAL HANDLING: GLGC / older consortia with no EA/NEA columns ─────
  if (!"EA" %in% names(dt) && "VAR_ID" %in% names(dt)) {
    message("  No EA column found → assuming GLGC/old-consortium format (effect = 2nd allele in VAR_ID)")
    dt[, c("ref_glgc", "alt_glgc") := tstrsplit(sub("^[^_]+_[^_]+_", "", VAR_ID), "_", fixed = FALSE)]
    dt[, EA  := alt_glgc]
    dt[, NEA := ref_glgc]
    dt[, c("ref_glgc", "alt_glgc") := NULL]
  }
  
  # Check row count using dt, not dt_merged
  if (nrow(dt) == 0) {
    message(sprintf("  No SNPs matched after key mapping for %s", trait))
    return(NULL)
  }
  
  
  # ---- 9. Ensure EA/NEA and BETA/SE present and consistent ----
  # At this stage, dt has EA and NEA (either originally or via VAR_ID extraction)
  if (!"EA" %in% names(dt)) {
    message(sprintf("  ERROR: No EA present for %s", trait))
    return(NULL)
  }
  if (!"NEA" %in% names(dt)) {
    message(sprintf("  ERROR: No NEA present for %s", trait))
    return(NULL)
  }
  
  # Remove rows where EA==NEA
  dt <- dt[toupper(EA) != toupper(NEA)]
  if (nrow(dt) == 0) return(NULL)
  
  # Standardize sample size column
  # Convert vector first — scalar := coerces to existing column type, so type must be fixed first
  if ("N_PH" %in% names(dt)) dt[, N_PH := as.numeric(N_PH)]
  if (!is.null(trait_ss_size) && trait %in% names(trait_ss_size) && !is.na(trait_ss_size[[trait]])) {
    dt[, N_PH := as.numeric(trait_ss_size[[trait]])]
  } else if (!"N_PH" %in% names(dt)) {
    dt[, N_PH := NA_real_]
  }
  
  # Standardize P_VALUE column numeric
  if ("P_VALUE" %in% names(dt)) {
    dt[, P_VALUE := as.numeric(P_VALUE)]
  } else if ("pval" %in% names(dt)) {
    dt[, P_VALUE := as.numeric(pval)]
  }
  
  # Standardize BETA/SE
  if (!("BETA" %in% names(dt)) || !("SE" %in% names(dt))) {
    message(sprintf("  Cannot compute z-score for %s without BETA/SE", trait))
    return(NULL)
  }
  dt[, SE := as.numeric(SE)]
  dt[, BETA := as.numeric(BETA)]
  
  # fix SE==0
  if (any(dt$SE == 0, na.rm = TRUE)) {
    min_pos_se <- min(dt$SE[dt$SE > 0], na.rm = TRUE)
    dt[SE == 0, SE := min_pos_se]
  }
  
  # ---- 10. Join with read_trait_input (gwas input alleles) to get Risk/Nonrisk alleles ----
  # ---- 10. Join with read_trait_input (gwas input alleles) to get Risk/Nonrisk alleles ----
  if (!all(c("Risk_Allele", "Nonrisk_Allele") %in% names(read_trait_input))) {
    stop("read_trait_input must have Risk_Allele and Nonrisk_Allele")
  }
  
  # FIX: Ensure both dt and read_trait_input have matching SNP formats
  dt_has_chr <- grepl("^chr", dt$SNP[1], ignore.case = TRUE)
  input_has_chr <- grepl("^chr", read_trait_input$SNP[1], ignore.case = TRUE)
  
  message(sprintf("  SNP format check: dt has chr=%s, input has chr=%s", 
                  dt_has_chr, input_has_chr))
  
  # Standardize to match read_trait_input format
  if (dt_has_chr && !input_has_chr) {
    message("  Removing 'chr' prefix from dt SNPs to match input")
    dt[, SNP := gsub("^chr", "", SNP, ignore.case = TRUE)]
  } else if (!dt_has_chr && input_has_chr) {
    message("  Adding 'chr' prefix to dt SNPs to match input")
    dt[, SNP := paste0("chr", SNP)]
  }
  
  # Convert to data.table and join
  if (!is.data.table(dt)) dt <- as.data.table(dt)
  if (!is.data.table(read_trait_input)) read_trait_input <- as.data.table(read_trait_input)
  
  dt <- dt[read_trait_input, on = "SNP", nomatch = 0]

  if (nrow(dt) == 0) {
    message(sprintf("  No SNPs after joining with input for %s", trait))
    return(NULL)
  }
  message(sprintf("  %d SNPs after joining with input", nrow(dt)))
  
  # ---- 11. Allele matching (unordered pair equality) ----
  dt[, trait_alleles := paste0(pmin(toupper(EA), toupper(NEA)), "_", pmax(toupper(EA), toupper(NEA)))]
  dt[, gwas_alleles  := paste0(pmin(toupper(Risk_Allele), toupper(Nonrisk_Allele)), "_", pmax(toupper(Risk_Allele), toupper(Nonrisk_Allele)))]
  
  message(sprintf("  Before allele matching: %d SNPs", nrow(dt)))
  message(sprintf("  Sample trait_alleles: %s", dt$trait_alleles[1]))
  message(sprintf("  Sample gwas_alleles: %s", dt$gwas_alleles[1]))
  matches <- sum(dt$trait_alleles == dt$gwas_alleles, na.rm = TRUE)
  message(sprintf("  Allele matches: %d", matches))
  
  dt <- dt[trait_alleles == gwas_alleles]
  if (nrow(dt) == 0) {
    message(sprintf("  No SNPs with matching alleles found for %s", trait))
    return(NULL)
  }
  
  # ---- 12. Compute z and orient to Risk_Allele ----
  dt[, z := BETA / SE]
  
  # If Effect allele equals Risk_Allele then z as-is, else flip sign
  dt[, z := ifelse(toupper(EA) == toupper(Risk_Allele), z,
                   ifelse(toupper(NEA) == toupper(Risk_Allele), -z, NA_real_))]
  
  if (!"P_VALUE" %in% names(dt)) {
    message(sprintf("  No P_VALUE in source — deriving from z-score for %s", trait))
    dt[, P_VALUE := 2 * pnorm(-abs(z))]
  }
  dt <- dt[, .(SNP, z, N_PH, P_VALUE)]
  dt[, P_VALUE := as.numeric(P_VALUE)]

  n_z_na <- sum(is.na(dt$z))
  if (n_z_na > 0)
    warning(sprintf("[%s]: Dropping %d rows where z is NA (missing BETA or SE).", trait, n_z_na))
  dt <- dt[!is.na(z)]

  if (nrow(dt) > 0 && all(is.na(dt$N_PH)))
    warning(sprintf(paste0("[%s]: N_PH is NA for all rows. Provide sample size via trait_ss_size ",
                           "in the config Excel, or add an 'N' column to the formatted summary ",
                           "stats file. The N matrix for this trait will use NA."), trait))

  message(sprintf("  Final: %d SNPs for %s", nrow(dt), trait))
  
  # Clean up temporaries if present
  cols_to_drop <- intersect(c("ref_from_varid", "alt_from_varid", "chr_from_varid", "pos_from_varid", "varid_chrpos_us", "trait_alleles", "gwas_alleles"), names(dt))
  if (length(cols_to_drop)) dt[, (cols_to_drop) := NULL]
  
  return(dt)
}

fetch_summary_stats <- function(df_input, gwas_ss_file, trait_ss_files, trait_ss_size = NULL, pval_cutoff = 1, checkpoint_dir = NULL) {
  
  # GOAL: Fetches and aligns summary statistics using "Chr:Pos" as the standard SNP identifier.
  #
  # ASSUMPTIONS:
  #   - df_input: A data frame with 'SNP' (Chr:Pos format), 'Risk_Allele', 'REF', and 'ALT' columns.
  #   - gwas_ss_file: Can be a file path or a data frame.
  #     - If a data frame, it must contain 'SNP', 'REF', 'ALT', 'BETA', and 'P_VALUE' columns.
  #     - If a file, it must be searchable by 'SNP' (Chr:Pos) and contain the same required columns.
  #   - trait_ss_files: A named vector of filepaths. Each file must have 'SNP', 'EA', 'NEA', 'BETA', etc.
  #
  # OUTPUT: A list containing:
  #   - df_z: A data frame of aligned Z-scores (SNPs x Traits).
  #   - df_N: A data frame of sample sizes (SNPs x Traits).
  #   - df_gwas: The filtered and formatted data from the primary GWAS.
  

  # -------------------------------------------------------------------------- #
  # Main Function Body
  # -------------------------------------------------------------------------- #
  
  # --- 1. Process Primary GWAS ---
  
  print("Writing SNP list (chr:pos) to file for grepping...")
  tmp_var_file <- "variants_to_query.tmp"
  writeLines(df_input$SNP, tmp_var_file)
  
  if (is.character(gwas_ss_file)) {
    print("Reading primary GWAS summary statistics from file...")
    headers <- as.character(fread(gwas_ss_file, nrows = 1, data.table = F, header = F))
    cmd_grep <- if (endsWith(gwas_ss_file, ".gz")) {
      sprintf("gzip -cd %s | grep -Fwf %s", gwas_ss_file, tmp_var_file)
    } else {
      sprintf("grep -Fwf %s %s", gwas_ss_file, tmp_var_file)
    }
    gwas_ss <- fread(cmd = cmd_grep, header = F, col.names = headers, data.table = FALSE, stringsAsFactors = FALSE)
  } else {
    print("Filtering pre-loaded primary GWAS data frame...")
    gwas_ss <- gwas_ss_file %>% filter(SNP %in% df_input$SNP)
    
  }

  print(sprintf("%i of %i SNPs found in primary GWAS.", nrow(gwas_ss), nrow(df_input)))

  # Standardize primary GWAS: Derive Risk/Non-risk Alleles from BETA, REF, and ALT
  if (!all(c("Risk_Allele", "Nonrisk_Allele") %in% names(gwas_ss))) {
    print("Standardizing primary GWAS: Deriving Risk/Non-risk Alleles...")
    if (!"BETA" %in% colnames(gwas_ss) && "ODDS_RATIO" %in% colnames(gwas_ss)) {
      gwas_ss <- gwas_ss %>% mutate(BETA = log(ODDS_RATIO))
    }


    required_cols <- c("REF", "ALT", "BETA")
    if(!all(required_cols %in% names(gwas_ss))){
      stop(paste("Primary GWAS is missing required columns for standardization:", paste(setdiff(required_cols, names(gwas_ss)), collapse=", ")))
    }
    
    gwas_ss <- gwas_ss %>%
      mutate(
        Risk_Allele = if_else(BETA > 0, ALT, REF),
        Nonrisk_Allele = if_else(BETA > 0, REF, ALT)
      )
  }
  
  # --- 2. Merge with Input and Filter SNPs ---
  df_input_wGWAS <- df_input %>%
    dplyr::rename(Risk_Allele_Orig = Risk_Allele) %>%
    select(SNP, Risk_Allele_Orig) %>%
    inner_join(gwas_ss, by = "SNP") %>%
    mutate(across(any_of("P_VALUE"), as.numeric)) %>%
    mutate(P_VALUE = if_else(P_VALUE == 0, .Machine$double.xmin, P_VALUE))
  
  print(paste0(nrow(df_input_wGWAS), " variants available after merging with primary GWAS."))
  
  pval_bonf <- 0.05 / nrow(df_input_wGWAS)
  opp_risk <- df_input_wGWAS %>% filter(between(P_VALUE, pval_bonf, pval_cutoff) & Risk_Allele != Risk_Allele_Orig)
  high_pval <- df_input_wGWAS %>% filter(P_VALUE > pval_cutoff)
  uniq_to_drop <- unique(c(opp_risk$SNP, high_pval$SNP))
  
  df_input_wGWAS_filtered <- df_input_wGWAS %>%
    filter(!SNP %in% uniq_to_drop)
  
  print(sprintf("%i variants remain after p-value and risk allele filtering.", nrow(df_input_wGWAS_filtered)))
  
  df_read_trait_input <- df_input_wGWAS_filtered %>%
    dplyr::select(SNP, Risk_Allele, Nonrisk_Allele) %>%
    as.data.table()
  
  writeLines(df_read_trait_input$SNP, tmp_var_file)
  
  # --- 3. Process All Trait Files ---
  if (!is.null(checkpoint_dir))
    dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. Run lapply on the trait names. This creates an UNNAMED list.
  list_of_results <- lapply(names(trait_ss_files), function(trait) {
    cp_file <- if (!is.null(checkpoint_dir))
      file.path(checkpoint_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", trait), ".rds"))
    else NULL

    if (!is.null(cp_file) && file.exists(cp_file)) {
      message(sprintf("  [checkpoint] %s", trait))
      return(readRDS(cp_file))
    }

    result <- tryCatch(
      read_single_trait_dt(
        trait = trait,
        read_trait_input = df_read_trait_input,
        trait_ss_files = trait_ss_files,
        trait_ss_size = trait_ss_size,
        tmp_var_file = tmp_var_file
      ),
      error = function(e) {
        message(sprintf("  SKIPPING %s (error): %s", trait, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(cp_file)) saveRDS(result, cp_file)

    result
  })
  
  # 2. *** THE FIX: Immediately set the names of the list ***
  #    This ensures the list is correctly named before any other step.
  names(list_of_results) <- names(trait_ss_files)
  
  # 3. Now, the rest of the logic will work correctly.
  #    Filter out any NULL results from failed traits.
  valid_results_list <- Filter(Negate(is.null), list_of_results)
  
  # 4. Bind the rows. It now receives a named list.
  trait_df_list <- bind_rows(valid_results_list, .id = "trait")
  
  print(head(trait_df_list))
  # --- 4. Consolidate and Return Output ---
  if (nrow(trait_df_list) > 0) {
    
    trait_order <- unique(names(trait_ss_files))
    trait_df_list$trait <- factor(trait_df_list$trait, levels = trait_order)
    
    # Pivot Z-scores
    z_df_wide <- trait_df_list %>%
      select(trait, SNP, z) %>%
      # ADD THIS BLOCK to ensure each SNP/trait is unique before pivoting
      group_by(trait, SNP) %>% 
      dplyr::slice(1) %>% 
      ungroup() %>%
      pivot_wider(names_from = "trait", values_from = "z")
    
    z_mat <- as.matrix(z_df_wide[, -1])
    rownames(z_mat) <- z_df_wide$SNP
    
    # Pivot Sample Sizes (your existing code for this is correct)
    N_df_wide <- trait_df_list %>%
      select(trait, SNP, N_PH) %>%
      group_by(trait, SNP) %>% 
      dplyr::slice(1) %>% 
      ungroup() %>%
      pivot_wider(names_from = "trait", values_from = "N_PH")
    
    N_mat <- as.matrix(N_df_wide[, -1])
    rownames(N_mat) <- N_df_wide$SNP
    
    df_z <- as.data.frame(z_mat)
    df_N <- as.data.frame(N_mat)
    
  } else {
    print("No overlapping SNPs with matching alleles found in any trait files. Returning empty data frames.")
    df_z <- data.frame(matrix(ncol = length(trait_ss_files), nrow = 0, dimnames = list(NULL, names(trait_ss_files))))
    df_N <- data.frame(matrix(ncol = length(trait_ss_files), nrow = 0, dimnames = list(NULL, names(trait_ss_files))))
  }
  
  file.remove(tmp_var_file)
  print("Done!")
  return(list(df_z = df_z, df_N = df_N, df_gwas = df_input_wGWAS_filtered))
}





prep_z_matrix <- function(z_mat, N_mat,
                          corr_cutoff=0.85,
                          keep_my_traits=NULL,
                          rm_traits=NULL,
                          pval_cutoff=NULL,
                          nonneg=T) {
  
  # Given a matrix of z-scores (N_variants x M_traits) and vector of median
  # sample sizes per trait:
  # 1) perform final pre-processing steps before bNMF clustering:
  # trait filtering by p-value, trait pruning based on correlation,
  # and z-score scaling based on sample size
  # 2) expand N x M matrix into N x 2M non-negative matrix
  
  
  # if no user input for pval_cutoff, use bonferroni definition
  z_mat_orig <- z_mat
  
  if (is.null(pval_cutoff)){
    pval_cutoff <- 0.05 / nrow(z_mat) 
  }
  print(paste0(sum(is.na(z_mat)), " missing values before pval pruning."))
  
  df_traits <- data.frame(trait = colnames(z_mat))
  
  df_traits_filtered <- data.frame(trait=as.character(),
                           result=as.character(),
                           note=as.character())
  
  # 1.) remove manually entered traits
  if (!is.null(rm_traits)) {
    print(sprintf("Removing %i manually entered traits!!!", length(rm_traits)))
    z_mat <- z_mat[, !colnames(z_mat) %in% rm_traits]
    
    df_remove_manual <- data.frame(trait=rm_traits,
                                   result="removed (manual)",
                                   note=NA)
    df_traits_filtered <- rbind(df_traits_filtered, df_remove_manual)
  }
  
  # 2.) Filter traits by p-value (min. p-value < 0.05/N_variants)
  print("Filtering traits w/ no pvalues below cutoff...")
  minP_vec <- apply(z_mat, 2, function(x) min(2 * pnorm(abs(x), lower.tail=F), na.rm=T))
  
  traits_removed <- colnames(z_mat)[minP_vec >= pval_cutoff]
  traits_kept <- colnames(z_mat)[minP_vec < pval_cutoff]
  
  df_lowP_vec <- data.frame(minPval=minP_vec[minP_vec >= pval_cutoff]) %>%
    rownames_to_column('trait') %>%
    mutate(result="removed (p-value)") %>% 
    mutate(note=paste("min pval=",format(minPval,scientific=T))) %>%
    dplyr::select(trait, result, note)
  
  df_traits_filtered <- rbind(df_traits_filtered, df_lowP_vec)
  cat(paste(sprintf("Removing traits with no variant having p < %.3e",pval_cutoff),
            paste(traits_removed,
                  collapse="\n"),sep=":\n"))
  z_mat <- z_mat[, traits_kept]

  
  # 4.) Prune traits by correlation (remove traits with Pearson |r| > 0.85)
  cat(sprintf("\n\nPrune traits by correlation (remove traits with Pearson |r| > %.2f)\n",
              corr_cutoff))
  trait_cor_mat <- cor(z_mat, use="pairwise.complete.obs")  # Trait-trait correlation matrix
  write.table(trait_cor_mat,"./trait_cor_mat.txt", sep="\t")
  
  # sort by max(z) instead of min(pval)
  remaining_traits <- names(sort(apply(z_mat, 2, max, na.rm=T),decreasing = T))
  print(paste("Initial number of traits:",length(remaining_traits)))
  
  keep_traits <- c()
  df_corr_removed <- c()
  
  while (length(remaining_traits) > 0) {
    # append to list of traits to keep
    keep_traits <- c(keep_traits, remaining_traits[1])
    if (length(remaining_traits)==1) {
      break
    }
    trait_cor_mat <- trait_cor_mat[remaining_traits, remaining_traits]
    print(dim(trait_cor_mat))
    to_remove <- rownames(trait_cor_mat)[abs(trait_cor_mat[, remaining_traits[1]]) >= corr_cutoff]

    # track results in dataframe
    if (length(to_remove)>1) {
      df_corr_removed_tmp <- data.frame(trait=to_remove[to_remove!=remaining_traits[1]],
                                        result="removed (correlation)",
                                        note=paste("correlated w/", remaining_traits[1]))
      df_traits_filtered <- rbind(df_traits_filtered, df_corr_removed_tmp)
    }
    
    if (length(to_remove) > 1) {
      cat(paste(sprintf("Correlated trait being removed for %s:",remaining_traits[1]),
                paste(to_remove[!to_remove %in% remaining_traits[1]],collapse="\n"),"\n",sep="\n"))
    }
    
    # remove all correlated traits from remaining list
    remaining_traits <- setdiff(
      remaining_traits, 
      to_remove
      )
  }
  
  keep_traits <- c(keep_traits, keep_my_traits)
  z_mat <- z_mat_orig[, c(keep_traits)]
  
  df_traits_filtered <- df_traits_filtered %>%
    right_join(df_traits, by="trait") %>%
    mutate(result = ifelse(is.na(result), "trait kept", result))

  pruned_traits <- df_traits_filtered %>%
    filter(result!="trait kept") %>%
    pull(trait)
  
  cat(paste("Traits removed in pruning process:", 
              paste(pruned_traits, collapse="\n"),sep = "\n"))
  
  if (!is.null(keep_my_traits)) {
    cat(paste("Manually kept:", 
              paste(keep_my_traits, collapse="\n"),sep = "\n"))
    
  }

  cat(paste0("\nNumber of remaining traits: ",ncol(z_mat),"\n"))
  

  # Adjust z-scores by sample size for each variant-trait combo
  # i.e. (z = z / sqrt(medN) * mean(sqrt(medN_all_traits)))
  cat("\n\n")
  print("Performing sample size adjustment...")
  medN_vec <- apply(N_mat[, colnames(z_mat)], 2, median, na.rm=T)
  z_mat <- z_mat / sqrt(N_mat[, colnames(z_mat)]) * mean(sqrt(medN_vec))
  
  # Replace missing values with zero
  print("Replacing remaining missing values with zero...")
  print(paste0(sum(is.na(z_mat)), " missing values were replaced."))
  z_mat[is.na(z_mat)] <- 0
  
  print("Save matrix before splitting into non-negative...")
  data.frame(z_mat) %>%
    rownames_to_column('variant') %>%
    write_csv("scaled_filtered_zmat.csv")
  
  # Expand into N x 2M non-negative matrix
  if (nonneg==T) {
    print("Expanding z-score matrix into non-negative matrix (N-variants x 2M-traits)...")
    z_mat_pos <- z_mat
    z_mat_pos[z_mat_pos < 0] <- 0
    colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
    z_mat_neg <- -z_mat
    z_mat_neg[z_mat_neg < 0] <- 0
    colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
    final_z_mat <- cbind(z_mat_pos, z_mat_neg)
  } else {
    final_z_mat <- z_mat
  }
  
  return(list(final_z_mat = final_z_mat,
              df_traits   = df_traits_filtered))
}


fill_missing_zscores <- function(initial_zscore_matrices,
                                df_snps,
                                trait_ss_files,
                                trait_ss_size,
                                main_ss_filepath,
                                rsID_map_file,
                                tool_path=NULL,
                                method_proxy="TopLD",
                                method_fill="median",
                                population=NULL) {
  
  z_mat <- initial_zscore_matrices$z_mat
  N_mat <- initial_zscore_matrices$N_mat
  
  # # find rows and cols of missing data
  na_ix <- data.frame(which(is.na(z_mat), arr.ind=T))
  na_ix_N <- data.frame(which(is.na(N_mat), arr.ind=T))

  # find rows (SNPs) where data was removed
  need_proxies <- rownames(z_mat)[unique(na_ix$row)]
  # find cols (traits) where data was removed
  trait_ss_files_covers <- trait_ss_files[unique(na_ix$col)]

  my_cols <- c("VAR_ID","rsID","PVALUE","Population")
  need_covers_df <- df_snps %>%
    separate(VAR_ID, into=c("CHR","POS","REF","ALT"),sep="_",remove = F) %>%
    mutate(ChrPos = paste(CHR, POS,sep=":")) %>%
    filter(ChrPos %in% need_proxies) %>%
    dplyr::select(any_of(my_cols))

  if ("Population" %in% names(need_covers_df)) {
    print(paste("Populations needing covers:",
                paste(unique(need_covers_df$Population),collapse = "\n"),
                sep = "\n"))
    df_list <- list()
    i = 1
    for (cur_population in sort(unique(need_covers_df$Population))) {
      
      print(sprintf("Search for cover proxies for %s...", cur_population))
    
      need_covers_df_tmp <- need_covers_df %>%
        filter(Population==cur_population)
      
      # force use of LDlinkR if population now available in TopLD
      method_pop <- ifelse(cur_population %in% c("EUR","SAS","EAS","AFR"), method_proxy, "LDlinkR")
      
      cover_proxy_set_tmp <- choose_proxies(need_proxies = need_covers_df_tmp,
                                        tabix_path = tool_path,
                                        ld_file = NULL,
                                        rsID_map_file = rsID_map_file,
                                        trait_ss_files = trait_ss_files_covers,
                                        pruned_variants = df_snps,
                                        method = method_pop,
                                        population = cur_population,
                                        frac_nonmissing_num = 0.5,
                                        r2_num = 0.5
      )
      df_list[[i]] <- cover_proxy_set_tmp[[2]]
      i <- i+1
    }
    for (k in 1:length(df_list)) {
      print(head(df_list[[k]]))
    }
    df_covers <- do.call(rbind, df_list) %>%
      filter(!duplicated(proxy_VAR_ID))
    print(sprintf("%i total covers found for all populations...", nrow(df_covers)))
    
    } else {
  
    print(sprintf("Searching for cover proxies for %s...", population))
    cover_proxy_set <- choose_proxies(need_proxies = need_covers_df,
                                        tabix_path = tool_path,
                                        ld_file = NULL,
                                        rsID_map_file = rsID_map_file,
                                        trait_ss_files = trait_ss_files_covers,
                                        pruned_variants = df_snps,
                                        method = method_proxy,
                                        population=population,
                                        frac_nonmissing_num = 0.5,
                                        r2_num = 0.5
    )
    final_snps_orig <- cover_proxy_set[[1]]
    df_covers <- cover_proxy_set[[2]]
    }
  
  print(sprintf("%i cover proxies found", nrow(df_covers)))
  print(head(df_covers))
  
  df_covers_merged <- df_covers %>%
    data.frame() %>%
    inner_join(df_snps[,c("VAR_ID","GWAS","full_path")], by="VAR_ID") %>%
    mutate(full_path = full_path) #full_path

  final_snps_covers <- unique(df_covers_merged$proxy_VAR_ID)
  # print(sprintf("Unique cover proxies found: %i", length(final_snps_covers)))
  
  print("Is cover proxy search null?")
  print(is.null(final_snps_covers))
  
  if (!is.null(final_snps_covers)){
    
    print("Getting cover proxies risk alleles...")
    
    
    df_cover_alleles <- data.frame(VAR_ID=df_covers_merged$proxy_VAR_ID,
                                   Risk_Allele=NA)
    

    
    print(sprintf("Fetching summary stats for %i cover SNPs...",nrow(df_cover_alleles)))
    trait_ss_size_covers <- trait_ss_size[unique(na_ix$col)]
    
    initial_zscore_matrices_cover <- fetch_summary_stats(
      df_cover_alleles,
      main_ss_filepath,
      trait_ss_files_covers,
      trait_ss_size_covers,
      pval_cutoff=1
    )
    print("Finished getting summary stats!")
    
    z_mat_covers <- initial_zscore_matrices_cover$z_mat
    N_mat_covers <- initial_zscore_matrices_cover$N_mat
    
  
    # match ChrPos for orig SNPs and cover proxies
    df_covers_merged <- df_covers_merged %>%
      filter(proxy_VAR_ID %in% df_cover_alleles$VAR_ID) %>%
      mutate(ChrPos_orig = gsub("_", ":", str_before_nth(VAR_ID, "_", 2))) %>%
      mutate(ChrPos_proxy = gsub("_", ":", str_before_nth(proxy_VAR_ID, "_", 2)))
    
    
    print(paste("Num missing data before:",
                sum(is.na(z_mat))))  # 123
    ix = 0
    jx = 0
    kx = 0
    df_z_fixed <- z_mat
    df_N_fixed <- N_mat
    for (i in 1:nrow(na_ix)){
      cur_snp <- rownames(z_mat)[na_ix$row[i]]
      cur_trait <- colnames(z_mat)[na_ix$col[i]]
      
      if (cur_snp %in% df_covers_merged$ChrPos_orig) {
        cur_proxy <- df_covers_merged$ChrPos_proxy[df_covers_merged$ChrPos_orig==cur_snp]
        
        if (cur_proxy %in% rownames(z_mat_covers) & cur_trait %in% colnames(z_mat_covers)){
          new_zscore <- z_mat_covers[cur_proxy, cur_trait]
          new_N <- N_mat_covers[cur_proxy, cur_trait]
          
          if (!is.na(new_zscore)) { # we have valid value to use
            df_z_fixed[na_ix$row[i], na_ix$col[i]] <- new_zscore
            df_N_fixed[na_ix$row[i], na_ix$col[i]] <- new_N
            
            ix <- ix + 1
          } else { # we have a proxy, but not for this trait
          jx <- jx + 1 }
      } else {  # we dont have a proxy
        kx <- kx + 1 }
      }
    }
    print(paste("No. missing filled:", ix))
    print(paste("No. proxy found, but not for specific trait:", jx))
    print(paste("No. proxy not found for orig SNP:", kx))
    print(paste("Num missing data after:", sum(is.na(df_z_fixed))))  
    # sum(is.na(df_N_fixed))  
  }
  else {
    df_z_fixed <- z_mat
    df_N_fixed <- N_mat
  }
  
  #----
  
  # update z_matrix
  if (method_fill=="median"){
    print("Filling in remaining missing values with trait median...")
    
    z_mat_cover_med <- apply(df_z_fixed, 2, function(x)
      replace(x, is.na(x), median(x, na.rm = TRUE)))
  } else if (method_fill=="zero") {
    print("Filling in remaining missing values with zero...")
    
    z_mat_cover_med <- apply(df_z_fixed, 2, function(x)
      replace(x, is.na(x), 0))
  } else if (method_fill=="remove") {
    print("Removing remaining SNPs with missing values...")
    
    z_mat_cover_med <- na.omit(df_z_fixed)
    
  } else {
    print("Input 'median' or 'zero' for method_fill argument...")
  }
  sum(is.na(z_mat_cover_med))
  
  print("Matching N_mat row names with z_mat and replacing missing values with median...")
  df_N_fixed <- df_N_fixed[rownames(z_mat_cover_med), ]
  N_mat_cover_med <- apply(df_N_fixed, 2, function(x)
    replace(x, is.na(x), median(x, na.rm = TRUE)))

  list(z_mat=z_mat_cover_med, N_mat=N_mat_cover_med, proxy_df=df_covers_merged, z_mat_covers)
}


prep_z_matrix_wProteins <- function(z_mat, N_mat,
                                    corr_cutoff=0.85,
                                    keep_my_traits=NULL,
                                    rm_traits=NULL,
                                    pval_cutoff=NULL,
                                    nonneg=T,
                                    protein_cols=NULL,          
                                    filter_concentrated_proteins=TRUE,  
                                    concentration_threshold=75,  
                                    block_reweight=TRUE,        
                                    target_protein_ratio=1.2,
                                    z_cap=12) {  # Added z_cap argument
  
  # Standard P-value filtering setup
  z_mat_orig <- z_mat
  if (is.null(pval_cutoff)){
    pval_cutoff <- 0.05 / nrow(z_mat) 
  }
  
  df_traits_filtered <- data.frame(trait=as.character(), result=as.character(), note=as.character())
  
  # 1. Remove manual traits
  if (!is.null(rm_traits)) {
    z_mat <- z_mat[, !colnames(z_mat) %in% rm_traits]
    df_traits_filtered <- rbind(df_traits_filtered, 
                                data.frame(trait=rm_traits, result="removed (manual)", note=NA))
  }
  
  # ---------------------------------------------------------
  # NEW STEP 1: WINSORIZATION (Apply immediately)
  # ---------------------------------------------------------
  if (!is.null(z_cap)) {
    cat(sprintf("\n\nWinsorizing Z-scores at +/- %.2f\n", z_cap))
    n_capped <- sum(abs(z_mat) > z_cap, na.rm = TRUE)
    cat(sprintf("Capping %d values (%.2f%% of matrix)\n", 
                n_capped, 100 * n_capped / prod(dim(z_mat))))
    
    # Cap values while preserving sign
    z_mat[z_mat > z_cap] <- z_cap
    z_mat[z_mat < -z_cap] <- -z_cap
  }
  # ---------------------------------------------------------
  
  # 2. Filter by P-value
  minP_vec <- apply(z_mat, 2, function(x) min(2 * pnorm(abs(x), lower.tail=F), na.rm=T))
  traits_kept <- colnames(z_mat)[minP_vec < pval_cutoff]
  
  # Log removed traits
  removed_p <- setdiff(colnames(z_mat), traits_kept)
  if(length(removed_p) > 0) {
    df_traits_filtered <- rbind(df_traits_filtered, 
                                data.frame(trait=removed_p, result="removed (p-value)", note=paste0("<", pval_cutoff)))
  }
  z_mat <- z_mat[, traits_kept]
  
  # 3. Prune by Correlation
  cat(sprintf("\nPruning traits by correlation (|r| > %.2f)\n", corr_cutoff))
  trait_cor_mat <- cor(z_mat, use="pairwise.complete.obs")
  
  # Sort priorities by max Z (magnitude of signal)
  remaining_traits <- names(sort(apply(z_mat, 2, max, na.rm=T), decreasing = T))
  keep_traits <- c()
  
  while (length(remaining_traits) > 0) {
    curr <- remaining_traits[1]
    keep_traits <- c(keep_traits, curr)
    
    # Find correlates
    cor_vals <- trait_cor_mat[curr, remaining_traits]
    to_remove <- names(cor_vals)[abs(cor_vals) >= corr_cutoff & names(cor_vals) != curr]
    
    if (length(to_remove) > 0) {
      df_traits_filtered <- rbind(df_traits_filtered,
                                  data.frame(trait=to_remove, result="removed (correlation)", 
                                             note=paste("correlated w/", curr)))
    }
    remaining_traits <- setdiff(remaining_traits, c(curr, to_remove))
  }
  
  # Restore any manually kept traits
  final_keep <- unique(c(keep_traits, intersect(keep_my_traits, colnames(z_mat))))
  z_mat <- z_mat[, final_keep]
  
  # 4. Sample Size Adjustment
  cat("\nPerforming sample size adjustment...\n")
  medN_vec <- apply(N_mat[, colnames(z_mat)], 2, median, na.rm=T)
  z_mat <- z_mat / sqrt(N_mat[, colnames(z_mat)]) * mean(sqrt(medN_vec))
  z_mat[is.na(z_mat)] <- 0
  
  # ---------------------------------------------------------
  # NEW STEP 2: PROTEIN CONCENTRATION FILTER
  # ---------------------------------------------------------
  protein_filtering_results <- NULL
  if (!is.null(protein_cols) && filter_concentrated_proteins) {
    cat("\nChecking for hyper-concentrated proteins...\n")
    
    current_prots <- intersect(protein_cols, colnames(z_mat))
    
    if (length(current_prots) > 0) {
      # Calculate concentration metrics
      prot_metrics <- data.frame(protein = current_prots, top1_pct = NA_real_)
      
      for (i in 1:nrow(prot_metrics)) {
        p <- prot_metrics$protein[i]
        vals <- z_mat[, p]^2 # Use squared signal
        total_signal <- sum(vals)
        if (total_signal > 0) {
          prot_metrics$top1_pct[i] <- 100 * max(vals) / total_signal
        } else {
          prot_metrics$top1_pct[i] <- 0
        }
      }
      
      to_remove <- prot_metrics$protein[prot_metrics$top1_pct > concentration_threshold]
      
      if (length(to_remove) > 0) {
        cat(sprintf("Removing %d proteins with >%d%% signal in 1 SNP\n", length(to_remove), concentration_threshold))
        z_mat <- z_mat[, !colnames(z_mat) %in% to_remove]
        
        df_traits_filtered <- rbind(df_traits_filtered,
                                    data.frame(trait=to_remove, result="removed (concentration)", 
                                               note=paste0("top1 > ", concentration_threshold, "%")))
      }
    }
  }
  
  # 5. Expand to Non-Negative
  if (nonneg) {
    z_mat_pos <- z_mat; z_mat_pos[z_mat_pos < 0] <- 0
    colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
    
    z_mat_neg <- -z_mat; z_mat_neg[z_mat_neg < 0] <- 0
    colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
    
    final_z_mat <- cbind(z_mat_pos, z_mat_neg)
  } else {
    final_z_mat <- z_mat
  }
  
  # ---------------------------------------------------------
  # NEW STEP 3: BLOCK REWEIGHTING (Post-Expansion)
  # ---------------------------------------------------------
  block_reweighting_results <- NULL
  if (!is.null(protein_cols) && block_reweight) {
    
    # Identify columns (handling _pos/_neg suffixes)
    all_cols <- colnames(final_z_mat)
    is_prot <- grepl(paste(protein_cols, collapse="|"), all_cols) & 
      !grepl(paste(setdiff(colnames(z_mat_orig), protein_cols), collapse="|"), all_cols)
    
    prot_cols_final <- all_cols[is_prot]
    trait_cols_final <- all_cols[!is_prot]
    
    if (length(prot_cols_final) > 0 && length(trait_cols_final) > 0) {
      frob_p <- sqrt(sum(final_z_mat[, prot_cols_final]^2))
      frob_t <- sqrt(sum(final_z_mat[, trait_cols_final]^2))
      
      # Calculate scalar weight
      w <- (frob_t * target_protein_ratio) / frob_p
      
      cat(sprintf("\nReweighting Proteins: Traits=%.1f, Prots=%.1f -> Applying weight %.2f\n", 
                  frob_t, frob_p, w))
      
      final_z_mat[, prot_cols_final] <- final_z_mat[, prot_cols_final] * w
      
      block_reweighting_results <- list(weight=w, frob_trait=frob_t, frob_prot_orig=frob_p)
    }
  }
  
  return(list(final_z_mat=final_z_mat, 
              df_traits=df_traits_filtered,
              block_reweighting=block_reweighting_results))
}
