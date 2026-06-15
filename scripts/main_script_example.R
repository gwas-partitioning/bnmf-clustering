# =============================================================================
# BMI bNMF Clustering Pipeline
# =============================================================================
# Description: Full pipeline for BMI GWAS clustering using Bayesian Non-negative
#              Matrix Factorization (bNMF). Identifies genetic subtypes of BMI
#              by clustering GWAS loci based on their associations across a panel
#              of metabolic traits.
#
# Usage: Adapt the USER CONFIGURATION section below for your environment, then
#        run sections sequentially. Steps 2-3 (LD pruning via LDlink) are
#        commented out by default since they require API access and may be
#        pre-computed; uncomment to run from scratch.
#
# Version history:
#  v1  - initial
#  v2  - LD pruning r2 cutoff from 0.1 → 0.05
#  v3  - updated BMI GWAS to multi-ancestry list
#  v4  - removed BMI from trait list
#  v5  - unisex traits
#  v6  - removed CIR, DI, Incr30
#  v7  - added OGTT & cardio traits
#  v8  - pipeline overhaul
#  v9  - updated pipeline & bNMF K0 setting
#  v10 - update WHR GWAS to UKBB-GIANT
#  v11 - update GWAS away from GIANT
#  v12 - use old Leptin GWAS as sensitivity analysis
#  v14 - fix proxy issue
#  v15 - update Leptin GWAS
#  v16 - actually fix proxy issue
#  v17 - update to match CAD-2
#  v18 - increase trait missingness filter to 30%
# =============================================================================

# =====================
# 0. Setup & Configuration
# =====================

library(data.table)
library(dplyr)
library(magrittr)
library(readxl)
library(softImpute)
library(strex)
library(tidyverse)
library(readr)

# ---- USER CONFIGURATION ----
# Set these variables before running the pipeline.

# Working directory (root of your project)
working_dir <- getwd()  # or set explicitly, e.g. "/path/to/project"

# Pipeline version label (used for output directory and checkpoint filenames)
version <- "my_trait_v1"

# Path to the GWAS manifest Excel file (see example_data/clustering_data_source_example.xlsx for the expected format)
gwas_file <- file.path(working_dir, "clustering_data_source_example.xlsx")

# Directory for helper R scripts (choose_variants, prep_bNMF, run_bNMF, post_bNMF)
scripts_dir <- working_dir

# LDlink API token — register at https://ldlink.nci.nih.gov/?tab=apiaccess
# Replace with your own token before running LD pruning steps.
my_LDlink_token <- Sys.getenv("LDLINK_TOKEN", unset = "YOUR_LDLINK_TOKEN_HERE")

# Variant selection thresholds
# PVCUTOFF:       genome-wide significance — defines sentinel variants and the LD-pruning input
# PVCUTOFF_PROXY: suggestive threshold — broader pool for proxy candidate fetch
#                 (set equal to PVCUTOFF to disable the wider pool)
# PROXY_WINDOW_KB: half-window around each sentinel (kb) for proxy candidate fetch;
#                  only variants within this window of a sentinel are fetched, keeping
#                  the z-score grid tractable even for highly polygenic traits (BMI, T2D)
PVCUTOFF          <- 5e-8
PVCUTOFF_PROXY    <- 5e-6
PROXY_WINDOW_KB   <- 500

# Directory containing per-chromosome rsID→position map files (used by choose_proxies)
# Each file should be named chr{N}.txt with columns: hg19_posID, rsID, ref_allele, alt_allele
my_rsid_map_dir <- "rsid_maps_by_chr"

# Column rename map for your primary GWAS summary stats file (if column names differ from
# the pipeline defaults: P_VALUE, BETA, SE, ODDS_RATIO). Set to NULL if no renaming needed.
# Example: rename_cols <- c(P_VALUE = "P_FIRTH_FE_IV", ODDS_RATIO = "OR_FIRTH_FE_IV")
rename_cols <- NULL

# Path to hg19→hg38 liftover chain file (required for post-hoc analysis and TOPMed audit)
# Download from: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
hg19_to_hg38_chain_file <- file.path(working_dir, "hg19ToHg38.over.chain")

# TOPMed presence filter (optional)
# When TRUE, variants not found in TOPMed are flagged for proxy replacement, and all
# proxy candidates are restricted to TOPMed-confirmed positions — ensuring the final
# variant set is fully genotyped in TOPMed (improves PRS portability).
# Requires: BRAVO per-chromosome VCF files (chr*.bravo.pub.vcf.gz + .tbi index).
# Download from: https://bravo.sph.umich.edu/freeze10/hg38/
USE_TOPMED_FILTER <- FALSE
bravo_dir         <- "/humgen/florezlab/users/ksmith/BRAVO"

# ---- END USER CONFIGURATION ----

setwd(working_dir)

# Derived output paths
main_dir      <- file.path(working_dir, paste0(version, "_results"))
main_dir_shrt <- basename(main_dir)
dir.create(main_dir, recursive = TRUE)

# Checkpoint file for saving R workspace between steps
df_save <- sprintf("my_workspace_%s.RData", version)

# Uncomment to resume from a saved checkpoint:
# load(df_save)


# =====================
# 1. Source Helper Scripts & Load GWAS Manifest
# =====================

source(file.path(scripts_dir, "choose_variants_2025.R"))  # LD pruning, clumping, proxy search
source(file.path(scripts_dir, "prep_bNMF_2025.R"))        # Summary stats fetch & z-matrix prep
source(file.path(scripts_dir, "run_bNMF_2025.R"))         # bNMF execution & summarization
source(file.path(scripts_dir, "post_bNMF_2025.R"))        # Post-hoc cluster analysis

select <- dplyr::select  # prevent namespace collision with other packages

# Load GWAS manifest

gwas <- read_excel(gwas_file, sheet = "main_gwas") %>%
  mutate(ID = paste(study, trait, population, sep = "_")) %>%
  drop_na(full_path)
# Note: gwas$full_path should point to summary statistic files on your system.

gwas_traits      <- read_excel(gwas_file, sheet = "trait_gwas") %>% drop_na(full_path)
main_ss_filepath <- gwas$full_path[gwas$largest == "Yes"]

# Named file vectors
gwas_ss_files      <- setNames(gwas$full_path, gwas$ID)
trait_ss_files     <- setNames(gwas_traits$full_path, gwas_traits$trait)
trait_ss_size      <- setNames(gwas_traits$sample_size, gwas_traits$trait)

# Sanity check: missing trait files
for (f in trait_ss_files) {
  if (!file.exists(f)) message("Missing trait file: ", f)
}

# =====================
# 2. Variant Selection & Clumping
# =====================
# NOTE: Steps 2-3 perform LD pruning via the LDlink API and may take several
# hours. Results are saved to disk via save.image() so you can resume from a
# checkpoint without re-running. If you already have pruned_vars, load the
# saved workspace and skip to Section 4.

start_time <- Sys.time()

# 2.1 Identify suggestive SNPs (PVCUTOFF_PROXY) as the broad proxy candidate pool.
#     The wider threshold gives the downstream proxy search more flexibility:
#     if a sentinel needs a proxy, any sub-threshold variant in its neighbourhood
#     is already fetched and evaluated.
vars_sig <- get_sig_snps(
  gwas         = gwas,
  rename_cols  = rename_cols,  # set in USER CONFIGURATION above; NULL if no renaming needed
  PVCUTOFF     = PVCUTOFF_PROXY    # broad pool; pruning step enforces PVCUTOFF
) %>%
  mutate(PVALUE = ifelse(PVALUE == 0, 1e-300, PVALUE))
message("Total suggestive SNPs (P < PVCUTOFF_PROXY): ", nrow(vars_sig))


# 2.2 Main-study variant set
tmp_main <- get_biggest_gwas(main_ss_filepath, vars_sig)
save.image(file = df_save)

# 2.3 Remove HLA region
vars_noHLA <- tmp_main %>%
  mutate(across(c(CHR, POS), as.integer)) %>%
  filter(!(CHR == 6 & between(POS, 28477797, 33448354)))
message("After HLA removal: ", nrow(vars_noHLA))
save.image(file = df_save)

# 2.4 SNP clumping — restrict to sentinel-level variants (PVCUTOFF) before clumping
#     so LDlinkR only processes genome-wide significant hits, not the full 5e-6 pool.
vars_noHLA_sentinel <- vars_noHLA %>% filter(PVALUE < PVCUTOFF)
message("Sentinel variants (P < PVCUTOFF) before clumping: ", nrow(vars_noHLA_sentinel))
clumped_ids <- snp_clump(vars_noHLA_sentinel, window = 100e3, id = 'VAR_ID')
vars_clumped <- vars_noHLA_sentinel %>% filter(VAR_ID %in% clumped_ids)
message("Variants after clumping: ", nrow(vars_clumped))
save.image(file = df_save)

# =====================
# 3. LD Pruning (LDlink SNPclip, multi-population)
# =====================

my_pops <- c("EUR","EAS","AFR","AMR","SAS")

for (LD_pop in my_pops) {
  print(cat(sprintf("\n\nLD-pruning using %s panel in LDlink::SNPclip...", LD_pop)))
  ld_pruning_SNP.clip(df_snps = vars_clumped,
                      pop = LD_pop,
                      output_dir = main_dir_shrt,
                      r2 = 0.05,
                      maf=0.001,
                      chr = 22:1,
                      token = my_LDlink_token)
}
print("Done!")
save.image(file = df_save)

# Combine LD-pruning results
num_pops <- length(my_pops)

print("Combining SNP.clip results...")
ld_files <- list.files(path=main_dir_shrt,
                       pattern = "^snpClip_results",
                       full.names = T)

my_ChrPos <- vars_clumped$ChrPos
print(sprintf("Starting with %i SNPs...", length(my_ChrPos)))

df_clipped_kept_all    <- data.frame()
df_clipped_removed_all <- data.frame()
df_missing_1000G       <- data.frame()
df_low_MAF_1000G       <- data.frame()

ld_cols <- c("RS_Number","Position","Alleles","Details","Population")
for (pop in my_pops) {
  print(pop)
  df_clipped_res <- data.frame(RS_Number=character(), Position=character(),
                                Alleles=character(), Details=character(),
                                Population=character())
  rename_cols_clipped <- c(RS_Number="RS Number")

  ld_files_tmp <- ld_files[grepl(pop, ld_files)]
  print(length(ld_files_tmp) == 22)

  for (ld_file in ld_files_tmp) {
    df <- fread(ld_file, stringsAsFactors = FALSE, data.table = FALSE) %>%
      dplyr::rename(any_of(rename_cols_clipped)) %>%
      mutate(Population = pop)
    if (!all(ld_cols %in% names(df)))
      print(sprintf("%s missing columns!", str_after_last(ld_file, "results_")))
    df_clipped_res <- rbind(df_clipped_res, df)
  }

  df_clipped_kept    <- df_clipped_res %>% filter(Details == "Variant kept.")
  df_clipped_removed <- df_clipped_res %>% filter(Details != "Variant kept.")
  print(sprintf("%i SNPs kept and %i removed for %s",
                nrow(df_clipped_kept), nrow(df_clipped_removed), pop))

  df_clipped_kept_all    <- rbind(df_clipped_kept_all,    df_clipped_kept)
  df_clipped_removed_all <- rbind(df_clipped_removed_all, df_clipped_removed)
}

# variant counts per reference panel
df_clipped_kept_all %>% count(Population)

# Keep only variants that passed pruning in ALL populations
num_independent_pops <- df_clipped_kept_all %>%
  count(Position) %>%
  group_by(n)

kept_in_all_pops <- num_independent_pops %>%
  filter(n == num_pops) %>%
  inner_join(df_clipped_kept_all[, c('RS_Number', 'Position')], by = 'Position') %>%
  distinct(RS_Number, .keep_all = TRUE)

pruned_vars <- vars_clumped %>%
  mutate(SNP = ChrPos, ChrPos = gsub("chr", "", ChrPos)) %>%
  inner_join(kept_in_all_pops, by = c('SNP' = 'Position'))

print(sprintf("Sig. SNPs pruned from %i to %i...", nrow(vars_clumped), nrow(pruned_vars)))
save.image(file = df_save)

# =====================
# 4. Summary Stats Fetch & Missingness
# =====================

# 4.1 Build proxy candidate pool: restrict vars_noHLA (P < PVCUTOFF_PROXY) to
#     variants within PROXY_WINDOW_KB of a LD-pruned sentinel. This keeps the
#     fetch tractable for highly polygenic traits (BMI, T2D) while ensuring any
#     proxy LDlinkR might return is already in the fetched grid.
fetch_input <- window_to_sentinels(
  candidates = vars_noHLA %>%
    separate(VAR_ID, into = c("CHR","POS","REF","ALT"), sep = "_", remove = FALSE) %>%
    mutate(ChrPos = paste(CHR, POS, sep = ":")) %>%
    arrange(PVALUE) %>%
    distinct(ChrPos, .keep_all = TRUE),
  sentinels  = pruned_vars,
  window_kb  = PROXY_WINDOW_KB
) %>%
  rename(SNP = ChrPos)
message("Proxy candidate pool (windowed fetch input): ", nrow(fetch_input), " variants")
#
# # 4.2 Single fetch across all traits — provides z-scores AND missingness for
# #     both sentinel variants and any proxy that might be needed.
# z_n_mats <- fetch_summary_stats(
#   df_input          = fetch_input,
#   gwas_ss_file      = main_ss_filepath,
#   trait_ss_files    = trait_ss_files,
#   trait_ss_size     = trait_ss_size,
#   pval_cutoff       = 0.05,
#   read_trait_method = 'datatable'
# )
# save.image(file = df_save)
# 
# # #----
#
# 
zmat0 <- z_n_mats$df_z
Nmat0 <- z_n_mats$df_N

# Filter traits based on sample size and missingness
med_N <- apply(Nmat0, 2, median, na.rm = TRUE)
big_N_traits <- names(med_N)[med_N > 5000]
zmat1 <- zmat0[, big_N_traits]

# FIND TRAITS MISSING >30% OF SNPs
traits_highMissing <- names(colSums(is.na(zmat1)))[colSums(is.na(zmat1)) > (nrow(zmat1) * 0.3)]
traits_lowMissing <- setdiff(big_N_traits, traits_highMissing)

df_zmat_fullset <- z_n_mats$df_z[, traits_lowMissing]
trait_ss_files_filtered <- trait_ss_files[traits_lowMissing]
trait_ss_size_filtered <- trait_ss_size[traits_lowMissing]

# Variant-level non-missingness calculation
var_counts <- df_zmat_fullset %>%
  rownames_to_column('ChrPos') %>%
  inner_join(pruned_vars[,c("ChrPos","VAR_ID")], by="ChrPos") %>%
  mutate(frac = rowSums(!is.na(dplyr::select(., all_of(traits_lowMissing)))) / length(traits_lowMissing))
nonmiss_frac <- setNames(var_counts$frac, var_counts$VAR_ID)

# QC CHECK 1: Verify all pruned variants got missingness calculated
message("QC Check 1: Missingness calculation coverage")
missing_frac_coverage <- length(nonmiss_frac) / nrow(pruned_vars)
message(sprintf("  ✓ Missingness calculated for %d/%d pruned variants (%.1f%%)",
                length(nonmiss_frac), nrow(pruned_vars), missing_frac_coverage * 100))

if (missing_frac_coverage < 0.95) {
  warning("Less than 95% of pruned variants have missingness data!")
}

save.image(file = df_save)

# =====================
# 5. TOPMed Audit (optional)
# =====================
# Checks pruned_vars and the proxy candidate pool against TOPMed BRAVO VCFs.
# Requires USE_TOPMED_FILTER = TRUE and bravo_dir pointing to chr*.bravo.pub.vcf.gz files.
# Produces:
#   topmed_fails        — VAR_IDs in pruned_vars absent from TOPMed (fed to find_variants_needing_proxies)
#   topmed_present_snps — hg19 ChrPos values in the proxy candidate pool confirmed in TOPMed
#                         (fed to choose_proxies to restrict all proxy selection to TOPMed variants)

if (USE_TOPMED_FILTER) {
  message("=== TOPMED AUDIT ===")

  # 5a. Check pruned_vars against TOPMed
  topmed_confirmed_pruned <- check_topmed_presence(
    variants_hg19 = pruned_vars %>%
      mutate(CHR = as.character(CHR), POS = as.integer(POS)) %>%
      select(VAR_ID, CHR, POS),
    chain_file    = hg19_to_hg38_chain_file,
    vcf_dir       = bravo_dir
  )
  topmed_fails <- pruned_vars$VAR_ID[!pruned_vars$VAR_ID %in% topmed_confirmed_pruned]
  message(sprintf("  %d / %d pruned variants absent from TOPMed (will be flagged for proxy search)",
                  length(topmed_fails), nrow(pruned_vars)))

  # 5b. Check proxy candidate pool (fetch_input) against TOPMed
  #     Only variants confirmed here are eligible as proxies for any variant
  message("  Checking proxy candidate pool against TOPMed...")
  proxy_candidates_hg19 <- fetch_input %>%
    mutate(CHR = as.character(CHR), POS = as.integer(POS)) %>%
    select(VAR_ID, CHR, POS)

  topmed_confirmed_proxies <- check_topmed_presence(
    variants_hg19 = proxy_candidates_hg19,
    chain_file    = hg19_to_hg38_chain_file,
    vcf_dir       = bravo_dir
  )
  # Convert confirmed VAR_IDs back to ChrPos (zmat_fullset rowname format)
  topmed_present_snps <- fetch_input$ChrPos[fetch_input$VAR_ID %in% topmed_confirmed_proxies]
  message(sprintf("  %d / %d proxy candidates confirmed in TOPMed",
                  length(topmed_present_snps), nrow(fetch_input)))

} else {
  topmed_fails        <- NULL
  topmed_present_snps <- NULL
}

# =====================
# 6. Identify Variants Needing Proxies
# =====================

message("=== IDENTIFYING VARIANTS NEEDING PROXIES ===")

proxies_needed <- find_variants_needing_proxies(
  gwas_variant_df    = pruned_vars,
  var_nonmissingness = nonmiss_frac,
  topmed_fails       = topmed_fails   # NULL if USE_TOPMED_FILTER = FALSE
)

# QC CHECK 2: Verify all high-missingness variants are captured
message("QC Check 2: High-missingness variant detection")
high_missing_vars <- names(nonmiss_frac[nonmiss_frac < 0.8])
missed_high_missing <- high_missing_vars[!high_missing_vars %in% proxies_needed$VAR_ID]

if (length(missed_high_missing) == 0) {
  message("  ✓ All high-missingness variants correctly identified for proxy search")
} else {
  warning(sprintf("  ✗ %d high-missingness variants NOT in proxies_needed:", length(missed_high_missing)))
  print(head(missed_high_missing))
}

# QC CHECK 3: Verify no problematic variants in non-proxy set
message("QC Check 3: Problematic variant exclusion")
non_proxy_vars <- pruned_vars %>% filter(!VAR_ID %in% proxies_needed$VAR_ID)

# Check for ambiguous variants not flagged
ambiguous_missed <- non_proxy_vars %>%
  filter(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")) %>%
  nrow()

# Check for multiallelic variants not flagged
multiallelic_missed <- non_proxy_vars %>%
  filter(grepl(",", ALT)) %>%
  nrow()

if (ambiguous_missed == 0 && multiallelic_missed == 0) {
  message("  ✓ No ambiguous or multiallelic variants in non-proxy set")
} else {
  warning(sprintf("  ✗ Found %d ambiguous and %d multiallelic variants not flagged for proxy search",
                  ambiguous_missed, multiallelic_missed))
}

save.image(file = df_save)

# =====================
# 7. Proxy Search
# =====================

message("=== PROXY SEARCH ===")

all_need_proxies <- proxies_needed %>%
  inner_join(pruned_vars[,c("VAR_ID","Population")], by="VAR_ID")

final_proxy_results <- choose_proxies(
  need_proxies        = all_need_proxies,
  rsid_map_dir        = my_rsid_map_dir,
  zmat_fullset        = df_zmat_fullset,
  pruned_variants     = pruned_vars,
  token               = my_LDlink_token,
  population          = "EUR",
  frac_nonmissing_num = 0.8,
  r2_num              = 0.8,
  topmed_present_snps = topmed_present_snps   # NULL if USE_TOPMED_FILTER = FALSE
)

save.image(file = df_save)

# =====================
# 8. Final Variant Assembly
# =====================
proxy_not_needed_ids <- final_proxy_results[[1]]

# final_proxy_results[[2]] is NULL when no proxies are found; produce empty frame in that case
proxy_df <- if (!is.null(final_proxy_results[[2]]) && nrow(final_proxy_results[[2]]) > 0) {
  final_proxy_results[[2]] %>%
    select(ChrPos = proxy_ChrPos, rsID = proxy_rsID, A1 = REF, A2 = ALT, original_SNP = VAR_ID) %>%
    inner_join(fetch_input, by = c('ChrPos' = 'SNP')) %>%
    filter((A1 == REF & A2 == ALT) | (A2 == REF & A1 == ALT)) %>%
    mutate(Variant_Type = "proxy") %>%
    select(VAR_ID, ChrPos, rsID, Variant_Type, CHR, POS, REF, ALT, PVALUE, Risk_Allele, GWAS, Population, original_SNP)
} else {
  message("  No proxies found; proxy_df will be empty.")
  data.frame(VAR_ID=character(), ChrPos=character(), rsID=character(), Variant_Type=character(),
             CHR=character(), POS=character(), REF=character(), ALT=character(),
             PVALUE=numeric(), Risk_Allele=character(), GWAS=character(),
             Population=character(), original_SNP=character())
}

# Combine final variant set
df_final <- pruned_vars %>%
  filter(!VAR_ID %in% all_need_proxies$VAR_ID) %>%
  mutate(Variant_Type="original", original_SNP=VAR_ID) %>%
  select(VAR_ID, ChrPos, rsID=RS_Number, Variant_Type, CHR, POS, REF, ALT, PVALUE, Risk_Allele, GWAS, Population, original_SNP) %>%
  rbind(proxy_df)

# Create and save rsID mapping
df_rsIDs <- df_final %>%
  select(VAR_ID, rsID)

write_delim(df_rsIDs, file.path(main_dir_shrt, "rsID_map.txt"), delim = "\t", col_names = TRUE)
message(sprintf("  Written rsID map with %d variants to rsID_map.txt", nrow(df_rsIDs)))



# QC CHECK 4: Verify variant accounting
message("QC Check 4: Variant accounting")
original_count <- nrow(pruned_vars)
proxies_needed_count <- nrow(proxies_needed)
proxies_found_count <- nrow(proxy_df)
proxies_not_found <- proxies_needed_count - proxies_found_count
final_count <- nrow(df_final)
expected_final <- original_count - proxies_not_found

message(sprintf("  Started with: %d pruned variants", original_count))
message(sprintf("  Needed proxies: %d variants", proxies_needed_count))
message(sprintf("  Found proxies: %d variants", proxies_found_count))
message(sprintf("  Expected final: %d variants", expected_final))
message(sprintf("  Actual final: %d variants", final_count))

if (final_count == expected_final) {
  message("  ✓ Variant accounting correct")
} else {
  warning(sprintf("  ✗ Variant accounting mismatch! Difference: %d", final_count - expected_final))
}

# =====================
# 9. GWAS Alignment & Beta Check
# =====================

message("=== GWAS ALIGNMENT ===")

variants_that_got_proxies <- proxy_df %>% distinct(original_SNP) %>% pull(original_SNP)
variants_needing_proxies_no_solution <- all_need_proxies %>%
  filter(!VAR_ID %in% variants_that_got_proxies)

gwas_final_filtered <- z_n_mats$df_gwas %>%
  select(SNP, Risk_Allele, Nonrisk_Allele, P_VALUE, BETA, SE) %>%
  inner_join(df_final, by = c('SNP'='ChrPos')) %>%
  mutate(BETA_aligned = case_when(
    ALT == Risk_Allele.x ~ BETA,
    REF == Risk_Allele.x ~ -BETA,
    TRUE ~ NA_real_
  )) %>%
  select(VAR_ID, ChrPos=SNP, rsID, REF, ALT, Risk_Allele=Risk_Allele.x, P_VALUE, SE, BETA, BETA_aligned)

write_csv(gwas_final_filtered, file.path(main_dir_shrt, "alignment_GWAS_summStats.csv"))
message(sprintf("  Written aligned GWAS summary stats with %d variants to alignment_GWAS_summStats.csv", nrow(gwas_final_filtered)))


# QC CHECK 5: Beta alignment verification
message("QC Check 5: Beta alignment")
positive_betas <- all(gwas_final_filtered$BETA_aligned > 0, na.rm = TRUE)
na_betas <- sum(is.na(gwas_final_filtered$BETA_aligned))

if (positive_betas && na_betas == 0) {
  message("  ✓ All BETA_aligned values are positive")
} else {
  warning(sprintf("  ✗ Beta alignment issues: %d NA values, all positive = %s",
                  na_betas, positive_betas))
}

# QC CHECK 6: Final problematic variant verification
message("QC Check 6: Final dataset quality")
final_ambiguous <- gwas_final_filtered %>%
  filter(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")) %>%
  nrow()

final_multiallelic <- gwas_final_filtered %>%
  filter(grepl(",", ALT)) %>%
  nrow()

if (final_ambiguous == 0 && final_multiallelic == 0) {
  message("  ✓ No ambiguous or multiallelic variants in final dataset")
} else {
  warning(sprintf("  ✗ Final dataset contains %d ambiguous and %d multiallelic variants",
                  final_ambiguous, final_multiallelic))
}

# =====================
# 10. Generate QC Summary Report
# =====================

qc_summary <- list(
  timestamp = Sys.time(),
  version = version,

  # Variant counts
  pruned_variants = nrow(pruned_vars),
  proxies_needed = nrow(proxies_needed),
  proxies_found = nrow(proxy_df),
  final_variants = nrow(gwas_final_filtered),

  # QC results
  missingness_coverage = missing_frac_coverage,
  missed_high_missing = length(missed_high_missing),
  ambiguous_missed = ambiguous_missed,
  multiallelic_missed = multiallelic_missed,
  accounting_correct = (final_count == expected_final),
  all_betas_positive = positive_betas,
  na_betas = na_betas,
  final_ambiguous = final_ambiguous,
  final_multiallelic = final_multiallelic,

  # Trait filtering
  total_traits = length(trait_ss_files),
  big_N_traits = length(big_N_traits),
  traits_high_missing = length(traits_highMissing),
  traits_final = length(traits_lowMissing)
)

saveRDS(qc_summary, file.path(main_dir, "qc_summary.rds"))

# Write human-readable QC report
qc_report_file <- file.path(main_dir, "quality_control_report.txt")
cat("=== QUALITY CONTROL REPORT ===\n", file = qc_report_file)
cat(sprintf("Generated: %s\n", Sys.time()), file = qc_report_file, append = TRUE)
cat(sprintf("Version: %s\n\n", version), file = qc_report_file, append = TRUE)

cat("VARIANT PROCESSING:\n", file = qc_report_file, append = TRUE)
cat(sprintf("  Pruned variants: %d\n", qc_summary$pruned_variants), file = qc_report_file, append = TRUE)
cat(sprintf("  Proxies needed: %d\n", qc_summary$proxies_needed), file = qc_report_file, append = TRUE)
cat(sprintf("  Proxies found: %d\n", qc_summary$proxies_found), file = qc_report_file, append = TRUE)
cat(sprintf("  Final variants: %d\n\n", qc_summary$final_variants), file = qc_report_file, append = TRUE)

cat("QUALITY CHECKS:\n", file = qc_report_file, append = TRUE)
cat(sprintf("  ✓ Missingness coverage: %.1f%%\n", qc_summary$missingness_coverage * 100), file = qc_report_file, append = TRUE)
cat(sprintf("  ✓ All high-missing captured: %s\n", ifelse(qc_summary$missed_high_missing == 0, "PASS", "FAIL")), file = qc_report_file, append = TRUE)
cat(sprintf("  ✓ Variant accounting correct: %s\n", ifelse(qc_summary$accounting_correct, "PASS", "FAIL")), file = qc_report_file, append = TRUE)
cat(sprintf("  ✓ All BETAs positive: %s\n", ifelse(qc_summary$all_betas_positive, "PASS", "FAIL")), file = qc_report_file, append = TRUE)
cat(sprintf("  ✓ No problematic variants in final set: %s\n", ifelse(qc_summary$final_ambiguous == 0 && qc_summary$final_multiallelic == 0, "PASS", "FAIL")), file = qc_report_file, append = TRUE)

message(sprintf("QC report written to: %s", qc_report_file))


# =====================
# 11. Imputation & Cross-Validation
# =====================

print("Imputing missing variant-trait data...")

# Simple approach - use gwas_final_filtered ChrPos to subset the matrices
zmat_preImputed <- z_n_mats$df_z[gwas_final_filtered$ChrPos, traits_lowMissing]
Nmat_preImputed <- z_n_mats$df_N[gwas_final_filtered$ChrPos, traits_lowMissing]

obs_mask <- !is.na(zmat_preImputed)
set.seed(123)

# Create CV holdouts
i_dx <- which(obs_mask, arr.ind = TRUE)
holdout_n <- floor(0.1 * nrow(i_dx))
cv_idx <- i_dx[sample(nrow(i_dx), holdout_n), ]

zmat_cv <- zmat_preImputed
zmat_cv[cv_idx] <- NA

# Lambda grid search
lam0 <- lambda0(zmat_preImputed)
lambda_grid <- seq(0.1 * lam0, lam0 + 0.2, length.out = 10)
errors <- sapply(lambda_grid, function(l) {
  fit <- softImpute(zmat_cv, rank.max = 50, lambda = l, type = 'svd')
  imp <- softImpute::complete(zmat_cv, fit)
  sqrt(mean((zmat_preImputed[cv_idx] - imp[cv_idx])^2, na.rm = TRUE))
})
best_lambda <- lambda_grid[which.min(errors)]
message('Best lambda:', best_lambda)

# Final imputation fits
final_fit <- softImpute(zmat_preImputed, rank.max = 50, lambda = best_lambda, type = 'svd')
z_imputed <- softImpute::complete(zmat_preImputed, final_fit)
N_imputed <- apply(Nmat_preImputed, 2, function(x) replace(x, is.na(x), median(x, na.rm = TRUE)))

save.image(file = df_save)

# =====================
# 12. Visualization of Imputation Results
# =====================

print("Visualizing softImpute results...")
dist_long <- as.data.frame(z_imputed) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = 'trait', values_to = 'value')
was_imp   <- as.data.frame(is.na(zmat_preImputed)) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = 'trait', values_to = 'imputed')
dist_data <- left_join(dist_long, was_imp, by = c('row_id','trait'))

dist_data %>%
  ggplot(aes(x = value, color = imputed)) +
  geom_density() +
  labs(
    x     = "Z-score",
    color = "Imputed?",
    title = "Density of Observed vs. Imputed Z-scores"
  ) +
  theme_minimal() +
  xlim(-5, 5)   # zoom in on the bulk of the distribution

# Build a little data.frame of CV holdouts
cv_df <- data.frame(
  true = zmat_preImputed[cv_idx],
  imp  = z_imputed[cv_idx]
)
ggplot(cv_df, aes(x = true, y = imp)) +
  geom_hex(bins = 50) +
  geom_abline(linetype = "dashed") +
  labs(x = "True Z", y = "Imputed Z",
       title = "SoftImpute: True vs. Imputed on CV Holdouts") +
  theme_minimal()

# =====================
# 13. bNMF Clustering
# =====================

# Prepare final z-matrix for bNMF
prepped <- prep_z_matrix(
  z_mat        = z_imputed,
  N_mat        = N_imputed,
  corr_cutoff  = 0.8
)
final_zscore_matrix <- as.matrix(prepped$final_z_mat)
df_traits_filtered     <- prepped$df_traits

# Optionally append traits that were removed upstream (e.g. low N) to the log.
# Set remove_traits <- c("trait1", "trait2") before this block if needed.
if (exists("remove_traits") && length(remove_traits) > 0) {
  df_traits_filtered <- rbind(df_traits_filtered,
                              data.frame(trait=remove_traits,
                                         result="removed (low N)",
                                         note=NA))
}

if (length(traits_highMissing)>0) {
  df_traits_filtered <- rbind(df_traits_filtered,
                              data.frame(trait=traits_highMissing,
                                         result="removed (high missing)",
                                         note=NA))
}
write_csv(x = df_traits_filtered,
          file = file.path(main_dir_shrt,"df_traits.csv"))


message(sprintf("Final matrix dimensions: %d SNPs x %d traits",
                nrow(final_zscore_matrix), ncol(final_zscore_matrix)/2))

#-----


# Run bNMF (parallel)
bnmf_start <- Sys.time()
my_n_reps <- 100
my_K <- 15
my_K0 <- 10
my_tolerance <- 1e-6
my_phi <- 1

my_bNMF_settings <- list(n_reps=my_n_reps,
                         K=my_K,
                         K0=my_K0,
                         tolerance=my_tolerance,
                         phi=my_phi)

saveRDS(my_bNMF_settings,file.path(main_dir, "bnmf_settings.rds"))


library(furrr)
plan(multisession, workers = pmax(1, parallelly::availableCores() - 2))

bnmf_out    <- run_bNMF_parallel(
  final_zscore_matrix,
  n_reps    = my_n_reps,
  K         = my_K,
  K0        = my_K0,
  tolerance = my_tolerance,
  phi       = my_phi,
  random_seed = 123
)
bnmf_end <- Sys.time()
message('bNMF runtime:', difftime(bnmf_end, bnmf_start, units = 'mins'))

# Summarize results
summarize_bNMF(bnmf_out, dir_save = main_dir)
plan(sequential)  # reset parallel backend

save.image(file = df_save)

# ====================


hg19_to_hg38_chain <- rtracklayer::import.chain(hg19_to_hg38_chain_file)

cluster_post_hoc <- do_post_analysis(main_dir = main_dir,
                                     df_gwas = gwas_final_filtered,
                                     my_chain = hg19_to_hg38_chain,
                                     make_excel = T,
                                     do_liftover = T,
                                     do_v2g = T)

save.image(file = sprintf("my_workspace_%s.RData", version))

saveRDS(cluster_post_hoc,file.path(main_dir, "cluster_post_hoc.rds"))
saveRDS(zmat_preImputed,file.path(main_dir, "zmat_preImputed.rds"))
saveRDS(z_imputed,file.path(main_dir, "zmat_imputed.rds"))



#---  HTML ----


k <- NULL
if (is.null(k)){
  html_filename <- paste0(version,"_maxK.html")
} else {
  html_filename <- sprintf("%s_K%i.html",version,k)
}


rmarkdown::render(
  file.path(scripts_dir, "format_bNMF_results_2025.Rmd"),  # update filename as needed
  output_file = html_filename,
  output_dir = main_dir,
  params = list(main_dir = main_dir)
)


#-----  trait table -----

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)   # for string matching

# Step 1: Put all your trait names into a vector
traits <- gwas_traits$file_name

# Create a data frame
df <- data.frame(trait = traits, stringsAsFactors = FALSE) %>%
  mutate(
    category = case_when(
      str_detect(trait, regex("BMI|HEIGHT|WHR|WAIST|HIP|WHRadj|HIPCadj|WAISTadj", ignore_case = TRUE)) ~ "Anthropometry / General Obesity",
      
      str_detect(trait, regex("vat|asat|gfat|fat.*(percentage|trunk)|body fat|trunk fat|adiponectin|leptin", ignore_case = TRUE)) ~ "Body Composition / Fat Distribution",
      
      str_detect(trait, regex("ISI|IFC|HOMA|HOMAIR|FG|FI|HbA1c|Glucose|2hGlu|Incr30|Ins30|DI|CIR|STUM|INS|BIG|Proins|2hr", ignore_case = TRUE)) ~ "Glycemic / Insulin Resistance",
      
      str_detect(trait, regex("HDL|LDL|Cholesterol|TG|Apolipoprotein|Lipoprotein|TC", ignore_case = TRUE)) ~ "Lipids / Lipoproteins",
      
      str_detect(trait, regex("ALT|AST|ALP|GGT|Gamma.*glutamyl|bilirubin|crp|reactive", ignore_case = TRUE)) ~ "Liver Enzymes / Function",
      
      str_detect(trait, regex("RBC|WBC|platelet|neutrophil|lymphocyte|monocyte|basophil|eosinophil|reticulocyte|MCV| corpuscular|sphered|rheum|red|count|creat|cyst|haem|heart", ignore_case = TRUE)) ~ "Hematology / Blood Counts",
      
      str_detect(trait, regex("SBP|DBP|HR|blood pressure|vent|duration|pulse", ignore_case = TRUE)) ~ "Blood Pressure / Cardiovascular",
      
      str_detect(trait, regex("testosterone|oestradiol|SHBG|IGF|vitamin D|leptin|adiponectin", ignore_case = TRUE)) ~ "Hormones / Adipokines",
      
      str_detect(trait, regex("CRP|albumin|calcium|phosphate|urate|urea|total protein", ignore_case = TRUE)) ~ "Inflammation / Other Biomarkers",
      
      str_detect(trait, regex("basal metabolic|metabolic rate", ignore_case = TRUE)) ~ "Metabolic Rate",
      
      TRUE ~ "Other / Unclassified"
    ),
    
    # Optional: extract base trait name without _female/_male/adjBMI for cleaner display
    base_trait = str_remove(trait, "_female|_male|adjBMI|adjbmi"),
    display_name = ifelse(str_detect(trait, "_female"), paste(base_trait, "(F)"),
                          ifelse(str_detect(trait, "_male"), paste(base_trait, "(M)"), trait)),
    display_name = str_wrap(gsub("_", " ", display_name), width = 20)
  ) %>%
  arrange(category)

# Step 3: Create a grid-like colored table with ggplot (tile style)
# Count per category to make reasonable number of columns
n_per_row <- 8   # adjust to your liking

df_plot <- df %>%
  mutate(
    category = factor(category),              # for ordering
    row_id = row_number(),
    col_id = (row_id - 1) %% n_per_row + 1,
    row_group = (row_id - 1) %/% n_per_row + 1
  )

ggplot(df_plot, aes(x = col_id, y = -row_group, fill = category, label = display_name)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(size = 5, color = "black", hjust = 0.5, vjust = 0.5) +
  scale_fill_brewer(palette = "Set3", name = "Category") +   # or use viridis::scale_fill_viridis(discrete = TRUE)
  theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(title = "Traits Grid Colored by Category",
       subtitle = "Each tile = one trait (sex-specific shown as F/M)") +
  coord_fixed(ratio = 0.5)   # adjust ratio so text fits nicely

# Optional: save as high-res image
ggsave("traits_grid_colored.png", width = 20, height = 16, dpi = 300)


