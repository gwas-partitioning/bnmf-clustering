# This script is intended to run the full pipeline for bNMF clustering based
# on summary statistics and test its agreement with scripts that are
# currently in place.

# NOTE: to run Tabix on UGER, must start with "ssh gsa4; use .zlib-1.2.6"

source("choose_variants.R")  # fld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("prep_bNMF.R")  # fetch_summary_stats & prep_z_matrix
source("run_bNMF.R")  # run_bNMF & summarize_bNMF

gwas_traits <- readxl::read_excel("../data/clustering_data_source.xlsx", sheet="gwas_traits")
trait_paths <- setNames(gwas_traits$full_path, gwas_traits$trait_name)
trait_paths <- trait_paths[!grepl("MAGIC", names(trait_paths))]  # Some MAGIC GWAS files don't have N_PH field

initial_t2d_snps <- read_tsv("../data/T2D_initial_vars_pval.txt")
set.seed(1)
initial_t2d_snps <- sample_n(initial_t2d_snps, size=1000) %>%
  select(VAR_ID=VAR_ID_hg19, PVALUE)

rsID_map_file <- "/humgen/diabetes2/users/clairekim/list_VARID_rsID_updated.txt"  # From dbSNP v1.38 -- maps positional IDs to rsIDs

# Variant choice steps

pruned_variants <- ld_pruning(initial_t2d_snps, rsID_map_file)

var_nonmissingness <- count_traits_per_variant(pruned_variants$VAR_ID, trait_paths)

proxies_needed_df <- find_variants_needing_proxies(pruned_variants, var_nonmissingness,
                                                   rsID_map_file)

tabix_path <- "/humgen/diabetes2/users/mvg/VariantClustering/tabix-0.2.6/tabix"
ld_file <- "/humgen/diabetes2/users/mvg/VariantClustering/LD_EUR.tsv.bgz"
final_variant_set <- choose_proxies(
  proxies_needed_df,
  tabix_path,
  ld_file,
  rsID_map_file,
  trait_paths,
  pruned_variants
)

# Prep bNMF steps

t2d_ss_filepath <- "/humgen/diabetes2/users/clairekim/Mahajan.NatGenet2018b.T2D.European_formatted.txt"
variant_vec <- final_variant_set
gwas_ss_file <- t2d_ss_filepath
trait_ss_files <- trait_paths
print("Retrieving z-scores and sample sizes for each trait...")
trait_df_long <- lapply(names(trait_ss_files), read_single_trait, variant_df) %>%
  bind_rows(.id="trait")  


initial_zscore_matrices <- fetch_summary_stats(
  final_variant_set,
  t2d_ss_filepath,
  trait_paths
)

final_zscore_matrix <- prep_z_matrix(initial_zscore_matrices$z_mat,
                                     initial_zscore_matrices$N_mat)

# Run bNMF steps

bnmf_reps <- run_bNMF(final_zscore_matrix, n_reps=10)

summarize_bNMF(bnmf_reps)



# # Run fetch_summary_stats function
# initial_trait_matrices <- fetch_summary_stats(trait_paths, t2d_snps)
# saveRDS(initial_trait_matrices, "../data/test_fetch_summary_stats_t2d_results.rds")
# 
# # Run prep_z_matrix function
# colnames(initial_trait_matrices$z_mat) <- make.names(colnames(initial_trait_matrices$z_mat))
# rownames(initial_trait_matrices$z_mat) <- gsub("^GWAS_", "", rownames(initial_trait_matrices$z_mat))
# colnames(initial_trait_matrices$z_mat) <- gsub("^GWAS_", "", colnames(initial_trait_matrices$z_mat))
# names(t2d_minP_vec) <- gsub("^GWAS_", "", names(t2d_minP_vec))
# rownames(initial_trait_matrices$N_mat) <- gsub("^GWAS_", "", rownames(initial_trait_matrices$N_mat))
# colnames(initial_trait_matrices$N_mat) <- gsub("^GWAS_", "", colnames(initial_trait_matrices$N_mat))
# colnames(initial_trait_matrices$N_mat) <- make.names(colnames(initial_trait_matrices$N_mat))
# base_t2d_p_df <- read_tsv("../data/T2D_v16_traits_pval.txt")
# t2d_minP_vec <- setNames(c(base_t2d_p_df$pval), base_t2d_p_df$trait)
# t2d_minP_vec["Lipoprotein A (quantile)_irnt_both"] <- 0  # Zeros were removed from LpA in Claire's pipeline
# names(t2d_minP_vec) <- gsub("-|\\,|\\(|\\)|\\ ", ".", names(t2d_minP_vec))
# names(t2d_minP_vec) <- gsub("(dv.)\\.", "\\1_", names(t2d_minP_vec))
# t2d_minP_vec <- t2d_minP_vec[!grepl("Creatinine|ystatin", names(t2d_minP_vec))]
# initial_trait_matrices$minP_vec <- t2d_minP_vec[colnames(initial_trait_matrices$z_mat)]
# bnmf_input_mat <- prep_z_matrix(initial_trait_matrices$z_mat,
#                                 initial_trait_matrices$N_mat,
#                                 initial_trait_matrices$minP_vec)


# ##### TEST IMMEDIATE PRE-PROCESSING AND bNMF PROCEDuRE #####
# 
# base_t2d_z_df <- read_tsv("../data/T2D_LDpruned_290snps_v16.txt", na=".")
# t2d_z_var_info_df <- select(base_t2d_z_df, 1:6)
# t2d_z_mat <- as.matrix(select(base_t2d_z_df, -(1:6)))
# rownames(t2d_z_mat) <- t2d_z_var_info_df$locus
# colnames(t2d_z_mat) <- gsub("\\.ZN", "", colnames(t2d_z_mat))
# keep_traits <- setdiff(colnames(t2d_z_mat), c("Creatinine..quantile._irnt_both",
#                                               "Cystatin.C..quantile._irnt_both"))
# t2d_z_mat <- t2d_z_mat[, keep_traits]
# 
# base_t2d_p_df <- read_tsv("../data/T2D_v16_traits_pval.txt")
# t2d_minP_vec <- setNames(c(base_t2d_p_df$pval), base_t2d_p_df$trait)
# t2d_minP_vec["Lipoprotein A (quantile)_irnt_both"] <- 0  # Zeros were removed from LpA in Claire's pipeline
# names(t2d_minP_vec) <- gsub("-|\\,|\\(|\\)|\\ ", ".", names(t2d_minP_vec))
# t2d_minP_vec <- t2d_minP_vec[keep_traits]
# 
# base_t2d_N_df <- read_tsv("../data/traits_samplesize_for_scaling_Jun2020.txt")
# t2d_medN_vec <- setNames(c(base_t2d_N_df$medN), base_t2d_N_df$trait)
# names(t2d_medN_vec) <- gsub("-|\\,|\\(|\\)|\\ ", ".", names(t2d_medN_vec))
# t2d_medN_vec <- t2d_medN_vec[keep_traits]
# 
# bnmf_input_mat <- prep_z_matrix(t2d_z_mat, t2d_minP_vec, t2d_medN_vec)
# saveRDS(bnmf_input_mat, "../data/bnmf_input_mat.rds")
# bnmf_input_mat <- bnmf_input_mat[, sort(colnames(bnmf_input_mat))]
# 
# set.seed(1)
# bnmf_res <- BayesNMF.L2EU(bnmf_input_mat, n.iter=1000)
# saveRDS(bnmf_res, "../data/bnmf_res.rds")
# 
# ##### CLAIRE RESULTS #####
# claire_bnmf_input_mat <- readRDS("../data/input.rds")
# rownames(claire_bnmf_input_mat) <- gsub("\\.ZN", "", rownames(claire_bnmf_input_mat))
# claire_bnmf_input_mat <- claire_bnmf_input_mat[sort(rownames(claire_bnmf_input_mat)), ]
# claire_bnmf_input_mat <- t(claire_bnmf_input_mat)
# 
# 
# 
# a <- read_tsv("../data/T2D_LDpruned_corrpruned_66traits_290snps_v16.2.txt")
# b <- select(a, -(1:4))
# colnames(b) <- gsub("\\.ZN", "", colnames(b))
# b <- b[, sort(colnames(b))]












# final_variant_ss <- read_tsv("../data/T2D_290snps_v16_ALLELES_BETA.txt") %>%
#   separate(VAR_ID_hg19, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove=F) %>%
#   select(VAR_ID=VAR_ID_hg19, REF, ALT, BETA=BETA_ALT)
# trait_filepath_df <- tibble(
#   root_path ="/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES",
#   prefix="GWAS",
#   study=c("MAGIC_dv1", "ADIPOGen_dv1", "BFpercent_dv1", "GiantUKBB_eu_dv1",
#           "MAGIC_dv1", "MAGIC_dv2", "MAGIC_eu_dv6", "GiantUKBB_eu_dv1",
#           "GIANT_eu_females_dv5", "GIANT_eu_males_dv5", "GIANT_eu_females_dv5",
#           "GIANT_eu_males_dv5", "MAGIC_dv1",  "MAGIC_dv1", "HRgene_dv1", "MAGIC_dv4",
#           "Leptin_eu_dv1", "Leptin_eu_dv1", "MAGIC_dv1", "GIANT_eu_females_dv5",
#           "GIANT_eu_males_dv5", "GIANT_eu_females_dv5", "GIANT_eu_males_dv5",
#           "GIANT_eu_females_dv5", "GIANT_eu_males_dv5", "GIANT_eu_females_dv5",
#           "GIANT_eu_males_dv5"),
#   trait=c("2hrG", "Adiponectin", "BFP", "BMI", "FG", "FI", "HBA1C", "HEIGHT", 
#           "HIPC", "HIPC", "HIPCadjBMI", "HIPCadjBMI", "HOMAB", "HOMAIR", "HR", 
#           "ISIadjAgeSexBMI", "LEP", "LEPadjBMI", "PI", "WAIST", "WAIST", 
#           "WAISTadjBMI", "WAISTadjBMI", "WHR", "WHR", "WHRadjBMI", "WHRadjBMI"),
#   dir="DATA"
# ) %>%
#   mutate(full_path=paste0(root_path, "/", prefix, "_", study, "/", trait, "/", dir, "/", 
#                           prefix, "_", study, ".", trait, ".1.txt"))
# trait_ss_list <- setNames(trait_filepath_df$full_path,
#                           paste(trait_filepath_df$study, trait_filepath_df$trait, sep="_"))

