# This script is intended to run the full pipeline for bNMF clustering based
# on summary statistics and test its agreement with scripts that are
# currently in place.

# NOTE: to run Tabix on UGER, must start with "ssh gsa4; use .zlib-1.2.6"

source("choose_variants.R")  # fld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("prep_bNMF.R")  # fetch_summary_stats & prep_z_matrix
source("run_bNMF.R")  # run_bNMF & summarize_bNMF

gwas_traits <- readxl::read_excel("../data/clustering_data_source.xlsx", sheet="gwas_traits")
trait_ss_files <- setNames(gwas_traits$full_path, gwas_traits$trait_name)
trait_ss_files <- trait_ss_files[!grepl("MAGIC", names(trait_ss_files))]  # Some MAGIC GWAS files don't have N_PH field

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
  trait_ss_files,
  pruned_variants
)

# Prep bNMF steps

t2d_ss_filepath <- "/humgen/diabetes2/users/clairekim/Mahajan.NatGenet2018b.T2D.European_formatted.txt"

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
