# This script shows how to run the bNMF clustering pipeline using toy datasets.

#----
start=Sys.time()

# load requires packages
install.packages("pacman")
pacman::p_load(tidyverse, data.table, readxl, magrittr, dplyr, strex,
               rstudioapi, DT, kableExtra, GenomicRanges)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("Homo.sapiens")

# load project scripts containing bNMF functions
source("../scripts/choose_variants.R")  # ld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("../scripts/prep_bNMF.R")  # fetch_summary_stats & prep_z_matrix
source("../scripts/run_bNMF.R")  # run_bNMF & summarize_bNMF

setwd(dirname(getActiveDocumentContext()$path))

#----

# USER INPUTS!!!
project_dir = './test_results' # path to where you want results saved
user_token = 'cb5457b210a6' # token for LDlinkR api

# create project folder 
dir.create(project_dir)

#----

# SECTION 1: PULL IN GWAS INFORMATION

data_dir = "../example_data/"
rsID_map_file <- file.path(data_dir, "rsID_map_example.txt")  # From dbSNP v1.38 -- maps positional IDs to rsIDs

# GWAS for main trait
gwas <- read_excel(file.path(data_dir, "clustering_data_sources_example.xlsx"),
                   sheet="main_gwas") %>%
  data.frame()

# GWAS for clustering traits
gwas_traits <- read_excel(file.path(data_dir, "clustering_data_sources_example.xlsx"),
                          sheet="trait_gwas")

# GWAS to be used for final allele alignment
main_ss_filepath <- gwas %>% filter(largest=="Yes") %>% pull(full_path)

gwas_ss_files <- setNames(gwas$full_path, gwas$study)
trait_ss_files <- setNames(gwas_traits$full_path, gwas_traits$trait)
trait_ss_size <- setNames(gwas_traits$sample_size, gwas_traits$trait)

#----

# SECTION 2: PULL SIGNIFICANT VARIANTS FROM MAIN TRAIT GWAS

# P-value threshold for variants in main trait
PVCUTOFF = 5e-8

n_gwas <- length(gwas_ss_files)

vars_sig = data.frame(VAR_ID = as.character(),
                      P_VALUE = as.numeric(),
                      Risk_Allele=as.character(),
                      GWAS=as.character())

print(sprintf("Pulling significant SNPs w/ pval<%.1e from %i T2D GWAS...", PVCUTOFF, n_gwas))

for(i in 1:n_gwas) {
  print(paste0("...Reading ", names(gwas_ss_files)[i], "..."))

  vars <- fread(gwas_ss_files[i], data.table = F, stringsAsFactors=F)

  if (!"BETA" %in% colnames(vars)){
    print("Converting Odds Ratio to Log Odds Ratio...")
    vars <- vars %>%
      mutate(BETA = log(as.numeric(ODDS_RATIO)))
  }
  vars <- vars %>%
    filter(as.numeric(P_VALUE) <= PVCUTOFF) %>%
    subset(grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]", VAR_ID)) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove = F) %>%
    mutate(Risk_Allele = ifelse(BETA>=0, ALT, REF)) %>%
    mutate(GWAS = gwas$study[i]) %>%
    select(VAR_ID, P_VALUE, Risk_Allele, GWAS)

  print(nrow(vars))
  vars_sig = rbind(vars_sig, vars)
}
print(paste("No. total SNPs below pval cutoff:",nrow(vars_sig)))

# remove duplicates
vars_sig_uniq <- vars_sig %>%
  arrange(VAR_ID, P_VALUE) %>%
  filter(!duplicated(VAR_ID)) %>% # so we remove duplicates with the higher pvalue
  rename(PVALUE = P_VALUE)
print(paste("No. unique SNPs:",nrow(vars_sig_uniq)))

# remove indels
vars_sig_noIndels <- vars_sig_uniq %>%
  separate(VAR_ID, into=c("CHR","POS","REF","ALT"),sep="_",remove = F) %>%
  mutate(alleles = paste0(REF,ALT)) %>%
  subset(nchar(alleles)==2 | (nchar(alleles)<=4 & grepl(",",alleles))) %>%
  select(VAR_ID, PVALUE, Risk_Allele, GWAS)
print(paste("No. SNPs excluding indels:",nrow(vars_sig_noIndels)))

save.image(file = file.path(project_dir, "pipeline_data.RData"))

# #----

# SECTION 3: VARIANT PRUNING (LD-BASED)

# LD pruning
print("LD-pruning using EUR panel in LDlinkR::SNPclip...")
ld_prune(df_snps = vars_sig_noIndels,
                    pop = "EUR",
                    output_dir = project_dir,
                    r2 = 0.05,
                    maf=0.001,
                    my_token = user_token,
                    chr = c(1:22))

#----

# combine LD-pruning results
print("Combining SNP.clip results...")
ld_files <- list.files(path = project_dir,
                       pattern = "^snpClip_results",
                       full.names = T)

df_clipped_res = data.frame("RS_Number"=as.character(),
                         "Position"=as.character(),
                         "Alleles"= as.character(),
                         "Details"=as.character())

rename_cols_clipped <- c(RS_Number="RS Number",
                         Position_grch37="Position")

for (ld_file in ld_files){
  df <- fread(ld_file, stringsAsFactors = F, data.table = F) %>%
    dplyr::rename(any_of(rename_cols_clipped))
  df_clipped_res <- rbind(df_clipped_res, df)
}

df_clipped_kept <- df_clipped_res %>%
    filter(Details=="Variant kept.")

pruned_vars <- vars_sig_noIndels %>%
    separate(VAR_ID, into=c("CHR","POS","REF","ALT"), sep = "_",remove = F) %>%
    mutate(ChrPos = paste0("chr", CHR, ":", POS)) %>%
    filter(ChrPos %in% df_clipped_kept$Position)
print(sprintf("T2D SNPs pruned from %i to %i...", nrow(vars_sig_noIndels), nrow(pruned_vars)))

save.image(file = file.path(project_dir, "pipeline_data.RData"))

#----
# SECTION 4: VARIANT MISSINGNESS

print("Searching for variants in trait GWAS...")
gwas_variants <- pruned_vars$VAR_ID
df_Ns <- count_traits_per_variant(gwas_variants,
                                  trait_ss_files)

# fix column names
df_Ns_rev <- df_Ns %>%
  column_to_rownames("VAR_ID") %>%
  set_colnames(names(trait_ss_files))

print("Calculating variant missingess in traits...")
variant_counts_df <- data.frame(VAR_ID=rownames(df_Ns_rev),
                                frac=rowSums(!is.na(df_Ns_rev[,names(trait_ss_files)]))/length(trait_ss_files))
var_nonmissingness <- ifelse(
  gwas_variants %in% variant_counts_df$VAR_ID,
  # if in counts data frame, take the non-missing fraction:
  variant_counts_df$frac[match(gwas_variants, variant_counts_df$VAR_ID)],
  # else not in data frame, so non-missing fraction is 0:
  0
)
var_nonmissingness <- setNames(var_nonmissingness, gwas_variants)

save.image(file = file.path(project_dir, "pipeline_data.RData"))

#----
# SECTION 5: DETERMINE VARIANTS NEEDING PROXIES

print("Identifying variants needing proxies...")
proxies_needed_df <- find_variants_needing_proxies(pruned_vars,
                                                   var_nonmissingness,
                                                   rsID_map_file,
                                                   missing_cutoff = 0.8)

#----
# SECTION 6: PROXY SEARCH

print("Searching for proxies with TopLD API...")
proxy_search_results <- choose_proxies(need_proxies = proxies_needed_df,
                                       method="LDlink",
                                       LDlink_token = user_token,
                                    topLD_path = api_path,
                                    rsID_map_file = rsID_map_file,
                                    trait_ss_files = trait_ss_files,
                                    pruned_variants = pruned_vars,
                                    population="EUR"
)

df_proxies <- proxy_search_results %>%
  dplyr::select(VAR_ID, proxy_VAR_ID) %>%
  dplyr::inner_join(pruned_vars[,c("VAR_ID","GWAS")], by="VAR_ID") %>%
  mutate(Risk_Allele=NA, PVALUE=NA)

save.image(file = file.path(project_dir, "pipeline_data.RData"))

#----

# SECTION 7: Fetch summary statistics for SNPs in trait GWAS

print("Prepping input for fetch_summary_stats...")
df_orig_snps <- pruned_vars %>%
  filter(!VAR_ID %in% proxies_needed_df$VAR_ID)

df_input_snps <- rbind(df_orig_snps[,c("VAR_ID","PVALUE", "Risk_Allele", "GWAS")],
                       df_proxies[,c("VAR_ID","PVALUE", "Risk_Allele", "GWAS")]) %>%
  arrange(PVALUE) %>%
  filter(!duplicated(VAR_ID))

cat(sprintf("\n%i original SNPs...\n", nrow(df_orig_snps)))
cat(sprintf("\n%i proxy SNPs...\n", nrow(df_proxies)))
cat(sprintf("\n%i total unique SNPs!\n", nrow(df_input_snps)))

initial_zscore_matrices <- fetch_summary_stats(
  df_input_snps,
  main_ss_filepath,
  trait_ss_files,
  trait_ss_size,
  pval_cutoff=0.05
)


save.image(file = file.path(project_dir, "pipeline_data.RData"))
system(sprintf("mv alignment_GWAS_summStats.csv %s", project_dir))

#----

# Section 8: get rsIDs for final variant set

print("Getting rsIDs for final snps and saving to results...")
z_mat <- initial_zscore_matrices$z_mat
N_mat <- initial_zscore_matrices$N_mat

df_var_ids <- df_input_snps %>%
  separate(VAR_ID, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
  mutate(ChrPos=paste(Chr,Pos,sep = ":")) %>%
  subset(ChrPos %in% rownames(z_mat))
write(df_var_ids$VAR_ID,'my_snps.tmp')

system(sprintf("grep -wFf my_snps.tmp %s > %s",
               rsID_map_file, file.path(project_dir, "rsID_map.txt")))

df_rsIDs <- fread(cmd=sprintf("grep -wFf my_snps.tmp %s",rsID_map_file),
                  header = F,
                  col.names = c("VAR_ID","rsID"))
print(sprintf("rsIDs found for %i of %i SNPs...", nrow(df_rsIDs), nrow(df_var_ids)))

df_rsIDs_final <- df_rsIDs %>%
  filter(VAR_ID %in% df_var_ids$VAR_ID)
write_delim(x=df_rsIDs_final,
            file = file.path(project_dir, "rsID_map.txt"),
            col_names = T)

save.image(file = file.path(project_dir, "pipeline_data.RData"))

#----

# Section 9: Fill missing data in z-score and N matrices

df_snps <- df_input_snps %>%
  inner_join(df_rsIDs_final, by="VAR_ID") %>%
  data.frame()

print("Searching for cover proxies for missing z-scores...")
initial_zscore_matrices_final <- fill_missing_zscores(initial_zscore_matrices,
                                                      df_snps,
                                                      trait_ss_files,
                                                      trait_ss_size,
                                                      main_ss_filepath,
                                                      rsID_map_file,
                                                     method_fill="median",
                                                     population="EUR")
save.image(file = file.path(project_dir, "pipeline_data.RData"))

#----

# Section 10.) Generate non-negative z-score matrix

prep_z_output <- prep_z_matrix(z_mat = initial_zscore_matrices_final$z_mat,
                               N_mat = initial_zscore_matrices_final$N_mat,
                               corr_cutoff = 0.8)

# prep_z_output has two outputs:

#   1.) The scaled, non-negative z-score matrix
final_zscore_matrix <- prep_z_output$final_z_mat

#   2.) Results from the trait filtering
df_traits_filtered <- prep_z_output$df_traits
write_csv(x = df_traits_filtered,
          file = file.path(project_dir,"df_traits.csv"))

# prep_z_matrix also save trait correlation matrix to working dir, so move to project dir
system(sprintf("mv trait_cor_mat.txt %s", project_dir))

print(sprintf("Final matrix: %i SNPs x %i traits",
      nrow(final_zscore_matrix),
      ncol(final_zscore_matrix)/2))

save.image(file = file.path(project_dir, "pipeline_data.RData"))

#----

# Section 11.) Run bNMF 
bnmf_reps <- run_bNMF(final_zscore_matrix,
                      n_reps=25,
                      tolerance = 1e-6)
summarize_bNMF(bnmf_reps, dir_save=project_dir)

save.image(file = file.path(project_dir, "pipeline_data.RData"))

end=Sys.time()
print("Total pipeline runtime:")
print(end-start)

#----

# format results
k <- NULL
if (is.null(k)){
  html_filename <- "results_for_maxK.html"
} else {
  html_filename <- sprintf("results_for_K_%i.html", k)
}

rmarkdown::render(
  './format_bNMF_results.Rmd',
  output_file = html_filename,
  params = list(main_dir = project_dir,
                k = k,
                loci_file="query",
                GTEx=F,
                my_traits=gwas_traits)
)


#----
