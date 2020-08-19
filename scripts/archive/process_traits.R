library(readxl)
library(tidyverse)


fetch_variant_fracs <- function(trait_filepath_list, gwas_varIDs) {
  
  # Given a named list of filepaths to trait GWAS summary stats and a list of
  # primary GWAS variant IDs, return a named vector indicating the fraction of 
  # those traits having each variant
  
  trait_variants <- lapply(trait_filepaths, function(f) {
    read_tsv(f, col_types=cols_only(VAR_ID="c"), n_max=100000) %>%
      filter(VAR_ID %in% gwas_varIDs)
  })
  var_presence_df <- do.call(bind_rows, c(trait_variants, .id="trait")) %>%
    mutate(var_present=1) %>%
    pivot_wider(names_from="trait", values_from="var_present", values_fill=0)  # 0/1 values: is variant (row) present in trait GWAS (column)?
  var_presence_mat <- as.matrix(var_presence_df[, 2:ncol(var_presence_df)])
  rownames(var_presence_mat) <- var_presence_df$VAR_ID
  var_fracs <- rowSums(var_presence_mat) / length(trait_variants)  # Fraction of traits having each variant
  var_fracs
}


traits_doc <- read_excel("../data/clustering_data_source.xlsx", sheet=3)
trait_filepaths <- setNames(traits_doc$full_path, traits_doc$trait_name)

gwas_varIDs <- sample(chosen_variants$varID, size=100)

variant_fracs <- fetch_variant_fracs(trait_filepaths, gwas_varIDs)




