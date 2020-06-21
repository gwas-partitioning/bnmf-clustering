library(tidyverse)
library(data.table)

fetch_summary_stats <- function(trait_ss_list, final_variant_ss) {
  
  # Given a dataset of final (pruned & proxied) variants to cluster on, 
  # fetch z-scores and sample size info from summary statistics for each of a 
  # series of traits
  # Final variant dataset should have: VAR_ID, REF, ALT (effect), BETA
  # FOR NOW: each trait summary statistic dataset must have the following 
  # columns: VAR_ID, Effect_Allele_PH, BETA, SE, P_VALUE, N_PH
  # VAR_IDs should be consistent with the format of IDs for variants to be 
  # clustered, and here will be: [CHR]_[POS]_[REF]_[ALT] (using hg19)
  
  read_single_trait <- function(trait_ss_path, final_variant_ss) {
    # Read/filter/process summary statistics for a single trait
    print(paste0("Processing ", trait_ss_path, "..."))
    df <- fread(trait_ss_path, data.table=F, stringsAsFactors=F)
    if (grepl("\\/UKBB\\/", trait_ss_path)) {
      df <- df %>%
        mutate(VAR_ID=gsub(":", "_", variant)) %>%
        filter(VAR_ID %in% final_variant_ss$VAR_ID) %>%
        separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove=F) %>%
        select(VAR_ID, Effect_Allele_PH=ALT, BETA=beta, SE=se, P_VALUE=pval, N_PH=n_complete_samples)
    }
    if (!("N_PH" %in% names(df))) df$N_PH <- as.integer(NA)
    df %>%
      right_join(final_variant_ss, by="VAR_ID", suffix=c(".gwas", ".trait")) %>%
      mutate(z=BETA.trait / SE,  # First, calculate z-score magnitude
             z=case_when(  # Next, align z-score sign with GWAS phenotype-raising allele
               (Effect_Allele_PH == ALT) & (BETA.gwas > 0) ~ z,
               (Effect_Allele_PH == REF) & (BETA.gwas > 0) ~ -z,
               (Effect_Allele_PH == ALT) & (BETA.gwas < 0) ~ -z,
               (Effect_Allele_PH == REF) & (BETA.gwas < 0) ~ z,
               TRUE ~ as.numeric(NA)  # For example, if trait effect allele matches neither REF nor ALT from GWAS
             )) %>%
      select(VAR_ID, z, N_PH, P_VALUE)
  }
  
  trait_df_long <- lapply(trait_ss_list, read_single_trait, final_variant_ss) %>%
    bind_rows(.id="trait")  # Bind all processed trait datasets into a single "long" data frame
  
  z_df_wide <- trait_df_long %>%
    select(trait, VAR_ID, z) %>%
    pivot_wider(names_from="trait", values_from="z")
  z_mat <- as.matrix(z_df_wide[, -1])
  rownames(z_mat) <- z_df_wide$VAR_ID
  
  N_df_wide <- trait_df_long %>%
    select(trait, VAR_ID, N_PH) %>%
    pivot_wider(names_from="trait", values_from="N_PH")
  N_mat <- as.matrix(N_df_wide[, -1])
  rownames(N_mat) <- N_df_wide$VAR_ID
  
  P_df <- trait_df_long %>%
    group_by(trait) %>%
    summarise(minP=min(P_VALUE, na.rm=T))
  
  list(z_mat=z_mat, N_mat=N_mat, minP_vec=setNames(P_df$minP, P_df$trait))
}


prep_z_matrix <- function(z_mat, N_mat, minP_vec) {
  
  # Given a matrix of z-scores (N_variants x M_traits) and vector of median
  # sample sizes per trait:
  # 1) perform final pre-processing steps before bNMF clustering:
  # trait filtering by p-value, trait pruning based on correlation,
  # and z-score scaling based on sample size
  # 2) expand N x M matrix into 2N x M non-negative matrix
  
  # Filter traits by p-value (min. p-value < 0.05/N_variants)
  stopifnot(all(colnames(z_mat) == names(minP_vec)))
  print(paste0("Removing traits with no variant having p < 0.05 / # variants: ",
               paste(colnames(z_mat)[minP_vec >= 0.05 / nrow(z_mat)], 
                     collapse=", ")))
  z_mat <- z_mat[, minP_vec < 0.05 / nrow(z_mat)]
  
  # Prune traits by correlation (remove traits with Pearson |r| > 0.85)
  trait_cor_mat <- cor(z_mat, use="pairwise.complete.obs")  # Trait-trait correlation matrix
  trait_min_pvals <- minP_vec[names(minP_vec) %in% colnames(z_mat)]  # Remove filtered traits
  remaining_traits <- names(sort(trait_min_pvals))
  keep_traits <- c()
  while (length(remaining_traits) > 0) {
    # print(remaining_traits[1])
    keep_traits <- c(keep_traits, remaining_traits[1])
    remaining_traits <- setdiff(
      remaining_traits, 
      rownames(trait_cor_mat)[abs(trait_cor_mat[, remaining_traits[1]]) >= 0.85]
    )
  }
  pruned_traits <- setdiff(colnames(z_mat), keep_traits)
  print(paste("Traits removed in pruning process:", 
              paste(pruned_traits, collapse=", ")))
  z_mat <- z_mat[, keep_traits]
  
  # Adjust z-scores by sample size for each variant-trait combo (z = z / sqrt(N))
  # z_mat <- t(t(z_mat) / sqrt(N_vec[match(pruned_traits, colnames(z_mat))]))
  # Multiply full matrix by mean(sqrt(median(N))) (a single number for the whole matrix)
  medN_vec <- apply(N_mat, 2, median, na.rm=T)
  z_mat <- z_mat / sqrt(N_mat) * mean(sqrt(medN_vec))
  # print(paste0("Multiplying sample size-adjusted z-score matrix by ", 
  #              round(mean(sqrt(medN_vec[pruned_traits]))), " (i.e. mean(sqrt(median(N))))"))
  
  # Replace missing values with zero
  z_mat[is.na(z_mat)] <- 0
  
  # Expand into 2N x M non-negative matrix
  z_mat_pos <- z_mat
  z_mat_pos[z_mat_pos < 0] <- 0
  colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
  z_mat_neg <- -z_mat
  z_mat_neg[z_mat_neg < 0] <- 0
  colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
  final_z_mat <- cbind(z_mat_pos, z_mat_neg)
  
  final_z_mat
}
