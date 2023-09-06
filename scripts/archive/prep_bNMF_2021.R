library(tidyverse)
library(data.table)

fetch_summary_stats <- function(variant_vec, gwas_ss_file, trait_ss_files) {
  
  # Given a final (pruned & proxied) set of variants to be clustered, 
  # fetch z-scores and sample size info from summary statistics for each of a 
  # series of traits
  # INPUTS:
  #   - variant_vec: vector of variants to be clustered
  #   - gwas_ss: filepath of  with VAR_IDs and betas from the original GWAS
  #   - trait_ss_vec: named vector of trait summary statistic filepaths
  # Final variant vector should be in VAR_ID format: [CHR]_[POS]_[REF]_[ALT] (using hg19)
  # GWAS summary statistic data frame must have at least the following 
  # columns: SNP (CHR:POS), REF, ALT, BETA
  # Each trait summary statistic dataset must have the following 
  # columns: VAR_ID, Effect_Allele_PH, BETA, SE, P_VALUE, N_PH
  
  # ISSUES:
  # - potential for strand-flip?
  
  read_single_trait <- function(trait, variant_df) {
    # Read/filter/process summary statistics for a single trait
    print(paste0("Processing ", trait, "..."))
    df <- fread(trait_ss_files[[trait]], data.table=F, stringsAsFactors=F)
    # if (grepl("\\/UKBB\\/", trait_ss_path)) {
    #   df <- df %>%
    #     mutate(VAR_ID=gsub(":", "_", variant)) %>%
    #     filter(VAR_ID %in% final_variant_ss$VAR_ID) %>%
    #     separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove=F) %>%
    #     select(VAR_ID, Effect_Allele_PH=ALT, BETA=beta, SE=se, P_VALUE=pval, N_PH=n_complete_samples)
    # }
    # if (!("N_PH" %in% names(df))) df$N_PH <- as.integer(NA)
    df %>%
      separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_") %>%
      mutate(SNP=paste(CHR, POS, sep=":")) %>%
      select(SNP, Effect_Allele_PH, N_PH, BETA, SE, P_VALUE) %>%
      right_join(variant_df, by="SNP", suffix=c(".gwas", ".trait")) %>%
      mutate(z=BETA / SE,  # First, calculate z-score magnitude
             z=case_when(  # Next, align z-score sign with GWAS phenotype-raising allele
               Effect_Allele_PH == Risk_Allele ~ z,
               Effect_Allele_PH == Nonrisk_Allele ~ -z,
               TRUE ~ as.numeric(NA)  # For example, if trait effect allele matches neither REF nor ALT from GWAS
             )) %>%
      select(SNP, z, N_PH, P_VALUE)
  }
  
  print("Retrieving risk alleles from the original GWAS summary statistics...")
  gwas_ss <- fread(gwas_ss_file, data.table=F, stringsAsFactors=F) %>%
    mutate(Risk_Allele=ifelse(BETA > 0, ALT, REF),
           Nonrisk_Allele=ifelse(BETA > 0, REF, ALT)) %>%
    select(SNP, Risk_Allele, Nonrisk_Allele)
  variant_df <- tibble(VAR_ID=variant_vec) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_") %>%
    mutate(SNP=paste(CHR, POS, sep=":")) %>%
    select(SNP) %>%
    inner_join(gwas_ss, by="SNP")
  print(paste0(nrow(variant_df), " of ", length(variant_vec),
               " variants are available in the primary GWAS."))
    
  
  print("Retrieving z-scores and sample sizes for each trait...")
  trait_df_long <- lapply(names(trait_ss_files), read_single_trait, variant_df) %>%
    setNames(names(trait_ss_files)) %>%
    bind_rows(.id="trait")  # Bind all processed trait datasets into a single "long" data frame
  
  z_df_wide <- trait_df_long %>%
    select(trait, SNP, z) %>%
    pivot_wider(names_from="trait", values_from="z")
  z_mat <- as.matrix(z_df_wide[, -1])
  rownames(z_mat) <- z_df_wide$SNP
  
  N_df_wide <- trait_df_long %>%
    select(trait, SNP, N_PH) %>%
    pivot_wider(names_from="trait", values_from="N_PH")
  N_mat <- as.matrix(N_df_wide[, -1])
  rownames(N_mat) <- N_df_wide$SNP
  
  # P_df <- trait_df_long %>%
  #   group_by(trait) %>%
  #   summarise(minP=min(P_VALUE, na.rm=T))
  
  list(z_mat=z_mat, N_mat=N_mat)
       # minP_vec=setNames(P_df$minP, P_df$trait))
}


prep_z_matrix <- function(z_mat, N_mat) {
  
  # Given a matrix of z-scores (N_variants x M_traits) and vector of median
  # sample sizes per trait:
  # 1) perform final pre-processing steps before bNMF clustering:
  # trait filtering by p-value, trait pruning based on correlation,
  # and z-score scaling based on sample size
  # 2) expand N x M matrix into N x 2M non-negative matrix
  
  # Filter traits by p-value (min. p-value < 0.05/N_variants)
  minP_vec <- apply(z_mat, 2, function(x) min(2 * pnorm(abs(x), lower.tail=F), na.rm=T))
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
  
  # Adjust z-scores by sample size for each variant-trait combo
  # i.e. (z = z / sqrt(medN) * mean(sqrt(medN_all_traits)))
  print("Performing sample size adjustment...")
  medN_vec <- apply(N_mat[, colnames(z_mat)], 2, median, na.rm=T)
  z_mat <- z_mat / sqrt(N_mat[, colnames(z_mat)]) * mean(sqrt(medN_vec))

  
  # Replace missing values with zero
  print("Replacing remaining missing values with zero...")
  print(paste0(sum(is.na(z_mat)), " missing values were replaced."))
  z_mat[is.na(z_mat)] <- 0
  
  # Expand into N x 2M non-negative matrix
  print("Expanding z-score matrix into non-negative matrix (N-variants x 2M-traits)...")
  z_mat_pos <- z_mat
  z_mat_pos[z_mat_pos < 0] <- 0
  colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
  z_mat_neg <- -z_mat
  z_mat_neg[z_mat_neg < 0] <- 0
  colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
  final_z_mat <- cbind(z_mat_pos, z_mat_neg)
  
  # Write N x M and N x 2M matrices
  saveRDS(z_mat, "z_score_mat.rds")
  saveRDS(final_z_mat, "z_score_mat_nonnegative.rds")
  
  final_z_mat
}
