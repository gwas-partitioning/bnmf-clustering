library(tidyverse)


##### PREPROCESSING #####

prep_z_matrix <- function(z_mat, minP_vec, medN_vec) {
  
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
  print(paste0("Multiplying sample size-adjusted z-score matrix by ", 
               round(mean(sqrt(medN_vec[pruned_traits]))), " (i.e. mean(sqrt(median(N))))"))
  z_mat <- z_mat * mean(sqrt(medN_vec[pruned_traits]))
  
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