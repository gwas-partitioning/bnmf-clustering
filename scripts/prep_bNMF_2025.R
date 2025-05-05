library(tidyverse)
library(data.table)
library(dplyr)


fetch_summary_stats <- function(df_input, gwas_ss_file, trait_ss_files, trait_ss_size=NULL, pval_cutoff=1, read_trait_method = "datatable") {
  
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
  

  read_single_trait_dt <- function(trait, read_trait_input) {
    message(sprintf("Processing %s...", trait))
    
    # Read the header line (as a character vector)
    header_line <- readLines(trait_ss_files[[trait]], n = 1)
    headers <- strsplit(header_line, "\t")[[1]]
    
    # Read the filtered data from the trait file using an external grep command.
    # (Assumes that "all_snps_pos.tmp" contains your list of VAR_ID values in the same format as in the file.)
    cmd_str <- if (endsWith(trait_ss_files[[trait]], ".gz")) {
      sprintf("gzip -cd %s | grep -Ff %s", trait_ss_files[[trait]], tmp_var_file)
    } else {
      sprintf("grep -Ff %s %s", tmp_var_file, trait_ss_files[[trait]])
    }
    
    # Use fread to read the results from the command into a data.table.
    dt <- fread(cmd = cmd_str, header = FALSE, col.names = headers, stringsAsFactors = FALSE)

        # If the trait_ss_size for this trait is available, assign it to column N_PH.
    if (!is.na(trait_ss_size[[trait]])) {
      dt[, N_PH := trait_ss_size[[trait]]]
    }
    
    # If no rows were found, create a dummy data.table.
    if (nrow(dt) == 0) {
      dt <- data.table(VAR_ID = NA, BETA = NA, SE = NA, P_VALUE = NA, N = NA)
    }
    
    # Rename columns as needed. (Note: data.table's setnames is used here.)
    # Our rename mapping: VAR_ID -> VAR_ID_hg19, N_PH -> N, BETA -> beta, SE -> se,
    # P_VALUE -> FTOP_p_value_meta, z -> Zscore. (The duplicate mapping for "z" is ignored.)

    rename_map <- list("VAR_ID_hg19"="VAR_ID", "N" = "N_PH", "beta" = "BETA",
                       "se" = "SE", "FTOP_p_value_meta" = "P_VALUE", "Zscore" = "z")
    for (old in names(rename_map)) {
      if (old %in% names(dt)) {
        setnames(dt, old, rename_map[[old]])
      }
    }
    
    # Remove duplicates based on VAR_ID_hg19.
    dt <- unique(dt, by = "VAR_ID")
    
    # Separate VAR_ID into CHR, POS, REF, ALT using tstrsplit.
    dt[, c("CHR", "POS", "REF", "ALT") := tstrsplit(VAR_ID, "_", fixed = TRUE)]
    
    # Create the SNP column as "CHR:POS".
    dt[, SNP := paste0(CHR, ":", POS)]
    
    # Create Effect_Allele_PH: if the column exists, leave it; otherwise, set to ALT.
    if (!("Effect_Allele_PH" %in% names(dt))) {
      dt[, Effect_Allele_PH := ALT]
    }
    
    # Keep only the necessary columns.
    keep_cols <- c("SNP", "Effect_Allele_PH", "REF", "ALT", "N_PH", "BETA", "SE", "P_VALUE", "z")
    keep_cols <- keep_cols[keep_cols %in% names(dt)]
    dt <- dt[, ..keep_cols]
    
    # Perform a right join with read_trait_input on "SNP".
    # Ensure that both dt and read_trait_input are data.tables keyed by "SNP".
    setkey(dt, SNP)
    setkey(read_trait_input, SNP)
    dt <- read_trait_input[dt]  # This is equivalent to a right join: all rows from dt, matched to read_trait_input.
    
    # Filter rows where concatenated alleles match:
    # Keep rows where paste0(REF, ALT) equals either paste0(Risk_Allele, Nonrisk_Allele)
    # or paste0(Nonrisk_Allele, Risk_Allele). (Risk_Allele and Nonrisk_Allele come from read_trait_input.)
    dt[, allele_combo := paste0(REF, ALT)]
    dt <- dt[(allele_combo == paste0(Risk_Allele, Nonrisk_Allele)) | 
               (allele_combo == paste0(Nonrisk_Allele, Risk_Allele))]
    
    # Adjust P_VALUE: if P_VALUE is zero, replace with minimum nonzero value.
    if (any(dt$P_VALUE == 0, na.rm = TRUE)) {
      min_val <- min(dt$P_VALUE[dt$P_VALUE > 0], na.rm = TRUE)
      dt[P_VALUE == 0, P_VALUE := min_val]
    }
    
    # If the Zscore column is missing or all NA, compute it as beta/se (after replacing zeros in se).
    if (!("Zscore" %in% names(dt)) || all(is.na(dt$Zscore))) {
      if (any(dt$SE == 0, na.rm = TRUE)) {
        min_se <- min(dt$SE[dt$SE > 0], na.rm = TRUE)
        dt[SE == 0, SE := min_se]
      }
      dt[, z := BETA / SE]
    }
    
    # Align Zscore sign based on effect allele:
    # If Effect_Allele_PH equals Risk_Allele then keep Zscore;
    # if equals Nonrisk_Allele then flip sign; otherwise set to NA.
    dt[, z := ifelse(Effect_Allele_PH == Risk_Allele, z,
                          ifelse(Effect_Allele_PH == Nonrisk_Allele, -z, NA_real_))]
    
    # Keep only the desired columns and rename as needed.
    dt <- dt[, .(SNP, z , N_PH, P_VALUE)]
    
    # Ensure P_VALUE is numeric.
    dt[, P_VALUE := as.numeric(P_VALUE)]
    
    # Remove rows with any missing values.
    dt <- na.omit(dt)
    
    return(dt)
  }
  
  
  
  print("Reading original GWAS summary statistics...")
  variant_vec <- df_input$VAR_ID
  tmp_var_file <- "variants_to_query.tmp"
  
  print("Writing initials SNPs to file for grepping...")
  writeLines(variant_vec, tmp_var_file)
  
  
  if (class(gwas_ss_file)=="character"){
    
    print("Reading headers...")
    headers <- as.character(fread(gwas_ss_file, nrows=1,
                                  data.table=F, stringsAsFactors=F, header=F))
    print(headers)
    
    print("Grepping file...")
    if (endsWith(gwas_ss_file,".gz")) {
      gwas_ss <- fread(cmd = sprintf("gzip -cd %s | grep -Ff variants_to_query.tmp", gwas_ss_file),
                       header=F, col.names=headers,
                       data.table = FALSE, stringsAsFactors = FALSE)
      
    } else {
      gwas_ss <- fread(cmd = sprintf("fgrep -Fwf variants_to_query.tmp %s", gwas_ss_file),
                       header=F, col.names=headers,
                       data.table = FALSE, stringsAsFactors = FALSE)
    }

    file.remove(tmp_var_file)
  } else {
    # input dataframe with VAR_ID, P_VALUE, SNP, Risk_Allele, Nonrisk_Allel
    gwas_ss <- gwas_ss_file %>%
      filter(VAR_ID %in% variant_vec)
  }

  print(sprintf("%i of %i SNPs found in main GWAS...", nrow(gwas_ss), nrow(df_input)))
  print(head(gwas_ss))
  
  if (!all(c("SNP","Risk_Allele","Nonrisk_Allele") %in% names(gwas_ss))){
    print("Separating VAR_ID...")
    gwas_ss <- gwas_ss %>%
      separate(VAR_ID, into=c("Chr","Pos","REF","ALT"), sep = "_", remove = F)
    
    if (!"BETA" %in% colnames(gwas_ss)){
      # print("Converting Odds Ratio to Log Odds Ratio...")
      gwas_ss <- gwas_ss %>%
        mutate(BETA = log(ODDS_RATIO))
    }
    # print("Creating SNP and Risk Allele columns...")
    gwas_ss <- gwas_ss %>%
      mutate(SNP = paste(Chr,Pos,sep=":")) %>%
      mutate(Risk_Allele=ifelse(BETA > 0, ALT, REF),
             Nonrisk_Allele=ifelse(BETA > 0, REF, ALT))
  } 
  
  my_vars <- c('SNP', 'ALT', 'REF', 'Risk_Allele', 'Nonrisk_Allele', 'P_VALUE', 'BETA', 'SE', 'z')
  gwas_ss <- gwas_ss %>%
    select(any_of(my_vars))
  
  print("Merging formatted GWAS with final variant vector...")
  
  print("Rename risk_allelel to Risk_Allele_Orig...")
  
  df_input_wGWAS <- df_input %>%
    dplyr::rename(Risk_Allele_Orig=Risk_Allele) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_") %>%
    mutate(SNP=paste(CHR, POS, sep=":")) %>%
    select(SNP, Risk_Allele_Orig) %>%
    inner_join(gwas_ss, by="SNP") %>%
    mutate_at(vars(P_VALUE),as.numeric)
  
  print(paste0(nrow(df_input_wGWAS), " of ", length(variant_vec)," variants are available in the primary GWAS."))
  print(sprintf("Max p-value in primary GWAS: %.3e", max(df_input_wGWAS$P_VALUE)))
  
  pval_bonf <- 0.05/nrow(df_input_wGWAS)
  
  opp_risk <- df_input_wGWAS %>%
    filter(between(P_VALUE, pval_bonf, pval_cutoff) & Risk_Allele != Risk_Allele_Orig)
  print(sprintf("%i variants with opposite risk allele and above Bonferroni cutoff...", nrow(opp_risk)))
  
  high_pval <- df_input_wGWAS %>%
    filter(P_VALUE > pval_cutoff)
  print(sprintf("%i variants above absolute p-value cutoff...", nrow(high_pval)))
  
  uniq_to_drop <- unique(c(opp_risk$SNP, high_pval$SNP))
  print(sprintf("%i unique variants being dropped due to risk allele or p-value...", length(uniq_to_drop)))
  
  df_input_wGWAS_filtered <- df_input_wGWAS %>%
    filter(P_VALUE < pval_bonf | between(P_VALUE, pval_bonf, pval_cutoff) & Risk_Allele == Risk_Allele_Orig)
  print(sprintf("%i variants after p-value and risk allele filter...", nrow(df_input_wGWAS_filtered)))
  
  df_read_trait_input <- df_input_wGWAS_filtered %>% 
    dplyr::select(SNP, Risk_Allele, Nonrisk_Allele) %>%
    as.data.table()
  print(paste(nrow(df_read_trait_input), "remaining SNPs after p-value filtering"))
  
  print("Writing final SNPs to file for grepping...")
  writeLines(gsub(":","_",df_read_trait_input$SNP), tmp_var_file)

  print("Retrieving z-scores and sample sizes for each trait...")

  # run read_single_trait
  trait_df_list <- lapply(names(trait_ss_files), read_single_trait_dt, df_read_trait_input) %>%
    setNames(names(trait_ss_files)) %>%
    bind_rows(.id="trait")  # Bind all processed trait datasets into a single "long" data frame

  z_df_wide <- trait_df_list %>%
    select(trait, SNP, z) %>%
    pivot_wider(names_from="trait", values_from="z")
  z_mat <- as.matrix(z_df_wide[, -1])
  rownames(z_mat) <- z_df_wide$SNP

  N_df_wide <- trait_df_list %>%
    select(trait, SNP, N_PH) %>%
    pivot_wider(names_from="trait", values_from="N_PH")
  N_mat <- as.matrix(N_df_wide[, -1])
  rownames(N_mat) <- N_df_wide$SNP

  system(paste("rm",tmp_var_file))

  df_z <- data.frame(z_mat) %>%
    set_colnames(names(trait_ss_files))
  df_N <- data.frame(N_mat) %>%
    set_colnames(names(trait_ss_files))

  print("Do all REF/ALT alleles match risk alleles?")
  print(all(with(df_input_wGWAS_filtered, paste0(Risk_Allele, Nonrisk_Allele)==paste0(REF, ALT) |
                   paste0(Risk_Allele, Nonrisk_Allele)==c(paste0(ALT, REF)))))
  print("Done!")
  return(list(df_z=df_z, df_N=df_N, df_gwas=df_input_wGWAS_filtered))
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

  # 3.) Remove traits with low variance and low SNR
  # Compute variance and median absolute z-score for each trait.
  trait_variances <- apply(z_mat, 2, function(x) var(x, na.rm = TRUE))
  trait_median_abs <- apply(z_mat, 2, function(x) median(abs(x), na.rm = TRUE))
  
  # Compute a simple SNR: median absolute / variance (you could use MAD instead)
  trait_snr <- trait_median_abs / trait_variances
  
  # Set thresholds—these thresholds need to be tuned based on your data.
  variance_threshold <- 1
  snr_threshold <- 1  # for example
  
  # Identify traits to keep: keep if variance is above the threshold or if SNR is high.
  traits_to_keep <- names(z_mat)[(trait_variances >= variance_threshold) | (trait_snr >= snr_threshold) ]
  traits_to_remove <- names(z_mat)[!names(z_mat) %in% traits_to_keep]
  z_mat <- z_mat[, traits_to_keep]
  
  cat("Traits removed for low variance/SNR:", paste(traits_to_remove, collapse = "\n"), "\n")
  df_traits_filtered <- rbind(df_traits_filtered,
                              data.frame(trait = traits_to_remove,
                                         result = "removed (low signal)",
                                         note = NA)
  )
  
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
  
  output <- list(final_z_mat=final_z_mat,
                 df_traits=df_traits_filtered)
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
