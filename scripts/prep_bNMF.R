library(tidyverse)
library(data.table)

fetch_summary_stats <- function(df_variants, gwas_ss_file, trait_ss_files, trait_ss_size=NULL, pval_cutoff=1) {
  
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
    headers <- as.character(fread(trait_ss_files[[trait]], nrows=1,
                                  data.table=F, stringsAsFactors=F, header=F))
    if (endsWith(trait_ss_files[[trait]],".gz")) {
      df <- fread(cmd=sprintf("gzip -cd %s | grep -Ff all_snps_pos.tmp ", trait_ss_files[[trait]]),
                          header=F, col.names=headers,
                          data.table=F, stringsAsFactors=F)
    } else {
      df <- fread(cmd=sprintf("grep -Ff all_snps_pos.tmp %s ", trait_ss_files[[trait]]),
                      header=F, col.names=headers,
                      data.table=F, stringsAsFactors=F)
    }
    if (!is.na(trait_ss_size[[trait]])) {     #(!("N_PH" %in% names(df)))
      cat("Using trait_ss_size for GWAS size!\n")
      df$N_PH <- trait_ss_size[[trait]]
    }
    rename_cols <- c(VAR_ID="VAR_ID_hg19",
                     N_PH="N",
                     BETA="beta",
                     SE="se",
                     P_VALUE="FTOP_p_value_meta",
                     z="Zscore",
                     z="ZSCORE")
    
    print(head(df))
    df <- df %>%
      rename(any_of(rename_cols)) %>%
      filter(!duplicated(VAR_ID)) %>%
      separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_") %>%
      mutate(SNP=paste(CHR, POS, sep=":")) %>%
      mutate(Effect_Allele_PH = if("Effect_Allele_PH" %in% names(.)) Effect_Allele_PH else ALT) %>%
      select(any_of(c('SNP', 'Effect_Allele_PH', 'REF', 'ALT', 'N_PH', 'BETA', 'SE', 'P_VALUE', 'z'))) %>%
      right_join(variant_df, by="SNP", suffix=c(".gwas", ".trait")) %>%
      subset(paste0(REF,ALT) %in% c(paste0(Risk_Allele, Nonrisk_Allele),
                                    paste0(Nonrisk_Allele,Risk_Allele))) %>% #if trait effect allele matches neither REF nor ALT from GWAS
      mutate(P_VALUE = ifelse(P_VALUE==0, min(P_VALUE[P_VALUE>0]), P_VALUE))
    if (!"z" %in% names(df)){
      df <- df %>%
        mutate(SE = ifelse(SE==0, min(SE[SE>0]), SE)) %>%
        mutate(z = BETA / SE)  # First, calculate z-score magnitude
    }
    df <- df %>%
      mutate(z = case_when(  # Next, align z-score sign with GWAS phenotype-raising allele
               Effect_Allele_PH == Risk_Allele ~ z,
               Effect_Allele_PH == Nonrisk_Allele ~ -z,
               TRUE ~ as.numeric(NA)  # For example, if trait effect allele matches neither REF nor ALT from GWAS
             )) %>%
      select(SNP, z, N_PH, P_VALUE) %>%
      drop_na()
  }

  variant_vec <- df_variants$VAR_ID
  
  if (class(gwas_ss_file)=="character"){
    print("Reading original GWAS summary statistics (file path)...")
    gwas_ss <- fread(gwas_ss_file, data.table=F, stringsAsFactors=F) %>%
      filter(VAR_ID %in% variant_vec)
  } else {
    # input dataframe with VAR_ID, P_VALUE, SNP, Risk_Allele, Nonrisk_Allel
    print("Reading original GWAS summary statistics (data.frame)...")
    gwas_ss <- gwas_ss_file %>%
      filter(VAR_ID %in% variant_vec)
  }

  print(sprintf("%i of %i SNPs found in main GWAS...", nrow(gwas_ss), nrow(df_variants)))
  
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
  
  print(head(gwas_ss))
  
  print("Merging formatted GWAS with final variant vector...")
  variant_df <- df_variants %>%
    rename(Risk_Allele_Orig=Risk_Allele) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_") %>%
    mutate(SNP=paste(CHR, POS, sep=":")) %>%
    select(SNP, Risk_Allele_Orig) %>%
    inner_join(gwas_ss, by="SNP")
  print(paste0(nrow(variant_df), " of ", length(variant_vec)," variants are available in the primary GWAS."))
  print(sprintf("Max p-value in primary GWAS: %.3e", max(variant_df$P_VALUE)))
  
  pval_bonf <- 0.05/nrow(variant_df)
  
  opp_risk <- variant_df %>%
    filter(between(P_VALUE, pval_bonf, pval_cutoff) & Risk_Allele != Risk_Allele_Orig)
  print(sprintf("%i variants with opposite risk allele and above Bonferroni cutoff...", nrow(opp_risk)))
  
  high_pval <- variant_df %>%
    filter(P_VALUE > pval_cutoff)
  print(sprintf("%i variants above absolute p-value cutoff...", nrow(high_pval)))
  
  uniq_to_drop <- unique(c(opp_risk$SNP, high_pval$SNP))
  print(sprintf("%i unique variants being dropped due to risk allele or p-value...", length(uniq_to_drop)))
  
  variant_df <- variant_df %>%
    filter(P_VALUE < pval_bonf | between(P_VALUE, pval_bonf, pval_cutoff) & Risk_Allele == Risk_Allele_Orig)
  print(sprintf("%i variants after p-value and risk allele filter...", nrow(variant_df)))
  
  print("Saving aligned & filtered GWAS to file...")
  write_csv(x = variant_df,
            file = "./alignment_GWAS_summStats.csv",
            col_names = T)
  
  variant_df <- variant_df %>% 
    select(SNP, Risk_Allele, Nonrisk_Allele)
  print(paste(nrow(variant_df), "remaining SNPs after p-value filtering"))
  
  print("Writing final SNPs to file for grepping...")
  write(gsub(":","_",variant_df$SNP), "all_snps_pos.tmp")

  print("Retrieving z-scores and sample sizes for each trait...")
  trait_df_long <- lapply(names(trait_ss_files), read_single_trait, variant_df) %>%
    setNames(names(trait_ss_files)) %>%
    bind_rows(.id="trait")  # Bind all processed trait datasets into a single "long" data frame
  # saveRDS(trait_df_long, file = "my_trait_df_long.rds")
  
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
  
  system("rm all_snps_pos.tmp")
  # P_df <- trait_df_long %>%
  #   group_by(trait) %>%
  #   summarise(minP=min(P_VALUE, na.rm=T))
  
  list(z_mat=z_mat, N_mat=N_mat)
       # minP_vec=setNames(P_df$minP, P_df$trait))
}


prep_z_matrix <- function(z_mat, N_mat,
                          corr_cutoff=0.85,
                          rm_traits=NULL,
                          pval_cutoff=NULL) {
  
  # Given a matrix of z-scores (N_variants x M_traits) and vector of median
  # sample sizes per trait:
  # 1) perform final pre-processing steps before bNMF clustering:
  # trait filtering by p-value, trait pruning based on correlation,
  # and z-score scaling based on sample size
  # 2) expand N x M matrix into N x 2M non-negative matrix
  
  
  # if no user input for pval_cutoff, use bonferroni definition
  if (is.null(pval_cutoff)){
    pval_cutoff <- 0.05 / nrow(z_mat) 
  }
  print(paste0(sum(is.na(z_mat)), " missing values before pval pruning."))
  
  df_traits <- data.frame(trait = colnames(z_mat))
  
  df_traits_filtered <- data.frame(trait=as.character(),
                           result=as.character(),
                           note=as.character())
  
  # remove manually entered traits
  if (!is.null(rm_traits)) {
    print(sprintf("Removing %i manually entered traits!!!", length(rm_traits)))
    z_mat <- z_mat[, !colnames(z_mat) %in% rm_traits]
    
    df_remove_manual <- data.frame(trait=rm_traits,
                                   result="removed (manual)",
                                   note=NA)
    df_traits_filtered <- rbind(df_traits_filtered, df_remove_manual)
  }
  
  # Filter traits by p-value (min. p-value < 0.05/N_variants)
  print("Filtering traits w/ no pvalues below cutoff...")
  minP_vec <- apply(z_mat, 2, function(x) min(2 * pnorm(abs(x), lower.tail=F), na.rm=T))
  traits_removed <- colnames(z_mat)[minP_vec >= pval_cutoff]
  df_lowP_vec <- data.frame(minPval=minP_vec[minP_vec >= pval_cutoff]) %>%
    rownames_to_column('trait') %>%
    mutate(result="removed (p-value)") %>% 
    mutate(note=paste("min pval=",format(minPval,scientific=T))) %>%
    dplyr::select(trait, result, note)
  
  df_traits_filtered <- rbind(df_traits_filtered, df_lowP_vec)
  cat(paste(sprintf("Removing traits with no variant having p < %.3e",pval_cutoff),
            paste(traits_removed,
                  collapse="\n"),sep=":\n"))
  z_mat <- z_mat[, minP_vec < pval_cutoff]

  # Prune traits by correlation (remove traits with Pearson |r| > 0.85)
  
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
  
  df_traits_filtered <- df_traits_filtered %>%
    right_join(df_traits, by="trait") %>%
    mutate(result = ifelse(is.na(result), "trait kept", result))

  pruned_traits <- df_traits_filtered %>%
    filter(result!="trait kept") %>%
    pull(trait)
  cat(paste("Traits removed in pruning process:", 
              paste(pruned_traits, collapse="\n"),sep = "\n"))
  
  cat(paste0("\nNumber of remaining traits: ",length(keep_traits),"\n"))
  
  # cat(paste("Remaing traits:", 
  #           paste(keep_traits, collapse="\n"),sep = "\n")) 
  print(paste0(sum(is.na(z_mat)), " missing values before corr pruning."))
  z_mat <- z_mat[, keep_traits]

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
  
  # Expand into N x 2M non-negative matrix
  print("Expanding z-score matrix into non-negative matrix (N-variants x 2M-traits)...")
  z_mat_pos <- z_mat
  z_mat_pos[z_mat_pos < 0] <- 0
  colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
  z_mat_neg <- -z_mat
  z_mat_neg[z_mat_neg < 0] <- 0
  colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
  final_z_mat <- cbind(z_mat_pos, z_mat_neg)
  
  output <- list(final_z_mat=final_z_mat,
                 df_traits=df_traits_filtered)
}


fill_missing_zcores <- function(initial_zscore_matrices,
                                df_snps,
                                trait_ss_files,
                                trait_ss_size,
                                main_ss_filepath,
                                rsID_map_file,
                                method_fill="median",
                                method_proxy="TopLD",
                                population="EUR",
                                LDlink_token=NULL,
                                topLD_path=NULL) {
  
  z_mat <- initial_zscore_matrices$z_mat
  N_mat <- initial_zscore_matrices$N_mat
  
  # update z_matrix
  if (method_fill=="median"){
    
    print("Filling missing zscores with trait medians...")
    z_mat_filled <- apply(z_mat, 2, function(x)
      replace(x, is.na(x), median(x, na.rm = TRUE)))
    N_mat_filled <- apply(N_mat, 2, function(x)
      replace(x, is.na(x), median(x, na.rm = TRUE)))
    
  } else if (method_fill=="zero") {
    
    print("Filling in missing values with zero...")
    z_mat_filled <- apply(z_mat, 2, function(x)
      replace(x, is.na(x), 0))
    N_mat_filled <- apply(N_mat, 2, function(x)
      replace(x, is.na(x), median(x, na.rm = TRUE)))
    
  } else if (method_fill=="remove") {
    
    print("Removing SNPs with missing values...")
    z_mat_filled <- na.omit(z_mat)
    N_mat_filled <- N_mat[rownames(z_mat_filled),]
    
    
  } else if (method_fill %in% c("proxy","proxies")){

    print("Searching for proxies for missing values/traits...")
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
        
        df_covers_tmp <- choose_proxies(need_proxies = need_covers_df_tmp,
                                              topLD_path = topLD_path,
                                              LDlink_token = LDlink_token,
                                              rsID_map_file = rsID_map_file,
                                              trait_ss_files = trait_ss_files_covers,
                                              pruned_variants = df_snps,
                                              method = method_pop,
                                              population = cur_population,
                                              frac_nonmissing_num = 0.5,
                                              r2_num = 0.5
        )
        df_list[[i]] <- df_covers_tmp
        i <- i+1
      }

      df_covers <- do.call(rbind, df_list) %>%
        filter(!duplicated(proxy_VAR_ID))
      print(sprintf("%i total covers found for all populations...", nrow(df_covers)))
      
      } else {
    
      print(sprintf("Searching for cover proxies for %s...", population))
      df_covers <- choose_proxies(need_proxies = need_covers_df,
                                        topLD_path = topLD_path,
                                        LDlink_token = LDlink_token,                                 
                                        rsID_map_file = rsID_map_file,
                                        trait_ss_files = trait_ss_files_covers,
                                        pruned_variants = df_snps,
                                        method = method_proxy,
                                        population=population,
                                        frac_nonmissing_num = 0.5,
                                        r2_num = 0.5)
        }
    
    print(sprintf("%i cover proxies found", nrow(df_covers)))
  
    if (nrow(df_covers)>0) {
      
      print("Fetching summary stats for cover SNPs...")
      df_covers <- df_covers %>%
        mutate(Risk_Allele=NA, PVALUE=NA) %>%
        rename(orig_VAR_ID=VAR_ID, VAR_ID=proxy_VAR_ID)

      trait_ss_size_covers <- trait_ss_size[unique(na_ix$col)]
      
      initial_zscore_matrices_cover <- fetch_summary_stats(
        df_covers,
        main_ss_filepath,
        trait_ss_files_covers,
        trait_ss_size_covers,
        pval_cutoff=1
      )
      print("Finished getting summary stats!")
      
      z_mat_covers <- initial_zscore_matrices_cover$z_mat
      N_mat_covers <- initial_zscore_matrices_cover$N_mat
      
    
      # match ChrPos for orig SNPs and cover proxies
      df_covers_merged <- df_covers %>%
        mutate(ChrPos_orig = gsub("_", ":", str_before_nth(orig_VAR_ID, "_", 2))) %>%
        mutate(ChrPos_proxy = gsub("_", ":", str_before_nth(VAR_ID, "_", 2)))
      
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
  
    print("Filling in remaining missing values with trait median...")
    
    z_mat_filled <- apply(df_z_fixed, 2, function(x)
        replace(x, is.na(x), median(x, na.rm = TRUE)))
    
    print("Matching N_mat row names with z_mat and replacing missing values with median...")
    df_N_fixed <- df_N_fixed[rownames(z_mat_filled), ]
    N_mat_filled <- apply(df_N_fixed, 2, function(x)
      replace(x, is.na(x), median(x, na.rm = TRUE)))
  } else {
    stop("Provide input for method_fill! (proxy, median, zero, remove)")
  }

  list(z_mat=z_mat_filled, N_mat=N_mat_filled)
}
