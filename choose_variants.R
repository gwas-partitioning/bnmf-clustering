library(tidyverse)
library(data.table)
library(LDlinkR)

# CURRENT ASSUMPTIONS ABOUT FORMATTING:
# - Genome build is hg19/GrCh37
# - Summary statistic datasets are whitespace-delimited with columns: VAR_ID, BETA, SE, N_PH
# - Variant IDs are all of the format: CHR_POS_REF_ALT 


prune <- function(var_df,
                  my_token,
                  population="CEU",
                  r2=0.1){
  pruned_vars <- c()
  if (nrow(var_df) == 1) {
    pruned_vars <- c(pruned_vars, var_df$rsID)
  }
  else if (nrow(var_df) > 1){

    print("Generating LD matrix using LDlinkR...")
    ld_mat <- LDlinkR::LDmatrix(snps=var_df$rsID,
                                pop=population,
                                r2d="r2",
                                token=my_token)  ## This should be replaced by each user's own token (retrieve at: https://ldlink.nci.nih.gov/?tab=apiaccess)
    if (!is.null(ld_mat$X.)) {
      cat("Error within LDlinkR! Chromosome not pruned!")
    }
    
    ld_mat <- ld_mat %>% column_to_rownames('RS_number')
    ld_mat <- as.matrix(ld_mat)
    ld_mat <- ld_mat[rowSums(is.na(ld_mat)) != ncol(ld_mat),
                     colSums(is.na(ld_mat)) != nrow(ld_mat)]
    remaining_snps <- var_df$rsID
    while(length(remaining_snps) > 0) {
      if (remaining_snps[1] %in% rownames(ld_mat)) {
        # only add to pruned_vars if found in LD matrix
        pruned_vars <- c(pruned_vars, remaining_snps[1])
        
        # remove current SNP and any other SNPs in LD with it
        remaining_snps <- setdiff(
          remaining_snps,
          rownames(ld_mat)[ld_mat[, remaining_snps[1]] >= r2]
        )
      } 
      else {
        # if current SNP not found in LD matrix, move on to next SNP
        remaining_snps <- setdiff(remaining_snps, remaining_snps[1])
      }
    }
  }
  pruning_output <- list(pruned_vars, ld_mat)
  return(pruning_output)
}

ld_pruning <- function(gwas_variants, rsID_map_file, my_token,
                       population="CEU", r2=0.1) {
  
  # Given a data frame of original GWAS variants (VAR_ID) and p-values (PVALUE), 
  # prune to a set of independent variants based on some LD threshold
  # Leverage the LDlinkR package to fetch LD relationships for a set of input SNPs
  
  write(gwas_variants$VAR_ID, "all_gwas_varid.tmp")
  
  print("Grepping for VAR_IDs in rsID map...")
  all_var_df <- fread(cmd=paste0("grep -wFf all_gwas_varid.tmp ", rsID_map_file),
                      header=F, col.names=c("VAR_ID", "rsID"),
                      data.table=F, stringsAsFactors=F) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove=F) %>%
    inner_join(gwas_variants, by="VAR_ID") %>%
    arrange(PVALUE)  # This ordering is important for the pruning steps below!
  
  print(paste("Num. SNPs mapped to rsID:",nrow(all_var_df)))
  
  pruned_vars <- c()
  ld_mats <- c()
  
  for (i in 1:22) {
    start=Sys.time()
    
    var_df_chr <- filter(all_var_df, CHR == i)
    n_snps <- nrow(var_df_chr)
    
    n_split <- 1000  # LDlinkR has query size limit of 1000
    mat_name <- as.character(i)
    
    # if n_snps above query size limit, break up the pruning 
    if (between(n_snps, 1 ,n_split)) { 
      
      print(paste0("Pruning chromosome ", i, "..."))
      var_df <- var_df_chr
  
      prune_output_chr <- prune(var_df=var_df,
                         my_token=my_token,
                         population=population,
                         r2=r2)
      
      # store pruned vars and LD matrix
      pruned_vars <- c(pruned_vars, prune_output_chr[[1]])
      ld_mats[[mat_name]] <- prune_output_chr[[2]]
      
      print(sprintf("Chr%i pruned from %i to %i variants",
                    i, n_snps, length(prune_output_chr[[1]])))
    }
    else if (n_snps > n_split) {
      var_df_list <- split(var_df_chr, (seq(nrow(var_df_chr))-1) %/% n_split)
      pruned_vars_list <- c()
      for (j in 1:length(var_df_list)){
        print(sprintf("Pruning subset %i for chromosome %i...", j, i))
        var_df <- var_df_list[[j]]
        
        prune_output_chr_j <- prune(var_df=var_df,
                                  my_token=my_token,
                                  population=population,
                                  r2=r2,
                                  method=method)
        
        # store temp pruned vars and LD matrix
        pruned_vars_list <- c(pruned_vars_list, prune_output_chr_j[[1]])
        ld_mats[[paste(mat_name, j, sep="_")]] <- prune_output_chr_j[[2]]
        
      }
      print(sprintf("Final pruning for chromosome %i...", i))
      var_df <- filter(var_df_chr, rsID %in% pruned_vars_list)
      
      prune_output_chr <- prune(var_df=var_df,
                               my_token=my_token,
                               population=population,
                               r2=r2,
                               method=method)
      
      # store final pruned vars and LD matrix
      pruned_vars <- c(pruned_vars, prune_output_chr[[1]])
      ld_mats[[mat_name]] <- prune_output_chr[[2]]
      
      print(sprintf("Chr%i pruned from %i to %i variants",
                    i, n_snps, length(prune_output_chr[[1]])))
    }
    else {
      print(sprintf("No SNPs on Chr%i!", i))
    }
    end=Sys.time()
    print(end-start)
  }
  print(paste0(length(pruned_vars), " VARIANTS REMAIN AFTER ALL PRUNING!"))
  
  # final outputs
  df_pruned <- filter(all_var_df, rsID %in% pruned_vars)
  return(list(df_pruned, ld_mats))
}

ld_pruning_topLD_api <-  function(gwas_variants, api_path, rsID_map_file, r2=0.1, population="EUR") {
  write(gwas_variants$VAR_ID, "all_gwas_varid.tmp")

  print("Grepping for VAR_IDs in rsID map...")
  all_var_df <- fread(cmd=paste0("grep -wFf all_gwas_varid.tmp ", rsID_map_file),
                      header=F, col.names=c("VAR_ID", "rsID"),
                      data.table=F, stringsAsFactors=F) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove=F) %>%
    inner_join(gwas_variants, by="VAR_ID") %>%
    arrange(PVALUE)

  system("touch outputLD.txt")
  system("touch outputInfo.txt")
  pruned_vars <- c()
  for (i in 1:22) {
    tmp <- all_var_df %>%
      filter(CHR==i)
    if (nrow(tmp)>1) {
      print(sprintf("Getting LD for CHR %i (%i variants)...", i, nrow(tmp)))
      rsID_pairs <- combn(tmp$rsID, 2, FUN=paste, collapse=',')
      cat(sprintf("%i rsID combinations...\n\n", length(rsID_pairs)))
      write(rsID_pairs, "to_prune_rsIDs.tmp")
      
      system(sprintf("%s -thres 0.0 -pop %s -maf 0.01 -inFile to_prune_rsIDs.tmp -outputLD outputLD_temp.txt -outputInfo outputInfo_temp.txt", api_path, population))
      system("awk 'FNR>1' outputLD_temp.txt >> outputLD.txt")
      system("awk 'FNR>1' outputInfo_temp.txt >> outputInfo.txt")

      ld_mat <- fread("outputLD.txt",
                      stringsAsFactors = F, data.table = F) %>%
        subset(rsID1 %like% "rs" & rsID1 != "rsID1")
      
      info_mat <- fread("outputInfo.txt",
                        stringsAsFactors = F, data.table = F) %>%
        subset(rsID %in% all_var_df$rsID) %>%
        subset(!duplicated(rsID))
      
      remaining_snps <- tmp %>%
        filter(rsID %in% info_mat$rsID) %>%
        arrange(PVALUE) %>%
        pull(rsID)
    
    pruned_vars_chr <- c()
    all_snps <- unique(c(ld_mat$rsID2, ld_mat$rsID2))
    while(length(remaining_snps) > 0) {
      if (remaining_snps[1] %in% all_snps) {
        # add to pruned_vars list
        pruned_vars_chr <- c(pruned_vars_chr, remaining_snps[1])
        # find rows where current SNP is present and R2>threshold
        tmp_mat <- ld_mat %>%
          subset(rsID1 %in% remaining_snps[1] | rsID2 %in% remaining_snps[1]) %>%
          subset(R2>r2)
        # remove any SNPs found in those rows from remaining_snps (this should also remove current SNP)
        remaining_snps <- setdiff(remaining_snps, unique(c(tmp_mat$rsID1, tmp_mat$rsID2, remaining_snps[1])))
      } else{
        # do  keep SNPs that are not found
        remaining_snps <- setdiff(remaining_snps, remaining_snps[1])
        }
      }
  } else {
    pruned_vars_chr <- tmp$rsID
  }
  print(sprintf("CHR %i pruned from %i to %i SNPs...",i, nrow(tmp), length(pruned_vars_chr)))
  pruned_vars <- c(pruned_vars, pruned_vars_chr)
  print(length(pruned_vars))
  }
    
  pruned_variants <- all_var_df %>%
    filter(rsID %in% pruned_vars)
}

snp_clump <- function(df_snps,
                      id="VAR_ID",
                      window=500000,
                      chr=1:22,
                      pos_range=c(1,Inf)) {
  
  clumped_snps <- c()
  uniq_chr <- sort(as.integer(unique(df_snps$CHR)))
  
  for (i in uniq_chr) {
    tmp <- df_snps %>%
      filter(CHR==i) %>%
      arrange(PVALUE) %>%
      mutate(POS = as.integer(POS)) %>%
      data.frame()
    
    if (i %in% chr) {
      print(sprintf("Clumping Chr. %i...",i))
      
      do_clump <- tmp %>%
        filter(between(POS, pos_range[1], pos_range[2]))
      dont_clump <- tmp %>%
        filter(!between(POS, pos_range[1], pos_range[2]))
      clumped_snps <- c(clumped_snps, dont_clump[,id])
      print(sprintf("Clumping %i variants, not clumping %i variants...",
                    nrow(do_clump), nrow(dont_clump)))
      
      remaining_snps <- do_clump[,id]
      
      j=0
      while (length(remaining_snps)>0){
        clumped_snps <- c(clumped_snps, remaining_snps[1])
        
        cur_pos <- do_clump$POS[do_clump[id]==remaining_snps[1]]
        
        close_snps <- do_clump[abs(do_clump$POS-cur_pos)<=window, id]

        remaining_snps <- setdiff(remaining_snps, close_snps)

        j=j+1
      }
    num_clumped <- j + nrow(dont_clump)
    } else {
      print(sprintf("No clumping for Chr %i!",i))
      num_clumped=nrow(tmp)
      clumped_snps <- c(clumped_snps, tmp[,id])
    }
    cat(sprintf("Chr %i clumped from %i to %i SNPs\n\n",i, nrow(tmp), num_clumped))
  }
  print(sprintf("No. SNPs after clumping: %i",length(clumped_snps)))
  return(clumped_snps)
}

ld_pruning_SNP.clip <- function(df_snps,
                                pop,
                                my_token,
                                r2=0.1,
                                maf=0.01,
                                chr=1:22,
                                output_dir="./") {
  
  snp_clip_input <- df_snps %>% 
    separate(VAR_ID, into=c("CHR","POS","REF","ALT"),sep="_",remove = F) %>%
    mutate(ChrPos = paste0("chr",CHR,":",POS)) %>%
    arrange(PVALUE)
  
  df_clipped <- data.frame(ChrPos=as.character(),
                           rsID=as.character())
  for (i in chr){
    start = Sys.time()
    cur_chr <- snp_clip_input %>%
      filter(CHR==i)
    print(sprintf("Chr %i (%i SNPs)",i, nrow(cur_chr)))

    
    if (nrow(cur_chr) == 0) {
      next
    }
    else if (nrow(cur_chr) == 1) {

    cur_snps <- cur_chr %>%
      pull(ChrPos)
    # if only one SNP, use LDhap to get the variant info
    clipped_res <- LDlinkR::LDhap(snps = cur_snps, 
          pop = pop, 
          token = my_token,
          genome_build = "grch37",
          table_type = "variant"
    ) %>%
      rename(Alleles=Allele_Frequency) %>%
      mutate(Details="Variant kept.")
  } else { # >1 SNP
    if (between(nrow(cur_chr), 1, 5000)) { 
      cur_snps <- cur_chr %>%
        pull(ChrPos)
      
    } else {
      print("Chromosome has >5000 SNPs; breaking into sections...")
      var_df_list <- split(cur_chr, (seq(nrow(cur_chr))-1) %/% 5000)
      cur_snps <- c()
      
      for (j in 1:length(var_df_list)){
        
        print(sprintf("Pruning subset %i for chromosome %i...", j, i))
        var_df <- var_df_list[[j]]
        cur_snps_j <- var_df %>%
          pull(ChrPos)
        
        clipped_res_split <- LDlinkR::SNPclip(
          cur_snps_j,
          pop = pop,
          r2_threshold = r2,
          maf_threshold = maf,
          token = my_token,
          file = FALSE,
          genome_build = "grch37")
        
        clipped_snps <- clipped_res_split %>%
          filter(Details=="Variant kept.") %>%
          pull(RS_Number)
        print(sprintf("Subset %i pruned to %i SNPs...", j, length(clipped_snps)))
        cur_snps <- c(cur_snps, clipped_snps)
      }
    }
    
    print(sprintf("Performing final chromosomal pruning for %i SNPs...", length(cur_snps)))
    clipped_res <- LDlinkR::SNPclip(
      cur_snps,
      pop = pop,
      r2_threshold = r2,
      maf_threshold = maf,
      token = my_token,
      file = FALSE,
      genome_build = "grch37")
  }
  
    fwrite(x=clipped_res,
           file = file.path(output_dir, sprintf("snpClip_results_%s_chr%i.txt", pop, i)),
           quote = F,
           sep = "\t")
    
    df_clipped_final <- clipped_res %>%
      filter(Details=="Variant kept.")
    print(sprintf("Chr%i pruned from %i to %i SNPs...",i, nrow(cur_chr), nrow(df_clipped_final)))
    
    end = Sys.time()
    print(end-start)
    
  }
  print("Done!")
  
}

count_traits_per_variant <- function(gwas_variants, ss_files) {

  # Given a vector of variants and a named vector of summary statistics files
  # for traits to be clustered, output a vector of non-missing trait fractions
  # per variant
  
  print("Assessing variant missingness across traits...")
  write(gwas_variants, "all_snps_varids.tmp")
  
  rename_cols <- c(N_PH="N")
  
  variant_df_list <- lapply(1:length(ss_files), function(i) {
    print(sprintf("...Reading %s...", names(ss_files)[i]))
    
    headers <- as.character(fread(ss_files[i], nrows=1,
                                  data.table=F, stringsAsFactors=F, header=F))
    
    if (endsWith(ss_files[i],".gz")) {
      df <- fread(cmd=sprintf("gzip -cd %s | fgrep -wf all_snps_varids.tmp ",ss_files[i]),
                  header=F,
                  col.names=headers,
                  data.table=F,
                  stringsAsFactors=F) %>%
        rename(any_of(rename_cols))
        
    } else {
      df <- fread(cmd=sprintf("fgrep -wf all_snps_varids.tmp %s ",ss_files[i]),
                  header=F,
                  col.names=headers,
                  data.table=F,
                  stringsAsFactors=F) %>%
        rename(any_of(rename_cols)) 
        
    }
    print(nrow(df))
    return(df)
  })
  
  # make dataframe of Ns
  df_N <- variant_df_list %>%
    setNames(names(ss_files)) %>%
    bind_rows(.id="trait") %>%
    select(trait, VAR_ID, N_PH) %>%
    pivot_wider(names_from="trait", values_from="N_PH") %>%
    data.frame()
}


find_variants_needing_proxies <- function(gwas_variant_df, var_nonmissingness,
                                          rsID_map_file, missing_cutoff=0.8) {
  
  # Given a data frame containing GWAS variants and alleles as well as a vector
  # of trait missingness fractions per variant (from count_traits_per_variant),
  # output a vector of variants that need proxies
  # Criteria (any of the following):
  #   Strand-ambiguous (AT or GC)
  #   Multi-allelic
  #   Low-count (available in < 80% of traits)
  # rsID_map_file should point to a whitespace-delimited file with columns
  # corresponding to VAR_ID and rsID
  
  print("Choosing variants in need of proxies...")
  
  gwas_variant_df <- gwas_variant_df %>%
      separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"),
               sep="_", remove=F)
  
  need_proxies_varid <- with(gwas_variant_df, {
    strand_ambig <- VAR_ID[paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")]
    print(paste0("...", length(strand_ambig), " strand-ambiguous variants"))
    
    multi_allelic <- grep("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", VAR_ID, value=T)  # i.e. ALT allele has a comma
    print(paste0("...", length(multi_allelic), " multi-allelic variants"))
    
    low_cnt <- VAR_ID[!(VAR_ID %in% names(var_nonmissingness)) |
                        var_nonmissingness[VAR_ID] < missing_cutoff]
    print(paste0("...", length(low_cnt), " variants with excessive missingness"))

    unique(c(strand_ambig, multi_allelic, low_cnt)) 
  })
  print(paste0("...", length(need_proxies_varid), " unique variants in total"))
  
  if (length(need_proxies_varid) == 0) return(tibble(VAR_ID=c(), rsID=c()))
  
  write(need_proxies_varid, "need_proxies_varid.tmp")
  varid_rsid_map <- fread(cmd=paste0("grep -wFf need_proxies_varid.tmp ",
                                     rsID_map_file),
                          header=F, col.names=c("VAR_ID", "rsID"),
                          data.table=F, stringsAsFactors=F)
  need_proxies_rsid <- varid_rsid_map$rsID[match(need_proxies_varid, 
                                                 varid_rsid_map$VAR_ID)]
  print(paste0("...", length(unique(varid_rsid_map$rsID)), 
               " of these are mapped to rsIDs"))
  system("rm need_proxies_varid.tmp")
  
  tibble(VAR_ID=need_proxies_varid) %>%
    left_join(varid_rsid_map, by="VAR_ID") %>%
    left_join(gwas_variant_df[,c("VAR_ID","PVALUE")], by="VAR_ID")
}


choose_proxies <- function(need_proxies, 
                           rsID_map_file,
                           trait_ss_files,
                           pruned_variants,
                           method="TopLD",
                           LDlink_token=NULL,
                           topLD_path=NULL,
                           population="EUR",
                           frac_nonmissing_num=0.8,
                           r2_num=0.8) {
  
  # Given a vector of variants (rsIDs) needing proxies
  # (from find_variants_needing_proxies) and an LD reference file,
  # output a data frame linking each variant to a data frame containing possible
  # proxies (variant ID + r^2 + alleles)
  # Criteria for eligibility:
  #   Not strand-ambiguous
  #   Trait fraction >= 80%
  #   r^2 >= 0.8 with the index variant
  # Choose based on first trait count, then r^2
  
  # First, run "/path/to/tabix /path/to/LDfile rsID_1 rsID_2 ...
  print(paste("Num rows need_proxies:",nrow(need_proxies)))
  if (method %in% c("LDlink","LDlinkR","LDproxy")) {
    
    print("Using LDlinkR:LDproxy_batch to find proxies...")
    need_proxies <- need_proxies %>%
      separate(VAR_ID, into=c("CHR","POS","REF","ALT"),sep = "_",remove = F) %>%
      mutate(query_snp = paste0("chr", CHR, ":", POS)) %>%
      select(-c(CHR, POS, REF, ALT))
    need_proxies_snps <- need_proxies$query_snp

   LDlinkR::LDproxy_batch(need_proxies_snps,
                          pop = population,
                          r2d = "r2",
                          token = LDlink_token,
                          append = T,
                          genome_build = "grch37")
   proxy_df <- read.table("./combined_query_snp_list_grch37.txt",sep = "\t",row.names = NULL) %>%
     filter(R2>r2_num) %>%
     filter(!Coord %in% need_proxies_snps) %>%
     inner_join(need_proxies, by = "query_snp") %>%
     arrange(PVALUE) %>%
     filter(!duplicated(RS_Number)) %>%
     dplyr::select(rsID, proxy_rsID=RS_Number, r2=R2) 
   need_proxies <- need_proxies %>%
     select(-c(query_snp))
    
  } else if (method=="TopLD") { # use TopLD
    print(sprintf("Using TopLD to find proxies for %s!", population))
    if (nrow(need_proxies)<100) {
      write(need_proxies$rsID, "need_proxies_rsIDs.tmp")
      system(sprintf("%s -thres %.1f -pop %s -maf 0.01 -inFile need_proxies_rsIDs.tmp -outputLD outputLD.txt -outputInfo outputInfo.txt", topLD_path, r2_num, population))
    } else { # need to split up
      print("Splitting proxy df into subsets (more than 100 SNPs)...")
      proxy_df_list <- split(need_proxies, (seq(nrow(need_proxies))-1) %/% 100)
      
      system("touch outputLD.txt")
      print("Running TopLD for proxy df segments...")
      for (j in 1:length(proxy_df_list)){
        print(sprintf("Querying LD subset %i/%i",j, length(proxy_df_list)))
        df <- proxy_df_list[[j]]
        write(df$rsID, "need_proxies_rsIDs.tmp")
        system(sprintf("%s -thres %.1f -pop %s -maf 0.01 -inFile need_proxies_rsIDs.tmp -outputLD outputLD_temp.txt -outputInfo outputInfo.txt", topLD_path, r2_num, population))
        system("cat outputLD_temp.txt >> outputLD.txt")
        }
      }
    proxy_df <- fread("outputLD.txt", stringsAsFactors = F, data.table = F) %>%
      select(rsID=rsID1, proxy_rsID=rsID2, r2=R2) %>%
      subset(proxy_rsID %like% "rs")
  } 
  else {
    stop("Enter appropriate proxy search method!") # Using stop function
    
  }
  print(paste("No. possible proxies found:",nrow(proxy_df))) # proxy_df should have columns (rsID, proxy_rsID, r2)
  write(proxy_df$proxy_rsID, "potential_proxies_rsid.tmp")
  
  if (nrow(proxy_df)>0) {
    print("Creating proxy rsID map...")
    potential_proxies_map <- fread(cmd=paste0("grep -wFf potential_proxies_rsid.tmp ",
                                              rsID_map_file),
                                   header=F, col.names=c("proxy_VAR_ID", "proxy_rsID"),
                                   data.table=F, stringsAsFactors=F)
    print(head(potential_proxies_map))
    
    proxy_variants <- potential_proxies_map$proxy_VAR_ID
    
    proxy_missingness <- count_traits_per_variant(
      proxy_variants,
      trait_ss_files
    )
    
    # get proxy missingness
    df_Ns_rev <- proxy_missingness %>%
      column_to_rownames("VAR_ID")
    df_Ns_rev[df_Ns_rev == 'NULL'] <- NA
  
    # get variant counts
    variant_counts_df <- data.frame(VAR_ID=rownames(df_Ns_rev),
                                    frac=rowSums(!is.na(df_Ns_rev))/length(trait_ss_files))

    proxy_missingness <- ifelse(
      proxy_variants %in% variant_counts_df$VAR_ID,
      variant_counts_df$frac[match(proxy_variants, variant_counts_df$VAR_ID)],  # If in counts data frame, take the non-missing fraction
      0  # If not in data frame, then the non-missing fraction is 0
    )
    proxy_missingness <- setNames(proxy_missingness, proxy_variants)
    
    proxy_missingness_df <- tibble(
      proxy_VAR_ID=names(proxy_missingness),
      frac_nonmissing=proxy_missingness
    )

    final_proxy_df <- proxy_df %>%
      inner_join(potential_proxies_map, by="proxy_rsID") %>%
      separate(proxy_VAR_ID, into=c("CHR", "POS", "REF", "ALT"), 
               sep="_", remove=F) %>%
      inner_join(proxy_missingness_df, by="proxy_VAR_ID") %>%
      filter(
        !(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")),  # Not strand-ambiguous
        !grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", proxy_VAR_ID),  # Not multi-allelic
        frac_nonmissing >= frac_nonmissing_num,  # Sufficient fraction of traits non-missing
        r2 >= r2_num  # Sufficient LD with the proxied variant
      ) %>%
      group_by(rsID) %>%
      arrange(desc(frac_nonmissing),
              desc(r2),
              CHR) %>%  # Arbitrary sort for reproducibility in case of missingness + r2 ties
      dplyr::slice(1) %>%
      ungroup() %>%
      inner_join(need_proxies, by="rsID") %>%  # added to include orig VAR_ID in output 
      data.frame()
    } else {
        final_proxy_df <- NULL
  }
  proxies_found <- final_proxy_df$rsID
  
  no_proxies_found <- setdiff(need_proxies$rsID, proxies_found)
  print(paste0("No proxies needed for ", 
               length(setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID)), 
               " variants."))
  print(paste0("Proxies found for ", length(proxies_found), " variants."))
  print(paste0("No adequate proxies found for ", length(no_proxies_found), 
               " variants."))
  
  return(final_proxy_df)

}

library(gwasrapidd)
use_gwas_catalog <- function(trait, my_pvalue=5e-8, min_N=10000, ancestry=NULL, rsID_map_file=NULL){
  # get all studies for breast cancer
  my_studies <- get_studies(efo_trait = trait)
  my_studies_big <- my_studies@studies %>% dplyr::filter(initial_sample_size>min_N)
  print(paste("# studies with > 10k samples:", nrow(my_studies_big)))
  
  if (!is.null(ancestry)){
    ancestry_studies <- my_studies@countries_of_recruitment %>%
      dplyr::filter(grepl(ancestry, region)) %>%
      pull(study_id)
    my_studies_big <- my_studies_big %>%
      dplyr::filter(study_id %in% ancestry_studies)
  }
  
  # how many studies
  # gwasrapidd::n(my_studies)
  print(paste("# studies with > 10k samples and in region of interest:", nrow(my_studies_big)))
  
  # get all associations from above studies
  # You could have also used get_associations(efo_trait = 'autoimmune disease')
  # my_associations <- get_associations(study_id = my_studies@studies$study_id,)
  my_associations_big <- get_associations(study_id = my_studies_big$study_id)
  
  # number of associations
  # gwasrapidd::n(my_associations) 
  print(paste("# associations:", gwasrapidd::n(my_associations_big)))
  
  # filter by p-value
  # Get association ids for which pvalue is less than 5e-8
  final_ids <- dplyr::filter(my_associations_big@associations, pvalue < my_pvalue) %>% # Filter by p-value
    tidyr::drop_na(pvalue) %>%
    dplyr::pull(association_id) 
  
  # Extract associations by association id
  final_associations <- my_associations_big[final_ids]
  # print(paste("# associations after pvalue filter:", gwasrapidd::n(final_associations)))
  
  # print var_IDs, risk alleles and risk_freq
  print("Getting final associations...")
  df_gwasRapidd <- final_associations@risk_alleles[c('association_id','variant_id', 'risk_allele')] %>%
    merge(final_associations@associations[c('association_id','pvalue')],by='association_id') %>%
    drop_na() %>%
    dplyr::select(rsID=variant_id, RiskAllele=risk_allele, PVALUE=pvalue) %>%
    filter(!duplicated(rsID))
  print(head(df_gwasRapidd))
  
  
  if (!is.null(rsID_map_file)){
    print("Grepping for rsIDs in rsID-to-VAR_ID map...")
    write(df_gwasRapidd$rsID,"gwas_rapidd.tmp")
    df_gwasRapidd <- fread(cmd=paste0("grep -wFf gwas_rapidd.tmp ", rsID_map_file),
                      header=F, col.names=c("VAR_ID", "rsID"),
                      data.table=F, stringsAsFactors=F) %>%
      merge(df_gwasRapidd, by="rsID") %>%
      dplyr::select(VAR_ID, PVALUE)
    print(head(df_gwasRapidd))
    }
  return(df_gwasRapidd)
}

