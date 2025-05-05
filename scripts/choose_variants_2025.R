library(tidyverse)
library(data.table)
library(LDlinkR)
library(tidyverse)
library(data.table)
library(future.apply)

# CURRENT ASSUMPTIONS ABOUT FORMATTING:
# - Genome build is hg19/GrCh37
# - Summary statistic datasets are whitespace-delimited with columns: VAR_ID, BETA, SE, N_PH
# - Variant IDs are all of the format: CHR_POS_REF_ALT 

# Set up a parallel plan – adjust based on your OS and desired number of workers
plan(multisession, workers = 2)

get_sig_snps <- function(gwas, PVCUTOFF = 5e-8, rename_cols = NULL) {
  # Prepare an empty list to store per-population results
  vars_sig_list <- list()
  my_gwas_populations <- unique(gwas$population)
  
  for (cur_pop in my_gwas_populations) {
    
    message(sprintf("Pulling significant SNPs from %s GWAS...", cur_pop))
    
    # Filter metadata for current population
    gwas_pop <- gwas %>% filter(population == cur_pop)
    gwas_ss_files_pop <- setNames(gwas_pop$full_path, gwas_pop$ID)
    
    # Process each file in parallel
    vars_pop_list <- lapply(seq_along(gwas_ss_files_pop), function(i) {
      cur_id <- names(gwas_ss_files_pop)[i]
      message(sprintf("...Reading %s...", cur_id))
      
      # Read summary stats and rename columns if necessary
      vars <- fread(gwas_ss_files_pop[i],
                    stringsAsFactors = FALSE,
                    data.table = FALSE) %>%
        dplyr::rename(any_of(rename_cols))
      
      # If BETA is missing, convert OR to log(OR)
      if (!"BETA" %in% colnames(vars)) {
        message("Converting Odds Ratio to Log Odds Ratio...")
        vars <- vars %>%
          mutate(BETA = log(as.numeric(ODDS_RATIO)))
      }
      
      # Process variant data:
      # 1. Convert BETA and P_VALUE to numeric.
      # 2. Filter well-formed VAR_ID and by p-value cutoff.
      # 3. Separate VAR_ID once (keeping REF and ALT).
      # 4. Compute Risk_Allele and add metadata.
      vars <- vars %>%
        mutate(across(c(BETA, P_VALUE), as.numeric)) %>%
        filter(grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]", VAR_ID),
               P_VALUE <= PVCUTOFF) %>%
        separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
        mutate(Risk_Allele = if_else(BETA >= 0, ALT, REF),
               Population = cur_pop,
               GWAS = cur_id) %>%
        select(VAR_ID, P_VALUE, Risk_Allele, GWAS, Population, REF, ALT)
      
      message(sprintf("...%s: %i significant SNPs...", cur_id, nrow(vars)))
      return(vars)
    })
    
    message(sprintf("Combining GWAS summary stats for %s...", cur_pop))
    vars_pop <- bind_rows(vars_pop_list)
    
    # Remove duplicates: keep the record with the smallest p-value per VAR_ID
    vars_pop <- vars_pop %>%
      arrange(VAR_ID, P_VALUE) %>%
      distinct(VAR_ID, .keep_all = TRUE) %>%
      dplyr::rename(PVALUE = P_VALUE)
    
    # Remove indels using REF and ALT columns
    vars_pop <- vars_pop %>%
      mutate(alleles = paste0(REF, ALT)) %>%
      filter(!(nchar(alleles) > 2 & !str_detect(ALT, ","))) %>%
      # Include REF and ALT in the final output if desired.
      select(VAR_ID, PVALUE, Risk_Allele, GWAS, Population, REF, ALT)
    
    message(sprintf("Final # SNPs for %s: %i", cur_pop, nrow(vars_pop)))
    
    vars_sig_list[[cur_pop]] <- vars_pop
  }
  
  # Combine results from all populations
  vars_sig <- bind_rows(vars_sig_list)
  message(sprintf("Final total # SNPs: %i", nrow(vars_sig)))
  
  return(vars_sig)
}


get_biggest_gwas <- function(main_ss_filepath, vars_sig) {
  
  # Ensure that main_ss_filepath exists
  if (!file.exists(main_ss_filepath)) {
    stop("The primary summary stats file does not exist at: ", main_ss_filepath)
  }
  
  # Prepare a variant data frame with separated columns and unique VAR_IDs
  vars_sig2 <- vars_sig %>%
    separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
    mutate(ChrPos = paste0("chr", CHR, ":", POS)) %>%
    arrange(PVALUE) %>%
    distinct(VAR_ID, .keep_all = TRUE)
  
  message("Getting p-values from largest GWAS...")
  
  # Write the list of VAR_IDs to a temporary file
  tmp_file <- "all_varid.tmp"
  write_lines(vars_sig2$VAR_ID, tmp_file)
  
  # Retrieve header from the compressed file.  
  # (This assumes that the first line of the file is a header.)
  headers <- fread(cmd = sprintf("gzip -cd %s | head -n1", main_ss_filepath),
                   header = FALSE,
                   data.table = FALSE,
                   stringsAsFactors = FALSE)
  colnames(headers) <- NULL  # Just to ensure no confusion if headers come in as a data frame
  header_line <- as.character(headers[1, ])  # use the header line for column names
  
  # Grep for VAR_IDs in the GWAS file using fgrep
  message("Grepping for VAR_IDs...")
  rename_cols <- c(PVALUE="P_VALUE")
                   
                     
  main_gwas <- fread(cmd = sprintf("gzip -cd %s | fgrep -wf %s", main_ss_filepath, tmp_file),
                        header = FALSE,
                        data.table = FALSE,
                        stringsAsFactors = FALSE,
                        col.names = header_line) %>%
    dplyr::rename(any_of(rename_cols))
  
  # Merge the p-values from MVP into our variant data frame
  vars_sig_with_big_gwas <- vars_sig2 %>%
    dplyr::rename(PVALUE.Pop = PVALUE) %>%
    inner_join(main_gwas[, c("VAR_ID", "PVALUE")], by = "VAR_ID") %>%
    mutate(PVALUE = as.numeric(PVALUE),
           PVALUE = if_else(PVALUE == 0, 1e-300, PVALUE))
  
  message(sprintf("%i of %i vars_sig found in primary GWAS...",
                  nrow(vars_sig_with_big_gwas), nrow(vars_sig2)))
  return(vars_sig_with_big_gwas)
  
}


prune <- function(var_df, my_token, population="CEU",
                  r2=0.1, method="curl"){
  pruned_vars <- c()
  if (nrow(var_df) == 1) {
    pruned_vars <- c(pruned_vars, var_df$rsID)
  }
  else if (nrow(var_df) > 1){
    if (method=="curl") {
      print(sprintf("Generating LD matrix using LDlinkRest (%i variants)...", nrow(var_df)))
      my_snps <- paste(var_df$rsID, collapse="%0A")
      
      my_snps <- paste(cur_snps, collapse="%0A")
      # my_snps <- paste(cur_snps, collapse="\n")

      ld_mat <- fread(cmd=sprintf("curl -k -X GET 'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?snps=%s&pop=%s&r2d=r2&token=%s'",
                          my_snps,
                          population,
                          my_token)
      )

    } else if (method=="LDlinkR") {
      print("Generating LD matrix using LDlinkR...")
      ld_mat <- LDlinkR::LDmatrix(snps=var_df$rsID,
                                  pop=population,
                                  r2d="r2",
                                  token=my_token)  ## This should be replaced by each user's own token (retrieve at: https://ldlink.nci.nih.gov/?tab=apiaccess)
      if (!is.null(ld_mat$X.)) {
        cat("Error within LDlinkR! Chromosome not pruned!")
      }
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
                       method="curl", population="CEU", r2=0.1) {
  
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
  
  print(head(all_var_df))
  
  bad_rsIDs <- subset(all_var_df, !startsWith(rsID,"rs"))

  print(paste("Num. SNPs mapped to rsID:",nrow(all_var_df)))
  
  # system("rm all_gwas_varid.tmp")
  
  pruned_vars <- c()
  ld_mats <- c()
  
  for (i in 1:22) {
    start=Sys.time()
    
    var_df_chr <- filter(all_var_df, CHR == i)
    n_snps <- nrow(var_df_chr)
    
    n_split <- ifelse(method=="curl", 250, 1000)
    mat_name <- as.character(i)
    
    if (between(n_snps, 1 ,n_split)) {  # LDlinkR has query size limit of 1000
      
      print(paste0("Pruning chromosome ", i, "..."))
      var_df <- var_df_chr
  
      prune_output_chr <- prune(var_df=var_df,
                         my_token=my_token,
                         population=population,
                         r2=r2,
                         method=method)
      
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
                      id="rsID",
                      window=500000,
                      chr=c(1:22),
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

library(tidyverse)
library(data.table)
library(LDlinkR)


# plan(multisession, workers = 4)

ld_pruning_SNP.clip <- function(df_snps,
                                pop,
                                r2 = 0.1,
                                maf = 0.01,
                                chr = 1:22,
                                output_dir = "./",
                                token = "cb5457b210a6") {
  
  # Check that required column exists
  if (!"VAR_ID" %in% colnames(df_snps)) {
    stop("df_snps must contain a VAR_ID column.")
  }
  
  # Prepare the SNP input data: separate VAR_ID and add ChrPos
  snp_clip_input <- df_snps %>% 
    separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
    mutate(ChrPos = paste0("chr", CHR, ":", POS)) %>%
    arrange(PVALUE)
  
  # Function to process a single chromosome
  process_chromosome <- function(i) {
    start_time <- Sys.time()
    out_file <- file.path(output_dir, sprintf("snpClip_results_%s_chr%i.txt", pop, i))
    
    # Check if output already exists to resume work
    if (file.exists(out_file)) {
      message(sprintf("Output file for Chr %i exists; skipping processing...", i))
      return(fread(out_file))
    }
    
    message(sprintf("Processing Chr %i...", i))
    cur_chr <- snp_clip_input %>% filter(CHR == i)
    
    if (nrow(cur_chr) == 0) {
      message(sprintf("No SNPs found for Chr %i. Skipping.", i))
      return(NULL)
    }
    
    # Initialize variable to hold pruned SNP results
    clipped_res <- NULL
    
    if (nrow(cur_chr) == 1) {
      # If only one SNP, use LDhap
      cur_snps <- cur_chr %>% pull(ChrPos)
      clipped_res <- tryCatch({
        LDlinkR::LDhap(snps = cur_snps, 
                       pop = pop, 
                       token = token,
                       genome_build = "grch37",
                       table_type = "variant") %>%
          rename(Alleles = Allele_Frequency) %>%
          mutate(Details = "Variant kept.")
      }, error = function(e) {
        message(sprintf("Error processing Chr %i (LDhap): %s", i, e$message))
        return(NULL)
      })
    } else {
      # More than one SNP: check if splitting is needed
      if (nrow(cur_chr) <= 5000) {
        cur_snps <- cur_chr %>% pull(ChrPos)
      } else {
        message("Chromosome has >5000 SNPs; breaking into sections...")
        var_df_list <- split(cur_chr, (seq_len(nrow(cur_chr)) - 1) %/% 5000)
        cur_snps <- c()
        
        for (j in seq_along(var_df_list)) {
          message(sprintf("Pruning subset %i for Chr %i...", j, i))
          var_df <- var_df_list[[j]]
          cur_snps_j <- var_df %>% pull(ChrPos)
          
          res_split <- tryCatch({
            LDlinkR::SNPclip(
              snps = cur_snps_j,
              pop = pop,
              r2_threshold = r2,
              maf_threshold = maf,
              token = token,
              file = FALSE,
              genome_build = "grch37")
          }, error = function(e) {
            message(sprintf("Error in subset %i for Chr %i: %s", j, i, e$message))
            return(NULL)
          })
          
          if (!is.null(res_split)) {
            clipped_snps <- res_split %>%
              filter(Details == "Variant kept.") %>%
              pull(RS_Number)
            message(sprintf("Subset %i pruned to %i SNPs...", j, length(clipped_snps)))
            cur_snps <- c(cur_snps, clipped_snps)
          }
        }
      }
      
      # Final pruning step for the chromosome
      message(sprintf("Performing final pruning for Chr %i on %i SNPs...", i, length(cur_snps)))
      clipped_res <- tryCatch({
        LDlinkR::SNPclip(
          snps = cur_snps,
          pop = pop,
          r2_threshold = r2,
          maf_threshold = maf,
          token = token,
          file = FALSE,
          genome_build = "grch37")
      }, error = function(e) {
        message(sprintf("Final pruning error for Chr %i: %s", i, e$message))
        return(NULL)
      })
    }
    
    # If pruning succeeded, write the output to file and return the result
    if (!is.null(clipped_res)) {
      fwrite(x = clipped_res, file = out_file, quote = FALSE, sep = "\t")
      df_clipped_final <- clipped_res %>% filter(Details == "Variant kept.")
      message(sprintf("Chr %i pruned from %i to %i SNPs.", i, nrow(cur_chr), nrow(df_clipped_final)))
    } else {
      message(sprintf("No results for Chr %i due to errors.", i))
    }
    
    message(sprintf("Time elapsed for Chr %i: %s", i, difftime(Sys.time(), start_time, units = "secs")))
    return(clipped_res)
  }
  
  # Process chromosomes sequentially using lapply instead of future_lapply.
  results_list <- lapply(chr, FUN = process_chromosome)
  
  message("LD pruning complete!")
  return(bind_rows(results_list))
}




library(tidyverse)
library(data.table)
library(future)
library(future.apply)


# Set up a parallel plan – adjust the number of workers as needed
rename_cols <- c(N_PH="N", N_PH="Neff")

count_traits_per_variant <- function(gwas_variants, ss_files) {
  
  # Given a vector of variants and a named vector of summary statistics files
  # for traits to be clustered, output a vector of non-missing trait fractions
  # per variant
  
  print("Assessing variant missingness across traits...")
  write(gwas_variants, "all_snps_varids.tmp")
  
  rename_cols <- c(N_PH="N", N_PH="Neff")
  
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
    if (nrow(df)==0) {
      df <- data.frame(VAR_ID="No matches", N_PH=0)
    }
    df <- df %>%
      filter(VAR_ID %in% gwas_variants) %>%
      select(VAR_ID, N_PH)
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


count_traits_per_variant_2025 <- function(gwas_variants, ss_files) {
  # Set up a parallel plan: use the minimum of the number of files and 8 workers
  plan(multisession, workers = 4)
  
  message("Assessing variant missingness across traits...")
  
  # Write variant IDs to a temporary file (used by the external grep command)
  tmp_file <- "all_snps_varids.tmp"
  writeLines(gwas_variants, con = tmp_file)
  
  # Use future_lapply to process each summary stats file in parallel
  variant_df_list <- future_lapply(seq_along(ss_files), function(i) {
    trait_name <- names(ss_files)[i]
    # Read header from the file (assume header is on the first line)
    headers <- as.character(fread(ss_files[i], nrows = 1,
                                  data.table = FALSE, header = FALSE, stringsAsFactors = FALSE))
    
    # Build the system command:
    cmd <- if (endsWith(ss_files[i], ".gz")) {
      sprintf("gzip -cd %s | fgrep -wf %s", ss_files[i], tmp_file)
    } else {
      sprintf("fgrep -wf %s %s", tmp_file, ss_files[i])
    }
    
    # Run the command and capture output using fread
    df <- fread(cmd = cmd,
                header = FALSE,
                col.names = headers,
                data.table = FALSE,
                stringsAsFactors = FALSE)
    
    # If no matches are found, create a placeholder data frame
    if (nrow(df) == 0) {
      df <- data.frame(VAR_ID = "No matches", N_PH = 0, stringsAsFactors = FALSE)
    }
    
    # Filter to keep only variants in our gwas_variants vector
    df <- df %>%
      filter(VAR_ID %in% gwas_variants) %>%
      select(VAR_ID, N_PH)
    
    return(df)
  })
  
  # Combine results into one data frame
  df_N <- variant_df_list %>%
    setNames(names(ss_files)) %>%
    bind_rows(.id = "trait") %>%
    select(trait, VAR_ID, N_PH) %>%
    pivot_wider(names_from = "trait", values_from = "N_PH") %>%
    as.data.frame()
  
  file.remove(tmp_file)
  return(df_N)
}# Example usage:
# result <- count_traits_per_variant(gwas_variants = your_variant_vector, ss_files = your_ss_files_named_vector)


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
                           tabix_path,
                           ld_file, 
                           rsID_map_file,
                           trait_ss_files,
                           pruned_variants,
                           zmat_fullset,           # Matrix (or data frame) with rownames in "Chr:Pos" format
                           method = "TopLD",
                           population = "EUR",
                           frac_nonmissing_num = 0.8, 
                           r2_num = 0.8) {
  
  message(sprintf("Number of rows in need_proxies: %d", nrow(need_proxies)))
  
  # ----- Proxy search via LDlinkR (or similar) or TopLD -----
  if (method %in% c("1000G", "LDlink", "LDlinkR", "LDproxy")) {
    message("Using LDlinkR:LDproxy_batch to find proxies...")
    need_proxies <- need_proxies %>%
      # Separate VAR_ID into its components so we can build a query in "chr:pos" format.
      separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
      mutate(query_snp = paste0("chr", CHR, ":", POS)) %>%
      select(-c(CHR, POS, REF, ALT))
    need_proxies_snps <- need_proxies$query_snp
    

    print(packageVersion("LDlinkR"))
    
    LDlinkR::LDproxy_batch(need_proxies_snps,
                           pop = population,
                           r2d = "r2",
                           token = "cb5457b210a6",
                           append = TRUE,
                           genome_build = "grch37")
    

    proxy_df <- read.table("./combined_query_snp_list_grch37.txt",sep = "\t",row.names = NULL) %>%
      filter(R2 > r2_num) %>% 
      filter(!Coord %in% need_proxies_snps) %>%
      inner_join(need_proxies, by = "query_snp") %>%
      arrange(PVALUE) %>%  # assumes the LD proxy output includes a PVALUE column
      filter(!duplicated(RS_Number)) %>%
      select(rsID, proxy_rsID = RS_Number, r2 = R2)
    print(sprintf("Found %i potential proxies for the %i queried SNPs...", nrow(proxy_df), nrow(need_proxies)))
    
    need_proxies <- need_proxies %>% select(-query_snp)
    
  } else if (method == "TopLD") {
    message(sprintf("Using TopLD to find proxies for %s!", population))
    if (nrow(need_proxies) < 100) {
      write(need_proxies$rsID, "need_proxies_rsIDs.tmp")
      system(sprintf("%s -thres %.1f -pop %s -maf 0.01 -inFile need_proxies_rsIDs.tmp -outputLD outputLD.txt -outputInfo outputInfo.txt", 
                     tabix_path, r2_num, population))
    } else { 
      message("Splitting proxy queries into subsets (more than 100 SNPs)...")
      proxy_df_list <- split(need_proxies, (seq(nrow(need_proxies)) - 1) %/% 100)
      system("touch outputLD.txt")
      for (j in seq_along(proxy_df_list)) {
        message(sprintf("Querying LD subset %d/%d", j, length(proxy_df_list)))
        df_sub <- proxy_df_list[[j]]
        write(df_sub$rsID, "need_proxies_rsIDs.tmp")
        system(sprintf("%s -thres %.1f -pop %s -maf 0.01 -inFile need_proxies_rsIDs.tmp -outputLD outputLD_temp.txt -outputInfo outputInfo.txt", 
                       tabix_path, r2_num, population))
        system("cat outputLD_temp.txt >> outputLD.txt")
      }

    }
    proxy_df <- fread("outputLD.txt", stringsAsFactors = FALSE, data.table = FALSE) %>%
      select(rsID = rsID1, proxy_rsID = rsID2, r2 = R2) %>%
      subset(proxy_rsID %like% "rs")
  } else {
    stop("Enter an appropriate proxy search method!")
  }
  
  message(sprintf("Number of possible proxies found: %d", nrow(proxy_df)))
  write(proxy_df$proxy_rsID, "potential_proxies_rsid.tmp")
  
  # ----- Use the rsID map and restrict proxies to those in the full summary stats -----
  
  if (nrow(proxy_df) > 0) {
    message("Creating proxy rsID map...")
    potential_proxies_map <- fread(cmd = sprintf("grep -wFf potential_proxies_rsid.tmp %s", rsID_map_file),
                                   header = FALSE, col.names = c("proxy_VAR_ID", "proxy_rsID"),
                                   data.table = FALSE, stringsAsFactors = FALSE)
        
    # Extract genomic coordinate from the proxy_VAR_ID (assumed to be in VAR_ID format, e.g. "chr1_12345_A_G")
    potential_proxies_map <- potential_proxies_map %>%
      separate(proxy_VAR_ID, into = c("CHR", "POS", "REF", "ALT"), 
               sep = "_", remove = FALSE) %>%
      mutate(proxy_ChrPos = paste(CHR, POS, sep = ":")) %>%
      # Restrict to proxies that are present in the z-score matrix (our final set)
      filter(proxy_ChrPos %in% rownames(zmat_fullset))
    print(sprintf("%i of %i potential proxies in the full z-matrix...", nrow(potential_proxies_map), nrow(proxy_df)))
    
    # ----- Assess missingness using zmat_fullset -----
    z_mat <- as.matrix(zmat_fullset)
    proxy_missingness <- rowSums(!is.na(z_mat)) / ncol(z_mat)
    proxy_missingness_df <- data.frame(
      proxy_ChrPos = names(proxy_missingness),
      frac_nonmissing = proxy_missingness,
      stringsAsFactors = FALSE
    )
  
    
    # ----- Merge missingness info with the proxy map and filter proxies -----
    final_proxy_df <- proxy_df %>%
      inner_join(potential_proxies_map, by = "proxy_rsID") %>%
      inner_join(proxy_missingness_df, by = "proxy_ChrPos")
    print(sprintf("Of the %i potential proxies after merging...", nrow(final_proxy_df)))
    
    final_proxy_df <- final_proxy_df %>%
      filter(r2 >= r2_num)
    print(sprintf("%i SNPs after LD r2 filter...", nrow(final_proxy_df)))
    
    final_proxy_df <- final_proxy_df %>%
      filter(!(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")) &
               !grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", proxy_VAR_ID))
    print(sprintf("%i SNPs after ambiguous and multi-allelic filter...", nrow(final_proxy_df)))
    
    final_proxy_df <- final_proxy_df %>%
      filter(frac_nonmissing >= frac_nonmissing_num)
    print(sprintf("%i SNPs after nonmissing filter...", nrow(final_proxy_df)))
    
    print("Selecting best remaining proxy for each query SNP...")
    final_proxy_df <- final_proxy_df %>% 
      group_by(rsID) %>%
      arrange(desc(frac_nonmissing), desc(r2), CHR) %>%  # Sort to pick the best proxy for each index variant
      dplyr::slice(1) %>%
      ungroup() %>%
      inner_join(need_proxies, by = "rsID")  # Retain original variant info (e.g., Risk Allele) for downstream use
  } else {
    final_proxy_df <- NULL
  }
  
  proxies_found <- if (!is.null(final_proxy_df)) final_proxy_df$rsID else character(0)
  no_proxies_found <- setdiff(need_proxies$rsID, proxies_found)
  
  message(sprintf("No proxies needed for %d variants.", length(setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID))))
  message(sprintf("Proxies found for %d variants.", length(proxies_found)))
  message(sprintf("No adequate proxies found for %d variants.", length(no_proxies_found)))
  
  if (length(no_proxies_found) > 0) {
    write(no_proxies_found, "no_proxies_found.txt")
    message("See no_proxies_found.txt for a list of these variants.")
  }
  
  final_variant_set <- list(
    setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID),  # Original pruned variants that don't need proxies
    final_proxy_df                                        # Proxy variants fulfilling the criteria
  )
  
  return(final_variant_set)
}





suzuki_pruning <- function(assoc_data_filtered, my_LDlinkR_token, populations = c("EUR"), chunk_size = 2500, output_file = "ld_pruning_results.rds") {
  
  # Check if a previous run exists
  if (file.exists(output_file)) {
    message(sprintf("Resuming from previous run: %s", output_file))
    results <- readRDS(output_file)
    
    chr_NAs <- results %>% group_by(CHR) %>%
      summarize(n = sum(!is.na(clump))) %>%
      arrange(desc(n)) %>%
      dplyr::rename(non_NA = n)
    chr_counts <- results %>% dplyr::count(CHR) %>% dplyr::rename(snps=n)
    chr_checks <- chr_NAs %>% inner_join(chr_counts, by='CHR')
    processed_chr <- chr_checks %>%
      filter(snps==non_NA) %>%
      pull(CHR)
    
    assoc_data_sorted <- results
  } else {
    message("Starting new LD pruning run.")
    assoc_data_sorted <- assoc_data_filtered %>%
      arrange(P, CHR, BP) %>%
      mutate(clump = NA_integer_, index_snp = FALSE)
    processed_chr <- character(0) # Initialize empty character vector
  }
  
  clump_index <- max(assoc_data_sorted$clump, na.rm = TRUE) + 1
  if (!is.finite(clump_index)) clump_index <- 1 #set clump index to 1 if no clump has been assigned yet.
  unique_chr <- unique(assoc_data_sorted$CHR)
  
  for (chr in unique_chr) {
    if (chr %in% processed_chr) {
      message(sprintf("Chromosome %s already processed. Skipping.", chr))
      next
    }
    
    message(sprintf("Processing chromosome %s", chr))
    
    assoc_chr <- assoc_data_sorted %>% filter(CHR == chr)
    
    # Important: We need to reset these vectors based on the current state of assoc_data_sorted
    # This ensures we're working with the latest data
    snp_list <- assoc_chr$rsID
    bp_list <- assoc_chr$BP
    p_list <- assoc_chr$P
    
    if (length(snp_list) < 2) {
      message(sprintf("  Skipping chromosome %s due to less than 2 SNPs.", chr))
      processed_chr <- c(processed_chr, chr)
      saveRDS(assoc_data_sorted, file = output_file) # Save progress
      next
    }
    
    snp_chunks <- split(snp_list, ceiling(seq_along(snp_list) / chunk_size))
    bp_chunks <- split(bp_list, ceiling(seq_along(bp_list) / chunk_size))
    
    combined_ld_chunks <- list()
    
    for (chunk_idx in seq_along(snp_chunks)) {
      chunk_snps <- snp_chunks[[chunk_idx]]
      message(sprintf("  Processing chunk %d of %d", chunk_idx, length(snp_chunks)))
      
      ld_list <- list()
      for (pop in populations) {
        message(sprintf("    Calling LDmatrix for population %s", pop))
        ld_mat <- tryCatch({
          LDlinkR::LDmatrix(snps = chunk_snps,
                            pop = pop,
                            token = my_LDlinkR_token,
                            genome_build = "grch37")
        }, error = function(e) {
          message(sprintf("      Error for pop %s: %s", pop, e$message))
          return(NULL)
        })
        if (!is.null(ld_mat)) {
          mat <- as.matrix(ld_mat[, -1])
          rownames(mat) <- ld_mat$RS_number
          ld_list[[pop]] <- mat
        }
      }
      
      combined_ld_chunk <- Reduce(function(x, y) pmax(x, y, na.rm = TRUE), ld_list)
      combined_ld_chunks[[chunk_idx]] <- combined_ld_chunk
    }
    
    if (length(combined_ld_chunks) > 1) {
      all_snps_in_chr = unique(unlist(lapply(combined_ld_chunks, rownames)))
      combined_ld = matrix(NA, nrow=length(all_snps_in_chr), ncol=length(all_snps_in_chr), dimnames = list(all_snps_in_chr, all_snps_in_chr))
      
      for (ld_chunk in combined_ld_chunks){
        combined_ld[rownames(ld_chunk), colnames(ld_chunk)] = ld_chunk
      }
    } else if (length(combined_ld_chunks) == 1) {
      combined_ld = combined_ld_chunks[[1]]
    } else {
      # No LD chunks were created - likely due to errors
      message(sprintf("  No valid LD data retrieved for chromosome %s. Skipping.", chr))
      processed_chr <- c(processed_chr, chr)
      saveRDS(assoc_data_sorted, file = output_file) # Save progress
      next
    }
    
    # Process SNPs by p-value order (most significant first)
    chr_indices <- which(assoc_data_sorted$CHR == chr & is.na(assoc_data_sorted$clump))
    
    # Sort by P-value - crucial for getting the right index SNPs
    sorted_indices <- chr_indices[order(assoc_data_sorted$P[chr_indices])]
    
    for (idx in sorted_indices) {
      # Skip if this SNP has already been assigned to a clump
      if (!is.na(assoc_data_sorted$clump[idx])) next
      
      snp_i <- assoc_data_sorted$rsID[idx]
      current_bp <- assoc_data_sorted$BP[idx]
      
      # Mark this SNP as an index SNP and assign it to a new clump
      assoc_data_sorted$clump[idx] <- clump_index
      assoc_data_sorted$index_snp[idx] <- TRUE
      
      message(sprintf("    Created new clump %d with index SNP %s", clump_index, snp_i))
      
      # Find SNPs in LD with this index SNP
      # Look for unassigned SNPs within 5Mb
      candidate_indices <- which(
        assoc_data_sorted$CHR == chr &
          is.na(assoc_data_sorted$clump) &
          abs(assoc_data_sorted$BP - current_bp) <= 5e6
      )
      
      clump_size <- 1  # Start with 1 (the index SNP)
      
      for (j in candidate_indices) {
        snp_j <- assoc_data_sorted$rsID[j]
        
        # Check if SNPs are in LD
        if (!is.na(snp_i) && !is.na(snp_j) && 
            snp_i %in% rownames(combined_ld) && 
            snp_j %in% colnames(combined_ld)) {
          
          ld_value <- combined_ld[snp_i, snp_j]
          
          if (!is.na(ld_value) && ld_value > 0.05) {
            # Assign to the current clump but mark as NOT an index SNP
            assoc_data_sorted$clump[j] <- clump_index
            assoc_data_sorted$index_snp[j] <- FALSE
            clump_size <- clump_size + 1
          }
        }
      }
      
      message(sprintf("    Clump %d has %d SNPs in total", clump_index, clump_size))
      
      # Move to the next clump
      clump_index <- clump_index + 1
    }
    
    processed_chr <- c(processed_chr, chr)
    saveRDS(assoc_data_sorted, file = output_file) # Save progress
  }
  
  # Final summary
  clump_summary <- assoc_data_sorted %>%
    filter(!is.na(clump)) %>%
    group_by(clump) %>%
    summarize(
      n_snps = n(),
      has_index = any(index_snp),
      chr = CHR[1],
      index_snp = ifelse(any(index_snp), rsID[which(index_snp)[1]], NA)
    )
  
  n_clumps <- nrow(clump_summary)
  n_index_snps <- sum(clump_summary$has_index)
  
  cat(sprintf("Clumping complete. Processed %i chromosomes.\n", length(processed_chr)))
  cat(sprintf("Found %i clumps with %i index SNPs.\n", n_clumps, n_index_snps))
  
  return(assoc_data_sorted)
}