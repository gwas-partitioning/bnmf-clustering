# =============================================================================
# choose_variants_2025.R
# =============================================================================
# Helper functions for GWAS variant selection, LD clumping, and proxy search.
#
# Functions:
#   get_sig_snps()               — Extract genome-wide significant SNPs across GWAS files
#   get_biggest_gwas()           — Cross-reference hits against a large reference GWAS
#   snp_clump()                  — Position-based clumping (no external API required)
#   ld_pruning_SNP.clip()        — LD pruning via LDlink SNPclip (multi-population)
#   prune()                      — Single-chromosome LD pruning via LDlinkR or curl
#   ld_pruning()                 — Full-genome LD pruning via LDlinkR (legacy)
#   ld_pruning_topLD_api()       — LD pruning via TopLD local binary (legacy)
#   count_traits_per_variant()   — (deprecated) Count non-missing traits per variant
#   count_traits_per_variant_2025() — Parallel version of above
#   window_to_sentinels()        — Filter variant pool to within kb-window of sentinel variants
#   check_topmed_presence()      — Liftover hg19 variants to hg38 and confirm presence in TOPMed VCFs
#   find_variants_needing_proxies() — Flag strand-ambiguous, multiallelic, high-missingness, or non-TOPMed variants
#   choose_proxies()             — Search for LD proxies via LDlinkR::LDproxy_batch
#
# Assumptions:
#   - Genome build: hg19/GRCh37
#   - Variant IDs formatted as: CHR_POS_REF_ALT
#   - Summary stat files are whitespace-delimited with at minimum: VAR_ID, BETA, SE
# =============================================================================

library(tidyverse)
library(data.table)
library(LDlinkR)
library(future)
library(future.apply)

# Set up a parallel plan — adjust workers based on your system
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
  return(pruned_variants)
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

ld_pruning_SNP.clip <- function(df_snps,
                                pop,
                                r2 = 0.1,
                                maf = 0.01,
                                chr = 1:22,
                                output_dir = "./",
                                token,
                                parallel = FALSE) {

  if (parallel) library(future.apply)
  
  #— 0) PREP: ensure VAR_ID, split into fields, add ChrPos & order by P
  if (!"VAR_ID" %in% colnames(df_snps)) {
    stop("df_snps must contain a VAR_ID column.")
  }
  
  snp_clip_input <- df_snps %>%
    separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
    mutate(
      ChrPos = paste0("chr", CHR, ":", POS),
      CHR = as.integer(CHR),
      POS = as.integer(POS)
    ) %>%
    arrange(PVALUE)
  
  #— 1) Find which chr‐files are already done
  pattern <- sprintf("^snpClip_results_%s_chr(\\\\d+)\\.txt$", pop)
  done_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  done_chr <- if (length(done_files) > 0) {
    as.integer(sub(pattern, "\\1", basename(done_files)))
  } else integer(0)
  
  #— 2) Subset to only pending chromosomes
  pending_chr <- setdiff(chr, done_chr)
  if (length(pending_chr) == 0) {
    message("All population–chromosome pairs already completed. Reading results from disk...")
    out_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
    return(bind_rows(lapply(out_files, data.table::fread)))
  }
  message("Will process chromosomes: ", paste(pending_chr, collapse = ", "))
  
  #— 3) Single‐chromosome handler
  process_chromosome <- function(i) {
    st <- Sys.time()
    out_file <- file.path(output_dir, sprintf("snpClip_results_%s_chr%i.txt", pop, i))
    if (file.exists(out_file)) {
      message(" Chr ", i, " already done; skipping.")
      return(NULL)
    }
    message(" Processing Chr ", i, " …")
    cur_chr <- snp_clip_input %>% filter(CHR == i)
    if (nrow(cur_chr) == 0) {
      message("  No SNPs on Chr ", i, "; skipping.")
      return(NULL)
    }
    
    clipped_res <- NULL
    
    if (nrow(cur_chr) == 1) {
      # Single SNP: use LDhap
      clipped_res <- tryCatch({
        LDhap(
          snps = cur_chr$ChrPos,
          pop = pop,
          token = token,
          genome_build = "grch37",
          table_type = "variant"
        ) %>%
          rename(Alleles = Allele_Frequency) %>%
          mutate(Details = "Variant kept.")
      }, error = function(e) {
        message("LDhap error: ", e$message)
        NULL
      })
      
    } else if (nrow(cur_chr) <= 5000) {
      # Moderate size: prune directly
      clipped_res <- tryCatch({
        LDlinkR::SNPclip(
          snps = cur_chr$ChrPos,
          pop = pop,
          r2_threshold = r2,
          maf_threshold = maf,
          token = token,
          file = FALSE,
          genome_build = "grch37"
        )
      }, error = function(e) {
        message("SNPclip error: ", e$message)
        NULL
      })
      
    } else {
      # Large: split into chunks
      message(" Chromosome ", i, " has >5000 SNPs; breaking into sections...")
      var_df_list <- split(cur_chr, (seq(nrow(cur_chr)) - 1) %/% 5000)
      kept_snps <- c()
      
      for (j in seq_along(var_df_list)) {
        message(sprintf(" Pruning subset %i for chromosome %i...", j, i))
        var_df <- var_df_list[[j]]
        cur_snps_j <- var_df %>% pull(ChrPos)
        
        clipped_res_split <- tryCatch({
          LDlinkR::SNPclip(
            snps = cur_snps_j,
            pop = pop,
            r2_threshold = r2,
            maf_threshold = maf,
            token = token,
            file = FALSE,
            genome_build = "grch37"
          )
        }, error = function(e) {
          message("Subset SNPclip error: ", e$message)
          NULL
        })
        
        if (!is.null(clipped_res_split)) {
          kept_snps_subset <- clipped_res_split %>%
            filter(Details == "Variant kept.") %>%
            pull(RS_Number)
          kept_snps <- c(kept_snps, kept_snps_subset)
          message(sprintf(" Subset %i pruned to %i SNPs...", j, length(kept_snps_subset)))
        }
      }
      
      # Final pruning pass over combined kept SNPs
      message(sprintf(" Performing final chromosomal pruning for %i SNPs...", length(kept_snps)))
      clipped_res <- tryCatch({
        LDlinkR::SNPclip(
          snps = kept_snps,
          pop = pop,
          r2_threshold = r2,
          maf_threshold = maf,
          token = token,
          file = FALSE,
          genome_build = "grch37"
        )
      }, error = function(e) {
        message("Final SNPclip error: ", e$message)
        NULL
      })
    }
    
    #— Write out if successful
    if (!is.null(clipped_res)) {
      data.table::fwrite(clipped_res, file = out_file, sep = "\t", quote = FALSE)
      msg <- sprintf("  Chr %i: %i → %i SNPs", 
                     i, nrow(cur_chr), nrow(clipped_res %>% filter(Details == "Variant kept.")))
      message(msg)
    } else {
      message("  No results for Chr ", i)
    }
    message("  elapsed: ", round(difftime(Sys.time(), st, units = "secs"), 1), "s")
    return(NULL)
  }
  
  #— 4) Run chromosomes (parallel or not)
  if (parallel) {
    future.apply::future_lapply(pending_chr, process_chromosome)
  } else {
    lapply(pending_chr, process_chromosome)
  }
  
  #— 5) Read all results
  final_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  all_res <- bind_rows(lapply(final_files, data.table::fread))
  message("LD pruning complete across populations & chromosomes.")
  return(all_res)
}
count_traits_per_variant <- function(gwas_variants, ss_files) {
  # DEPRECATED: Use count_traits_per_variant_2025() for parallel execution.
  #
  # Given a vector of variants and a named vector of summary statistics files
  # for traits to be clustered, output a data frame of sample sizes per variant-trait.
  
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


# -----------------------------------------------------------------------------
#' Filter a variant pool to within a genomic window of sentinel variants
#'
#' Used to create a tractable proxy candidate pool for fetch_summary_stats():
#' instead of fetching all P < PVCUTOFF_PROXY variants genome-wide (which can
#' be 100k+ for polygenic traits), restrict to variants near LD-pruned sentinels.
#' Any proxy returned by LDlinkR will by definition be in this neighbourhood.
#'
#' @param candidates  data.frame with CHR (int/char) and POS (int) columns
#' @param sentinels   data.frame with CHR and POS columns (LD-pruned sentinel variants)
#' @param window_kb   numeric; half-window size in kilobases (default 500)
#' @return Subset of candidates (same class/columns) within window_kb of at
#'   least one sentinel. Rows are unique; order is preserved.
window_to_sentinels <- function(candidates, sentinels, window_kb = 500) {
  window_bp <- as.integer(round(window_kb * 1000))

  cand_dt <- copy(as.data.table(candidates))
  sent_dt <- as.data.table(sentinels)

  # Point intervals for candidates
  cand_dt[, `:=`(chr_key = as.character(CHR),
                 start   = as.integer(POS),
                 end     = as.integer(POS))]

  # Sentinel windows
  sent_windows <- sent_dt[, .(
    chr_key = as.character(CHR),
    start   = pmax(0L, as.integer(POS) - window_bp),
    end     = as.integer(POS) + window_bp
  )]

  setkey(sent_windows, chr_key, start, end)
  setkey(cand_dt,      chr_key, start, end)

  hits <- foverlaps(cand_dt, sent_windows, type = "any", nomatch = NULL)

  # Return unique rows from candidates (VAR_ID preferred; fall back to row index)
  if ("VAR_ID" %in% names(candidates)) {
    out <- candidates[candidates$VAR_ID %in% unique(hits$VAR_ID), ]
  } else {
    out <- candidates[candidates$ChrPos %in% unique(hits$ChrPos), ]
  }

  message(sprintf(
    "  window_to_sentinels: %d / %d candidates retained within %g kb of a sentinel",
    nrow(out), nrow(candidates), window_kb
  ))
  out
}


# -----------------------------------------------------------------------------
#' Check which hg19 variants are present in TOPMed (BRAVO VCF files)
#'
#' Lifts variant positions from hg19 to hg38 using a liftOver chain file, then
#' queries per-chromosome BRAVO VCF files (tabix-indexed) via Rsamtools to confirm
#' presence by position. Variants that fail liftover are treated as absent.
#' Called twice in the pipeline: once for pruned_vars (to flag "not_in_topmed")
#' and once for the proxy candidate pool (to restrict proxy selection to TOPMed).
#'
#' @param variants_hg19  data.frame with columns VAR_ID (CHR_POS_REF_ALT),
#'                       CHR (character or integer), POS (integer) — hg19 coordinates
#' @param chain_file     path to hg19ToHg38.over.chain file
#' @param vcf_dir        directory containing chr*.bravo.pub.vcf.gz (+ .tbi index) files
#' @return character vector of VAR_IDs confirmed present in TOPMed at the
#'         lifted-over hg38 position; absent or unliftable variants are excluded
check_topmed_presence <- function(variants_hg19, chain_file, vcf_dir) {

  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("rtracklayer required. Install via: BiocManager::install('rtracklayer')")
  if (!requireNamespace("Rsamtools", quietly = TRUE))
    stop("Rsamtools required. Install via: BiocManager::install('Rsamtools')")

  message(sprintf("  [TOPMed] Lifting over %d variants (hg19 -> hg38)...", nrow(variants_hg19)))
  chain <- rtracklayer::import.chain(chain_file)

  gr_hg19 <- GenomicRanges::GRanges(
    seqnames = paste0("chr", variants_hg19$CHR),
    ranges   = IRanges::IRanges(start = as.integer(variants_hg19$POS),
                                end   = as.integer(variants_hg19$POS)),
    VAR_ID   = variants_hg19$VAR_ID
  )

  lifted_list <- rtracklayer::liftOver(gr_hg19, chain)

  # Discard variants that map to 0 or multiple hg38 positions
  one_to_one <- lengths(lifted_list) == 1L
  n_dropped  <- sum(!one_to_one)
  if (n_dropped > 0)
    message(sprintf("  [TOPMed] %d variant(s) failed liftover — treated as absent.", n_dropped))

  gr_hg38 <- unlist(lifted_list[one_to_one])

  df_hg38 <- data.frame(
    VAR_ID  = gr_hg38$VAR_ID,
    chr_num = sub("^chr", "", as.character(GenomicRanges::seqnames(gr_hg38))),
    pos     = GenomicRanges::start(gr_hg38),
    stringsAsFactors = FALSE
  )

  confirmed <- character(0)

  for (chr_num in unique(df_hg38$chr_num)) {
    vcf_path <- file.path(vcf_dir, sprintf("chr%s.bravo.pub.vcf.gz", chr_num))

    if (!file.exists(vcf_path)) {
      message(sprintf("  [TOPMed] VCF not found: chr%s.bravo.pub.vcf.gz — %d variant(s) marked absent.",
                      chr_num, sum(df_hg38$chr_num == chr_num)))
      next
    }

    chr_df <- df_hg38[df_hg38$chr_num == chr_num, ]

    ranges <- GenomicRanges::GRanges(
      paste0("chr", chr_num),
      IRanges::IRanges(start = chr_df$pos, end = chr_df$pos)
    )

    tbx  <- Rsamtools::TabixFile(vcf_path)
    hits <- Rsamtools::scanTabix(tbx, param = ranges)

    # hits: named list, one element per queried range; non-empty = position found in VCF
    found      <- vapply(hits, function(x) length(x) > 0L, logical(1L))
    confirmed  <- c(confirmed, chr_df$VAR_ID[found])
  }

  message(sprintf("  [TOPMed] %d / %d variants confirmed present in TOPMed.",
                  length(confirmed), nrow(variants_hg19)))
  confirmed
}


find_variants_needing_proxies <- function(gwas_variant_df,
                                          var_nonmissingness,
                                          missing_cutoff = 0.8,
                                          topmed_fails   = NULL) {

  print("Choosing variants in need of proxies...")

  # Track reasons for needing proxies
  proxy_reasons_list <- list()

  # Find strand-ambiguous variants
  strand_ambig <- with(gwas_variant_df, VAR_ID[paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")])
  print(paste0("...", length(strand_ambig), " strand-ambiguous variants"))
  if (length(strand_ambig) > 0) {
    proxy_reasons_list[["ambiguous"]] <- data.frame(VAR_ID = strand_ambig, reason = "ambiguous", stringsAsFactors = FALSE)
  }

  # Find multi-allelic variants
  multi_allelic <- grep("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", gwas_variant_df$VAR_ID, value = TRUE)
  print(paste0("...", length(multi_allelic), " multi-allelic variants"))
  if (length(multi_allelic) > 0) {
    proxy_reasons_list[["multiallelic"]] <- data.frame(VAR_ID = multi_allelic, reason = "multiallelic", stringsAsFactors = FALSE)
  }

  # Find variants with excessive missingness
  low_cnt <- with(gwas_variant_df, VAR_ID[!(VAR_ID %in% names(var_nonmissingness)) |
                                            var_nonmissingness[VAR_ID] < missing_cutoff])
  print(paste0("...", length(low_cnt), " variants with excessive missingness"))
  if (length(low_cnt) > 0) {
    proxy_reasons_list[["high_missingness"]] <- data.frame(VAR_ID = low_cnt, reason = "high_missingness", stringsAsFactors = FALSE)
  }

  # Find variants absent from TOPMed (optional — only when topmed_fails is provided)
  if (!is.null(topmed_fails)) {
    not_topmed <- intersect(gwas_variant_df$VAR_ID, topmed_fails)
    print(paste0("...", length(not_topmed), " variants not found in TOPMed"))
    if (length(not_topmed) > 0) {
      proxy_reasons_list[["not_in_topmed"]] <- data.frame(VAR_ID = not_topmed, reason = "not_in_topmed", stringsAsFactors = FALSE)
    }
  }
  
  # Combine all reasons and handle overlaps
  if (length(proxy_reasons_list) > 0) {
    all_reasons <- bind_rows(proxy_reasons_list)
    
    # For variants with multiple reasons, combine them
    reasons_summary <- all_reasons %>%
      group_by(VAR_ID) %>%
      summarise(reason = paste(reason, collapse = ","), .groups = 'drop')
    
    need_proxies_varid <- unique(all_reasons$VAR_ID)
  } else {
    need_proxies_varid <- character(0)
    reasons_summary <- data.frame(VAR_ID = character(0), reason = character(0))
  }
  
  print(paste0("...", length(need_proxies_varid), " unique variants in total"))
  
  # Simply return the variants with their reasons and existing rsIDs
  result <- gwas_variant_df %>%
    filter(VAR_ID %in% need_proxies_varid) %>%
    left_join(reasons_summary, by = "VAR_ID") %>%
    select(VAR_ID, rsID = RS_Number, reason, PVALUE)
  
  return(result)
}

choose_proxies <- function(need_proxies,
                           rsid_map_dir = "rsid_maps_by_chr",
                           pruned_variants,
                           zmat_fullset,           # matrix/data frame with rownames in "CHR:POS" format
                           token,                  # LDlink API token (Sys.getenv("LDLINK_TOKEN"))
                           population = "EUR",
                           frac_nonmissing_num = 0.8,
                           r2_num = 0.8,
                           topmed_present_snps = NULL) {  # optional: hg19 ChrPos (e.g. "1:12345") confirmed in TOPMed

  message(sprintf("Number of rows in need_proxies: %d", nrow(need_proxies)))

  # Initialize candidate storage
  candidates_with_data <- NULL

  # ----- Proxy search via LDlinkR::LDproxy_batch -----
  message(sprintf("Using LDlinkR to find proxies for %s population...", population))

  need_proxies <- need_proxies %>%
    separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
    mutate(query_snp = paste0("chr", CHR, ":", POS)) %>%
    select(-c(CHR, POS, REF, ALT))
  need_proxies_snps <- need_proxies$query_snp

  LDlinkR::LDproxy_batch(need_proxies_snps,
                         pop            = population,
                         r2d            = "r2",
                         token          = token,
                         append         = TRUE,
                         genome_build   = "grch37")

  # LDproxy_batch writes combined_query_snp_list_grch37.txt to the working directory
  proxy_out_file <- "./combined_query_snp_list_grch37.txt"
  if (file.exists(proxy_out_file)) {
    proxy_df <- read.table(proxy_out_file, sep = "\t", header = TRUE, row.names = NULL) %>%
      filter(R2 > 0.1) %>%
      filter(!Coord %in% need_proxies_snps) %>%
      inner_join(need_proxies, by = "query_snp") %>%   # join on query_snp col present in LDproxy batch output
      filter(!duplicated(RS_Number)) %>%
      select(rsID, proxy_rsID = RS_Number, r2 = R2)
  } else {
    proxy_df <- data.frame(rsID = character(), proxy_rsID = character(), r2 = numeric())
  }

  message(sprintf("Found %i potential proxies for the %i queried SNPs.", nrow(proxy_df), nrow(need_proxies)))
  need_proxies <- need_proxies %>% select(-query_snp)
  
  message(sprintf("Number of possible proxies found: %d", nrow(proxy_df)))
  
  # ----- Use chromosome-specific rsID maps to find proxy information -----
  
  if (nrow(proxy_df) > 0) {
    message("Creating proxy rsID map using chromosome-specific files...")
    
    proxy_rsids <- unique(proxy_df$proxy_rsID)
    potential_proxies_list <- list()
    chr_files <- list.files(rsid_map_dir, pattern = "chr.*\\.txt$", full.names = TRUE)
    
    if (length(chr_files) == 0) stop(sprintf("No chromosome files found in %s", rsid_map_dir))
    
    for (chr_file in chr_files) {
      chr_name <- basename(chr_file)
      # message(sprintf("Searching %s...", chr_name))
      chr_rsid_map <- fread(chr_file, 
                            col.names = c("hg19_posID", "proxy_rsID", "ref_allele", "alt_allele"),
                            key = "proxy_rsID")
      chr_matches <- chr_rsid_map[proxy_rsID %in% proxy_rsids]
      if (nrow(chr_matches) > 0) potential_proxies_list[[chr_name]] <- chr_matches
    }
    
    # Combine results from all chromosomes
    if (length(potential_proxies_list) > 0) {
      potential_proxies_map <- rbindlist(potential_proxies_list)
      potential_proxies_map <- as.data.frame(potential_proxies_map)
      
      # Process Map: Extract Coords and Fix Orientation
      potential_proxies_map <- potential_proxies_map %>%
        separate(hg19_posID, into = c("CHR", "POS"), sep = ":", remove = FALSE) %>%
        mutate(CHR = gsub("chr", "", CHR), 
               proxy_VAR_ID_orig = paste(CHR, POS, ref_allele, alt_allele, sep = "_"),
               proxy_VAR_ID_flip = paste(CHR, POS, alt_allele, ref_allele, sep = "_"),
               proxy_ChrPos = paste(CHR, POS, sep = ":")) %>%
        filter(proxy_ChrPos %in% rownames(zmat_fullset)) %>%
        mutate(
          proxy_VAR_ID = case_when(
            proxy_VAR_ID_orig %in% rownames(zmat_fullset) ~ proxy_VAR_ID_orig, 
            proxy_VAR_ID_flip %in% rownames(zmat_fullset) ~ proxy_VAR_ID_flip,
            TRUE ~ proxy_VAR_ID_orig 
          ),
          REF_final = case_when(proxy_VAR_ID == proxy_VAR_ID_flip ~ alt_allele, TRUE ~ ref_allele),
          ALT_final = case_when(proxy_VAR_ID == proxy_VAR_ID_flip ~ ref_allele, TRUE ~ alt_allele)
        ) %>%
        select(-proxy_VAR_ID_orig, -proxy_VAR_ID_flip) %>%
        dplyr::rename(REF = REF_final, ALT = ALT_final) %>%
        select(-ref_allele, -alt_allele)
      
      print(sprintf("%i of %i potential proxies in the full z-matrix...", nrow(potential_proxies_map), nrow(proxy_df)))

      # ----- TOPMed filter (optional) -----
      # Restrict all proxy candidates to those confirmed present in TOPMed so that
      # the final variant set contains only TOPMed-genotyped variants (better PRS coverage).
      if (!is.null(topmed_present_snps) && nrow(potential_proxies_map) > 0) {
        n_before <- nrow(potential_proxies_map)
        potential_proxies_map <- potential_proxies_map %>%
          filter(proxy_ChrPos %in% topmed_present_snps)
        message(sprintf("  TOPMed filter: %d -> %d proxy candidates (removed %d not in TOPMed)",
                        n_before, nrow(potential_proxies_map), n_before - nrow(potential_proxies_map)))
      }

    } else {
      potential_proxies_map <- data.frame()
    }

    if (nrow(potential_proxies_map) > 0) {
      # ----- Assess missingness -----
      z_mat <- as.matrix(zmat_fullset)
      proxy_missingness <- rowSums(!is.na(z_mat)) / ncol(z_mat)
      proxy_missingness_df <- data.frame(
        proxy_ChrPos = names(proxy_missingness),
        frac_nonmissing = proxy_missingness,
        stringsAsFactors = FALSE
      )
      
      # ----- SAVE CANDIDATES WITH DATA (Before Filtering) -----
      candidates_with_data <- proxy_df %>%
        inner_join(potential_proxies_map, by = "proxy_rsID") %>%
        inner_join(proxy_missingness_df, by = "proxy_ChrPos")
      
      # ----- Apply Strict Filters for "Success" -----
      final_proxy_df <- candidates_with_data %>%
        filter(r2 >= r2_num) %>%
        filter(!(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")) &
                 !grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", proxy_VAR_ID)) %>%
        filter(frac_nonmissing >= frac_nonmissing_num)
      
      # Select best proxy
      final_proxy_df <- final_proxy_df %>% 
        group_by(rsID) %>%
        arrange(desc(frac_nonmissing), desc(r2), CHR) %>% 
        dplyr::slice(1) %>%
        ungroup() %>%
        inner_join(need_proxies, by = "rsID") 
      
    } else {
      final_proxy_df <- NULL
    }
  } else {
    final_proxy_df <- NULL
  }
  
  # Calculate Stats
  proxies_found <- if (!is.null(final_proxy_df)) final_proxy_df$rsID else character(0)
  no_proxies_found <- setdiff(need_proxies$rsID, proxies_found)
  
  message(sprintf("No proxies needed for %d variants.", length(setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID))))
  message(sprintf("Proxies found for %d variants.", length(proxies_found)))
  message(sprintf("No adequate proxies found for %d variants.", length(no_proxies_found)))
  
  if (length(no_proxies_found) > 0) {
    write(no_proxies_found, "no_proxies_found.txt")
    message("See no_proxies_found.txt for a list of these variants.")
  }

  # ==========================================================
  # FINAL SELECTION LOGIC: FALLBACK WITH GUARDRAILS (FIXED)
  # ==========================================================
  
  # 1. Identify who got a proxy
  if (!is.null(final_proxy_df) && nrow(final_proxy_df) > 0) {
    replaced_originals <- unique(final_proxy_df$VAR_ID)
  } else {
    replaced_originals <- character(0)
  }
  
  # 2. Identify who FAILED to get a proxy
  failed_search_vars <- setdiff(need_proxies$VAR_ID, replaced_originals)
  
  # 3. Assess "Rescue" Potential for Failed Variants
  rescued_vars <- character(0)
  
  if (length(failed_search_vars) > 0) {
    message("Assessing failed variants for rescue (Strict QC)...")
    
    check_df <- data.frame(VAR_ID = failed_search_vars, stringsAsFactors = FALSE) %>%
      separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
      mutate(
        ChrPos_prefix = paste0("chr", CHR, ":", POS),  # Format: chr1:100
        ChrPos_noprefix = paste0(CHR, ":", POS)        # Format: 1:100
      )
    
    # Calculate missingness using BOTH formats to be safe
    if (!is.null(zmat_fullset)) {
      z_mat <- as.matrix(zmat_fullset)
      
      # Check which format matches the Z-matrix keys
      keys_prefix <- intersect(check_df$ChrPos_prefix, rownames(z_mat))
      keys_noprefix <- intersect(check_df$ChrPos_noprefix, rownames(z_mat))
      
      # Decide which key column to use
      if (length(keys_prefix) >= length(keys_noprefix)) {
        use_col <- "ChrPos_prefix"
        valid_keys <- keys_prefix
      } else {
        use_col <- "ChrPos_noprefix"
        valid_keys <- keys_noprefix
      }
      # message(sprintf("Matching Z-matrix using format: %s (Matches: %d)", use_col, length(valid_keys)))
      
      # Calculate fractions
      miss_stats <- rowSums(!is.na(z_mat[valid_keys, , drop=FALSE])) / ncol(z_mat)
      
      check_df$frac_nonmissing <- 0 
      match_idx <- match(check_df[[use_col]], names(miss_stats))
      check_df$frac_nonmissing[!is.na(match_idx)] <- miss_stats[na.omit(match_idx)]
      
    } else {
      warning("zmat_fullset missing; cannot calculate missingness for rescue. Defaulting to 0.")
      check_df$frac_nonmissing <- 0
    }
    
    # APPLY THE FILTERS
    rescue_candidates <- check_df %>%
      mutate(
        is_ambiguous = paste0(REF, ALT) %in% c("AT", "TA", "GC", "CG"),
        is_multiallelic = grepl(",", ALT),
        is_missing_too_high = frac_nonmissing < 0.5
      ) 
    
    # Corrected Logging: Use 'rescue_candidates' (the mutated df), not 'check_df'
    message(sprintf("  - Dropped %d ambiguous variants", sum(rescue_candidates$is_ambiguous)))
    message(sprintf("  - Dropped %d multi-allelic variants", sum(rescue_candidates$is_multiallelic)))
    message(sprintf("  - Dropped %d variants with <50%% data", sum(rescue_candidates$is_missing_too_high)))
    
    # Final Selection
    rescued_final <- rescue_candidates %>%
      filter(!is_ambiguous & !is_multiallelic & !is_missing_too_high)
    
    rescued_vars <- rescued_final$VAR_ID
    message(sprintf("  - RESCUED %d variants (Clean & >50%% data)", length(rescued_vars)))
  }
  
  # 4. Construct Final Sets
  vars_never_needed_proxy <- setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID)
  final_originals_to_keep <- c(vars_never_needed_proxy, rescued_vars)
  
  final_variant_set <- list(
    final_originals_to_keep,  # 1. Originals
    final_proxy_df           # 2. Proxies
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