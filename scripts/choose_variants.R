library(tidyverse)
library(data.table)


# CURRENT ASSUMPTIONS ABOUT FORMATTING:
# - Genome build is hg19/GrCh37
# - Summary statistic datasets are whitespace-delimited with columns: VAR_ID, BETA, SE, N_PH
# - Variant IDs are all of the format: CHR_POS_REF_ALT 


ld_pruning <- function(gwas_variants, rsID_map_file, r2=0.1) {
  
  # Given a data frame of original GWAS variants (VAR_ID) and p-values (PVALUE), 
  # prune to a set of independent variants based on some LD threshold
  # Leverage the LDlinkR package to fetch LD relationships for a set of input SNPs
  
  write(gwas_variants$VAR_ID, "all_gwas_varid.tmp")
  all_var_df <- fread(cmd=paste0("grep -wFf all_gwas_varid.tmp ",
                                 rsID_map_file),
                      header=F, col.names=c("VAR_ID", "rsID"),
                      data.table=F, stringsAsFactors=F) %>%
    separate(VAR_ID, into=c("CHR", "POS", "REF", "ALT"), 
             sep="_", remove=F) %>%
    inner_join(gwas_variants, by="VAR_ID") %>%
    arrange(PVALUE)  # This ordering is important for the pruning steps below!
  system("rm all_gwas_varid.tmp")
  
  pruned_vars <- c()
  
  for (i in 1:22) {
    print(paste0("Pruning chromosome ", i, "..."))
    var_df <- filter(all_var_df, CHR == i)
    
    if (nrow(var_df) == 0) {
      next
    } else if (nrow(var_df) == 1) {
      pruned_vars <- c(pruned_vars, var_df$rsID)
      next
    }
    
    ld_mat <- LDlinkR::LDmatrix(snps=var_df$rsID,
                                pop="CEU",
                                r2d="r2",
                                token="20ff5a8454d7")  ## This should be replaced by each user's own token (retrieve at: https://ldlink.nci.nih.gov/?tab=apiaccess)
    rownames(ld_mat) <- ld_mat$RS_number
    ld_mat$RS_number <- NULL
    ld_mat <- as.matrix(ld_mat)
    ld_mat <- ld_mat[rowSums(is.na(ld_mat)) != ncol(ld_mat),
                     colSums(is.na(ld_mat)) != nrow(ld_mat)]

    remaining_snps <- var_df$rsID
    
    while(length(remaining_snps) > 0) {
      pruned_vars <- c(pruned_vars, remaining_snps[1])
      if (remaining_snps[1] %in% rownames(ld_mat)) {
        remaining_snps <- setdiff(
          remaining_snps,
          rownames(ld_mat)[ld_mat[, remaining_snps[1]] >= r2]
        )
      } else {
        remaining_snps <- setdiff(remaining_snps, remaining_snps[1])
      }
    }
  }
  
  print(paste0(length(pruned_vars), " variants remain after pruning."))
  filter(all_var_df, rsID %in% pruned_vars)
}


count_traits_per_variant <- function(gwas_variants, ss_files) {

  # Given a vector of variants and a named vector of summary statistics files 
  # for traits to be clustered, output a vector of non-missing trait fractions
  # per variant
  
  print("Assessing variant missingness across traits...")
  
  variant_df_list <- lapply(1:length(ss_files), function(i) {
    print(paste0("...Reading ", names(ss_files)[i], "..."))
    fread(ss_files[i], data.table=F, stringsAsFactors=F) %>%
      filter(VAR_ID %in% gwas_variants)
  })
  variant_counts_df <- bind_rows(variant_df_list) %>%
    group_by(VAR_ID) %>%
    summarise(frac=n() / length(ss_files))
  variant_counts <- ifelse(
    gwas_variants %in% variant_counts_df$VAR_ID,
    variant_counts_df$frac[match(gwas_variants, variant_counts_df$VAR_ID)],  # If in counts data frame, take the non-missing fraction
    0  # If not in data frame, then the non-missing fraction is 0
  )
  setNames(variant_counts, gwas_variants)
}


find_variants_needing_proxies <- function(gwas_variant_df, var_nonmissingness,
                                          rsID_map_file) {
  
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
                        var_nonmissingness[VAR_ID] < 0.8]
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
  
  left_join(tibble(VAR_ID=need_proxies_varid), varid_rsid_map, by="VAR_ID")
}


choose_proxies <- function(need_proxies, 
                           tabix_path, ld_file, 
                           rsID_map_file, ss_files,
                           pruned_variants) {
  
  # Given a vector of variants (rsIDs) needing proxies
  # (from find_variants_needing_proxies) and an LD reference file,
  # output a data frame linking each variant to a data frame containing possible
  # proxies (variant ID + r^2 + alleles)
  # Criteria for eligibility:
  #   Not strand-ambiguous
  #   Trait fraction >= 80%
  #   r^2 >= 0.8 with the index variant
  # Choose based on first trait count, then r^2
  
  # First, run "/path/to/tabix /path/to/LDfile rsID_1 rsID_2 ..."
  system(paste0(tabix_path, " ", ld_file, " ",
                paste(need_proxies$rsID, collapse=" "),
                " > ld_ref.tmp"))
  
  proxy_df <- read_tsv("ld_ref.tmp", col_names=c("rsID", "proxy_data")) %>%
    separate_rows(proxy_data, sep=";") %>%
    separate(proxy_data, into=c("proxy_rsID", "r2", "D"), sep=",")
  
  write(proxy_df$proxy_rsID, "potential_proxies_rsid.tmp")
  potential_proxies_map <- fread(cmd=paste0("grep -wFf potential_proxies_rsid.tmp ",
                                            rsID_map_file),
                                 header=F, col.names=c("proxy_VAR_ID", "proxy_rsID"),
                                 data.table=F, stringsAsFactors=F)

  system("rm ld_ref.tmp potential_proxies_rsid.tmp")
  
  proxy_missingness <- count_traits_per_variant(
    potential_proxies_map$proxy_VAR_ID,
    ss_files
  )
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
      frac_nonmissing >= 0.8,  # Sufficient fraction of traits non-missing
      r2 >= 0.8  # Sufficient LD with the proxied variant
    ) %>%
    group_by(rsID) %>%
    arrange(desc(frac_nonmissing),
            desc(r2)) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  proxies_found <- final_proxy_df$rsID
  no_proxies_found <- setdiff(need_proxies$rsID, proxies_found)
  print("No proxies needed for ", 
        length(setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID)), 
        " variants.")
  print(paste0("Proxies found for ", length(proxies_found), " variants."))
  print(paste0("No adequate proxies found for ", length(no_proxies_found), 
               " variants."))
  if (length(no_proxies_found) > 0) {
    write(no_proxies_found, "no_proxies_found.txt")
    print("See no_proxies_found.txt for a list of these variants.")
  }
  
  final_variant_set <- c(
    setdiff(pruned_variants$VAR_ID, need_proxies$VAR_ID),  # Original pruned variants that don't need proxies
    final_proxy_df$proxy_VAR_ID  # Proxy variants fulfilling the necessary criteria
  )
  unique(final_variant_set) 
  # NOTE: the unique() above simply discards duplicate proxies 
  # (same proxy for multiple variants and/or proxy variant that is already a primary GWAS variant).
  # There may be a better way to deal with this.

}


