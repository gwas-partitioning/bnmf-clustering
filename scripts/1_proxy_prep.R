library(tidyverse)


yengo_top941_ss <- read_tsv(  # Read in initial set of 941 COJO variants from Yengo 2018
  paste0("../data/raw/sum_stats/yengo2018/", 
         "Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt")
) %>%
  rename(rsid=SNP)
write(yengo_top941_ss$rsid, file="../data/processed/yengo_941_rsIDs.txt")  # Write rsIDs to file

#### Code to be run interactively on UGER -- retrieve potential LD proxies for Yengo variants ####
# SOURCE_PATH=/humgen/diabetes2/users/mvg/VariantClustering
# rsids_space_delim="$(cat yengo_941_rsIDs.txt | tr '\n' ' ')"
# $SOURCE_PATH/tabix-0.2.6/tabix $SOURCE_PATH/LD_EUR.tsv.bgz $rsids_space_delim > yengo_LD_ref.txt  # LD ref. is from HaploReg v4, which uses 1KG


find_potential_proxies <- function(ld_file, id_map, ld_thresh, base_trait) {
  # Process and filter an LD reference panel file to generate a set of potential 
  # variant proxies
  # ld_file: from LD reference (formatted as in 1000 Genomes reference panel)
  #   rs123   rs1,0.86-0.95;rs2,0.5,-0.4;etc.
  # id_map: data frame with two columns (snpid=chr:pos:A2:A1, rsid)
  # ld_thresh: minimum LD with index variant required to keep a potential proxy
  # base_trait: trait whose variants are being proxied (for file naming)
  potential_proxies_df <- map_dfr(ld_file, function(sap) {  # Loop through index SNPs
    tibble(rsid=sap[1],  # Index SNP (sap = "snp and proxies")
           proxies=strsplit(sap[2], split=";")[[1]]) %>%  # Set of semicolon-delimited LD SNPs
      separate(proxies, c("proxy_rsid", "ld", "unsure"), ",")  # Get rsID and LD for each potential proxy
  })
  potential_proxies_filtered_df <- potential_proxies_df %>%
    left_join(id_map, by=c("proxy_rsid"="rsid")) %>%  # Join proxy rsIDs with positional IDs to add locus and allele info
    separate(snpid, c("chr", "pos", "A2", "A1"), "_", remove=F) %>%
    filter(ld >= ld_thresh,  # Retain only SNPs having LD > specified threshold with the index SNP
           !is.na(A1) & !is.na(A2),
           !(paste0(A2, A1) %in% c("AT", "TA", "CG", "GC"))) %>%  # Remove strand-ambiguous variants
    mutate(proxy_loc=paste(chr, pos, sep=":")) %>%  # Construct a chr:pos ID per SNP
    select(rsid, proxy_rsid, proxy_snpid=snpid, proxy_loc, proxy_ld=ld)
  write_tsv(potential_proxies_filtered_df, 
            paste0("../data/processed/", base_trait, "_potential_proxies.txt"))
  potential_proxies_filtered_df
}


all_snps_and_proxies <- strsplit(readLines("../data/processed/yengo_LD_ref.txt"), 
                                 split="\t")  # Result is list of length 941

snpid_rsid_map <- read_tsv("../data/raw/list_VARID_rsID.txt",  # From dbSNP v1.38 -- maps positional IDs to rsIDs
                           col_names=c("snpid", "rsid")) %>%
  distinct(rsid, .keep_all=T)  # Remove duplicate rsIDs

potential_proxies_filtered_df <- find_potential_proxies(all_snps_and_proxies, snpid_rsid_map, 0.6, "bmi")
