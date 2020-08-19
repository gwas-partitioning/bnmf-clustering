library(tidyverse)


diamante_gwas <- read_tsv(
  "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/GWAS_DIAMANTE_eu_UKBB_dv2/T2D/DATA/GWAS_DIAMANTE_eu_UKBB_dv2.T2D.1.txt",
) %>%
  select(varID=VAR_ID, p=P_VALUE, N=N_PH) %>%
  filter(p < 5e-8) %>%
  separate(varID, into=c("chr", "pos", "EA", "NEA"), sep="_", remove=F)

exchip_gwas <- read_tsv(
  "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/ExChip_ExTexT2D_dv1/T2D/DATA/ExChip_ExTexT2D_dv1.T2D.1.txt",
) %>%
  select(varID=VAR_ID, p=P_VALUE, N=Neff) %>%
  filter(p < 5e-8) %>%
  separate(varID, into=c("chr", "pos", "EA", "NEA"), sep="_", remove=F)
  
diagram_gwas <- read_tsv(
  "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/GWAS_DIAGRAM_eu_dv2/T2D/DATA/GWAS_DIAGRAM_eu_dv2.T2D.1.txt",
) %>%
  select(varID=VAR_ID, p=P_VALUE, N=N_PH) %>%
  filter(p < 5e-8) %>%
  separate(varID, into=c("chr", "pos", "EA", "NEA"), sep="_", remove=F)

wgs_got2d_gwas <- read_tsv(
  "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/WGS_GoT2Dimputed_dv1/T2D/DATA/WGS_GoT2Dimputed_dv1.T2D.1.txt",
) %>%
  select(varID=VAR_ID, p=P_VALUE, N=N_PH) %>%
  filter(p < 5e-8) %>%
  separate(varID, into=c("chr", "pos", "EA", "NEA"), sep="_", remove=F)

mahajan_gwas <- read_tsv(
  "/humgen/diabetes2/users/clairekim/Mahajan.NatGenet2018b.T2D.European.txt",
  col_types=cols(SNP="c")
) %>%
  select(varID=SNP, chr=Chr, pos=Pos, p=Pvalue, N=Neff, EA, NEA) %>%
  mutate(varID=paste(chr, pos, EA, NEA, sep="_")) %>%
  filter(p < 5e-8)

t2d_eur_gwas <- read_tsv(
  "/humgen/diabetes2/users/clairekim/T2D_European.BMIunadjusted.txt",
  col_types=cols(SNP="c")
) %>%
  select(varID=SNP, chr=CHR, pos=BP, p=Pvalue, N=Neff, 
         EA=EFFECT_ALLELE, NEA=OTHER_ALLELE) %>%
  mutate(varID=paste(chr, pos, EA, NEA, sep="_")) %>%
  filter(p < 5e-8)


choose_gwas_variants <- function(ss_df_list) {
  
  # Given a list of summary statistic data frames, choose variants to take
  # forward for clustering
  # Summary statistic data frames should contain the following fields:
  # chr, pos, N, p, EA, NEA
  
  # Filter out studies with N < 10k
  low_N <- sapply(ss_df_list, function(ss_df) median(ss_df$N) < 10000)
  ss_df_list <- ss_df_list[which(!low_N)]
  
  # Standardize column types for binding
  ss_df_list <- lapply(ss_df_list, function(df) {
    mutate(df, 
           chr=as.character(chr),
           pos=as.integer(pos))
  })
  
  # Bind summary stats from each GWAS and select variants
  do.call(bind_rows, c(ss_df_list, .id="gwas")) %>%
    group_by(chr, pos) %>%
    dplyr::slice(which.max(N)) %>%
    ungroup()
}

gw_variants_list <- list(
  diamante=diamante_gwas,
  exchip=exchip_gwas,
  diagram=diagram_gwas,
  wgs_got2d=wgs_got2d_gwas,
  mahajan=mahajan_gwas,
  t2d_eur=t2d_eur_gwas
)

initial_gwas_variants_df <- choose_gwas_variants(gw_variants_list)

# input_file_list <- list(
#   c(varID="VAR_ID", p="P_VALUE", N="N_PH"),
#   c(varID="VAR_ID", p="P_VALUE", N="Neff"),
#   c(varID="VAR_ID", p="P_VALUE", N="N_PH"),
#   c(varID="VAR_ID", p="P_VALUE", N="N_PH"),
#   c(varID="VAR_ID", p="Pvalue", n="Neff"),
#   c(varID="VAR_ID", p=)
# )
# 
# names(input_file_list) <- c(
#   "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/GWAS_DIAMANTE_eu_UKBB_dv2/T2D/DATA/GWAS_DIAMANTE_eu_UKBB_dv2.T2D.1.txt",
#   "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/ExChip_ExTexT2D_dv1/T2D/DATA/ExChip_ExTexT2D_dv1.T2D.1.txt",
#   "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/GWAS_DIAGRAM_eu_dv2/T2D/DATA/GWAS_DIAGRAM_eu_dv2.T2D.1.txt",
#   "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/WGS_GoT2Dimputed_dv1/T2D/DATA/WGS_GoT2Dimputed_dv1.T2D.1.txt",
#   "/humgen/diabetes2/users/clairekim/Mahajan.NatGenet2018b.T2D.European.txt",
#   "/humgen/diabetes2/users/clairekim/T2D_European.BMIunadjusted.txt"
# )


