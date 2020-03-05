library(tidyverse)


### This script will generally require manually-specified functions to accomodate
### the varying summary statistic formats from different GWAS. For each trait
### or study, a single .RData file should be saved. This will contain two objects:
### a named list of summary statistic data frames (for all index and potential 
### proxy variants, with fields: rsid, A1, A2, z) and a sample size data frame
### (fields: trait, N).

################################################################################

# Specify:
# index_ss: Data frame containing index variants (rsID and chr:pos ID)
# potential_proxies: Data frame containing potential proxy IDs and LD information (from 1_proxy_prep.R)

yengo_top941_ss <- read_tsv(  # Read in initial set of 941 COJO variants from Yengo 2018
  paste0("../data/raw/sum_stats/yengo2018/",
         "Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt")
) %>%
  mutate(loc=paste(CHR, POS, sep=":")) %>%
  rename(rsid=SNP)
index_ss <- yengo_top941_ss
potential_proxies <- read_tsv("../data/processed/yengo_potential_proxies.txt",
                              col_types="ccccd")

################################################################################

#### Assemble the set of variants to be retrieved for each trait ####

keep_snps <- bind_rows(select(index_ss, rsid, loc),  # Data frame with positional and rsIDs for all index and potential proxy SNPs
                       select(potential_proxies, rsid=proxy_rsid, loc=proxy_loc))


#### Gather and preprocess summary statistics for each group of traits ####

# YENGO 2018 -- LOCKE + UKB META-ANALYSIS FOR BMI
process_yengo_bmi <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/yengo2018/", ss_file)) %>%
    mutate(z=BETA / SE) %>%
    select(rsid=SNP, A1=Tested_Allele, A2=Other_Allele, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
yengo_bmi_files <- list(bmi="Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
yengo_bmi_ss <- lapply(yengo_bmi_files, process_yengo_bmi)
yengo_bmi_N <- tibble(trait="bmi", N=716500)  # Based on mean from sum stats
save("yengo_bmi_ss", "yengo_bmi_N", 
     file="../data/processed/cleaned_sum_stats/yengo_bmi.RData")

# PULIT 2018 -- LOCKE + UKB META-ANALYSIS FOR BODY FAT DISTRIBUTION
process_pulit_whr <- function(ss_file) {
  read_delim(paste0("../data/raw/sum_stats/pulit2018/", ss_file), delim=" ") %>%
    mutate(z=BETA / SE,
           rsid=gsub(":.*", "", SNP),
           loc=paste(CHR, POS, sep=":")) %>%
    select(rsid, A1=Tested_Allele, A2=Other_Allele, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
pulit_whr_files <- list(
  whr_adjBMI="Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
)
pulit_whr_ss <- lapply(pulit_whr_files, process_pulit_whr)
pulit_whr_N <- tibble(trait="whr_adjBMI", N=694700)
save("pulit_whr_ss", "pulit_whr_N", 
     file="../data/processed/cleaned_sum_stats/pulit_whr.RData")

# YANG 2012 -- BMI VARIABILITY
process_yang_variability <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/yang2012/", ss_file)) %>%
    mutate(z=b / se,
           A1=toupper(Allele1), 
           A2=toupper(Allele2)) %>%
    select(rsid=MarkerName, A1, A2, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
yang_variability_files <- list(
  bmi_variability="GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_BMI.txt.gz"
)
yang_variability_ss <- lapply(yang_variability_files, process_yang_variability)
yang_variability_N <- tibble(trait="bmi_variability", N=133200)  # From publication 
save("yang_variability_ss", "yang_variability_N", 
     file="../data/processed/cleaned_sum_stats/yang_variability.RData")

# LEAN BODY MASS (ZILLIKENS 2018; GEFOS)
process_lbm <- function(ss_file) {
  read_table2(paste0("../data/raw/sum_stats/gefos/", ss_file)) %>%
    mutate(z=Effect / StdErr,
           A1=toupper(Allele1),
           A2=toupper(Allele2)) %>%
    select(rsid=MarkerName, A1, A2, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
lbm_files <- list(lbm="wholebodyleanmass.results.metal_.txt.gz")
lbm_ss <- lapply(lbm_files, process_lbm)
lbm_N <- tibble(trait="lbm", N=37500)  # From s.s.
save("lbm_ss", "lbm_N", file="../data/processed/cleaned_sum_stats/lbm.RData")

# BONE AREA (STYRKARSDOTTIR 2019)
process_ba <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/styrkarsdottir2019/", ss_file)) %>%
    mutate(z=sign(BetaA1) * qnorm(GCcorrP / 2, lower.tail=F),
           loc=gsub("chr", "", paste(Chrom, Pos, sep=":"))) %>%
    select(rsid=rsids, loc, A1, A2=A0, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
ba_files <- list(
  hip_bone_area="Hip_Total_Area.txt.gz"
)
ba_ss <- lapply(ba_files, process_ba)
ba_N <- tibble(trait="hip_bone_area", N=29000)  # From publication
save("ba_ss", "ba_N", file="../data/processed/cleaned_sum_stats/ba.RData")

# NEALE LAB MEGA-GWAS
process_neale <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/neale_mega/", ss_file)) %>%
    separate(variant, c("CHR", "POS", "A2", "A1")) %>%
    mutate(z=beta / se,
           loc=paste(CHR, POS, sep=":")) %>%
    select(loc, A1, A2, z) %>%
    inner_join(keep_snps, by="loc") %>%
    select(rsid, A1, A2, z)
}
neale_files <- c(
  trunk_fat_pct="23127_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  leg_fat_pct="23115_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  bmr="23105_irnt.gwas.imputed_v3.both_sexes.tsv.gz"
)
neale_ss <- lapply(neale_files, process_neale)
neale_N <- tibble(trait=names(neale_files),
                  N=c(354600, 354600, 354800))
save("neale_ss", "neale_N", 
     file="../data/processed/cleaned_sum_stats/neale.RData")

# DOHERTY 2018 -- DEVICE-MEASURED PHYSICAL ACTIVITY
process_doherty <- function(ss_file) {
  read_csv(paste0("../data/raw/sum_stats/doherty2018/", ss_file)) %>%
    mutate(z=BETA / SE,
           loc=paste(CHR, BP, sep="_")) %>%
    select(rsid=SNP, loc, A1=ALLELE1, A2=ALLELE0, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
doherty_files <- list(
  act_moderate="Doherty-2018-NatureComms-moderate.csv.gz",
  act_sedentary="Doherty-2018-NatureComms-sedentary.csv.gz",
  act_walking="Doherty-2018-NatureComms-walking.csv.gz",
  act_overall="Doherty-2018-NatureComms-overall-activity.csv.gz"
)
doherty_ss <- lapply(doherty_files, process_doherty)
doherty_N <- tibble(trait=names(doherty_files), 
                    N=rep(91100, length(doherty_files)))  # From publication
save("doherty_ss", "doherty_N", 
     file="../data/processed/cleaned_sum_stats/doherty.RData")

# MERINO 2019 -- MACRONUTRIENTS

process_macros <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/macros/", ss_file)) %>%
    mutate(z=Effect / StdErr,
           A1=toupper(Allele1)) %>%
    select(rsid=MarkerName, A1, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
macros_files <- list(
  carbohydrate="UKB_CHARGE_Carb_MA1_overlapSNP.txt",
  fat="UKB_CHARGE_Fat_MA1_overlapSNP.txt",
  protein="UKB_CHARGE_Protein_MA1_overlapSNP.txt",
  carbohydrate_adjBMI="UKB_CHARGE_CarbBMI_MA1_overlapSNP.txt",
  fat_adjBMI="UKB_CHARGE_FatBMI_MA1_overlapSNP.txt",
  protein_adjBMI="UKB_CHARGE_ProteinBMI_MA1_overlapSNP.txt"
)
macros_ss <- lapply(macros_files, process_macros)
macros_N <- tibble(trait=names(macros_files), 
                   N=rep(91100, length(macros_files)))  # From publication
save("macros_ss", "macros_N", 
     file="../data/processed/cleaned_sum_stats/macros.RData")

# NEALE LAB MEGA-GWAS
process_ukb_dietary <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/neale_mega/", ss_file)) %>% 
    separate(variant, c("CHR", "POS", "A2", "A1")) %>%
    mutate(z=beta / se,
           loc=paste(CHR, POS, sep=":")) %>%
    select(loc, A1, A2, z) %>%
    inner_join(keep_snps, by="loc") %>%
    select(rsid, A1, A2, z)
}
ukb_dietary_files <- c(
  # neale_food_weight="100001_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  neale_energy="100002_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  # neale_protein="100003_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  # neale_fat="100004_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  # neale_carbohydrate="100005_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  neale_sfa="100006_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  neale_pufa="100007_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  neale_sugar="100008_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  neale_fiber="100009_irnt.gwas.imputed_v3.both_sexes.tsv.gz",
  neale_portion_size="100010_gwas.imputed_v3.both_sexes.tsv.gz"
)
ukb_dietary_ss <- lapply(ukb_dietary_files, process_ukb_dietary)
ukb_dietary_N <- tibble(trait=names(ukb_dietary_files),
                        N=rep(51500, length(ukb_dietary_files)))  # From s.s.
save("ukb_dietary_ss", "ukb_dietary_N", 
     file="../data/processed/cleaned_sum_stats/ukb_dietary.RData")

# ACCELEROMETER-BASED SLEEP PHENOTYPES (JONES 2019)
process_sleep <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/jones2019/", ss_file)) %>%
    separate(VAR_ID, c("CHR", "POS", "A2", "A1"), "_") %>%
    mutate(z=BETA / SE,
           loc=paste(CHR, POS, sep=":")) %>%
    select(loc, A1, A2, z) %>%
    inner_join(keep_snps, by="loc") %>%
    select(rsid, A1, A2, z)
}
sleep_files <- list(
  # sleep_diurn_inact="GWAS_Accel_eu_dv1.MeanSleepDiurnalInactRN.1.txt",
  sleep_duration="GWAS_Accel_eu_dv1.MeanSleepDurationRN.1.txt",
  sleep_eff="GWAS_Accel_eu_dv1.MeanSleepEfficiencyRN.1.txt",
  # sleep_L5time="GWAS_Accel_eu_dv1.MeanSleepL5timeRN.1.txt",
  # sleep_M10time="GWAS_Accel_eu_dv1.MeanSleepM10timeRN.1.txt",
  sleep_midpoint="GWAS_Accel_eu_dv1.MeanSleepMidPointRN.1.txt"
  # sleep_num_eps="GWAS_Accel_eu_dv1.MeanSleepNumEpisodesRN.1.txt",
  # sleep_duration_sdev="GWAS_Accel_eu_dv1.StDevSleepDurationRN.1.txt"
)
sleep_ss <- lapply(sleep_files, process_sleep)
sleep_N <- tibble(trait=names(sleep_files), 
                  N=rep(85700, length(sleep_files)))  # From publication
save("sleep_ss", "sleep_N", 
     file="../data/processed/cleaned_sum_stats/sleep.RData")

# SSGAC GWAS: RISK BEHAVIORS (KARLSON LINNER 2019), 
# EDUCATIONAL & COGNITIVE (LEE 2018), NEUROTICISM & DEPRESSION (TURLEY 2018) 
process_ssgac <- function(ss_file) {
  read_tsv(paste0("../data/raw/sum_stats/ssgac/", ss_file)) %>%
    mutate(z=Beta / SE) %>%
    select(rsid=MarkerName, A1, A2, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
ssgac_files <- list(
  # ssgac="RISK_GWAS_MA_UKB+replication.txt",
  alcohol_drinks_per_wk="DRINKS_PER_WEEK_GWAS.txt",
  ever_smoke="EVER_SMOKER_GWAS_MA_UKB+TAG.txt",
  risk_PC1="RISK_PC1_GWAS.txt",
  edu_attainment="GWAS_EA_excl23andMe.txt",
  cognitive_performance="GWAS_CP_all.txt",
  neuroticism="GWAS_NEUR_full.txt",
  depressive_symptoms="Depression_DS_Full.txt.gz"
)
ssgac_ss <- lapply(ssgac_files, process_ssgac)
ssgac_N <- tibble(  # From publication
  trait=names(ssgac_files),
  N=c(
    # 975400, 
    414300, 518600, 315900, 1131900, 257900, 168100,
    354900  # Note: this is an "effective sample size" (MTAG analysis)
    )
)
save("ssgac_ss", "ssgac_N", 
     file="../data/processed/cleaned_sum_stats/ssgac.RData")

# ANXIETY FACTOR SCORE GWAS (OTOWA 2016)
process_anxiety <- function(ss_file) {
  read_table2(paste0("../data/raw/sum_stats/otowa2016/", ss_file)) %>%
    mutate(z=Effect / StdErr,
           A1=toupper(Allele1),
           A2=toupper(Allele2)) %>%
    select(rsid=SNPID, A1, A2, z) %>%
    filter(rsid %in% keep_snps$rsid)
}
anxiety_files <- list(anxiety="anxiety.meta.full.fs.tbl.gz")
anxiety_ss <- lapply(anxiety_files, process_anxiety)
anxiety_N <- tibble(trait="anxiety", N=18200)  # From publication
save("anxiety_ss", "anxiety_N", 
     file="../data/processed/cleaned_sum_stats/anxiety.RData")
