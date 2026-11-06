# Shared setup for the pipeline tests: source the helpers once and provide
# small synthetic summary-statistics files.

repo_dir <- Sys.getenv("BNMF_REPO_DIR", unset = normalizePath(".."))

suppressWarnings(suppressPackageStartupMessages({
  source(file.path(repo_dir, "scripts", "choose_variants_2025.R"))
  source(file.path(repo_dir, "scripts", "prep_bNMF_2025.R"))
}))

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(
    out <- suppressMessages(suppressWarnings(force(expr)))
  )
  out
}

write_tsv_fixture <- function(df, dir, name) {
  path <- file.path(dir, name)
  data.table::fwrite(df, path, sep = "\t", na = "NA", quote = FALSE)
  path
}

# -----------------------------------------------------------------------------
# Discovery / primary GWAS fixtures
# -----------------------------------------------------------------------------
# primary : VAR_ID layout in reference REF/ALT order, effect allele column = ALT
# discB   : SNP / EA / NEA layout, allele order deliberately flipped relative to
#           the primary so no VAR_ID string could match
# discC   : VAR_ID layout with NO effect-allele column (effect = ALT is assumed)
# discD   : VAR_ID layout whose effect-allele column is REF, not ALT
make_discovery_fixtures <- function(dir = tempfile("bnmf_fixture_")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  primary <- data.frame(
    VAR_ID = c("1_1000_A_G", "1_2000_C_T", "1_3000_A_C", "1_4000_G_T",
               "1_5000_A_G", "1_6000_C_T", "2_1000_A_G", "2_3000_A_G",
               "2_4000_A_G", "3_1000_A_G", "3_2000_A_G"),
    Effect_Allele_PH = c("G", "T", "C", "T", "G", "T", "G", "G", "G", "G", "G"),
    BETA = c(0.2, -0.1, 0.05, 0.1, NA, 0, 0.3, 0.1, 0.1, 0.1, -0.1),
    SE = 0.01,
    P_VALUE = c(1e-20, 1e-3, 0.20, 1e-4, 1e-5, 1e-5, 1e-30, NA, 1e-3, 1e-3, 1e-3),
    stringsAsFactors = FALSE
  )

  discB <- data.frame(
    SNP = c("1:1000", "1:2000", "1:3000", "1:4000", "1:5000", "1:6000",
            "2:2000", "2:3000", "2:4000"),
    EA  = c("A", "C", "C", "G", "G", "T", "C", "G", "G"),
    NEA = c("G", "T", "A", "T", "A", "C", "A", "A", "A"),
    BETA = c(-0.3, 0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, NA),
    SE = 0.02,
    P_VALUE = c(1e-9, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-9),
    stringsAsFactors = FALSE
  )

  discC <- data.frame(
    VAR_ID = "3_1000_A_G", BETA = 0.2, SE = 0.02, P_VALUE = 1e-9,
    stringsAsFactors = FALSE
  )

  discD <- data.frame(
    VAR_ID = "3_2000_A_G", Effect_Allele_PH = "A", BETA = 0.2, SE = 0.02,
    P_VALUE = 1e-9, stringsAsFactors = FALSE
  )

  paths <- c(
    primary = write_tsv_fixture(primary, dir, "primary.txt"),
    discB   = write_tsv_fixture(discB,   dir, "discB.txt"),
    discC   = write_tsv_fixture(discC,   dir, "discC.txt"),
    discD   = write_tsv_fixture(discD,   dir, "discD.txt")
  )

  gwas <- data.frame(
    ID = c("primary_EUR", "discB_EAS", "discC_AFR", "discD_SAS"),
    population = c("EUR", "EAS", "AFR", "SAS"),
    full_path = unname(paths),
    largest = c("Yes", "No", "No", "No"),
    stringsAsFactors = FALSE
  )

  list(dir = dir, paths = paths, gwas = gwas)
}

# -----------------------------------------------------------------------------
# Trait fixtures
# -----------------------------------------------------------------------------
# trait1 : VAR_ID layout with an effect-allele column
# trait2 : SNP / EA / NEA layout, reported on the opposite strand for 1:1000
# trait3 : VAR_ID layout with no effect-allele column
make_trait_fixtures <- function(dir = tempfile("bnmf_traits_")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  trait1 <- data.frame(
    VAR_ID = c("1_1000_A_G", "1_2000_C_T", "2_1000_A_G"),
    Effect_Allele_PH = c("G", "T", "G"),
    BETA = c(0.5, 0.2, -0.3), SE = 0.1, P_VALUE = c(1e-6, 0.05, 0.003),
    stringsAsFactors = FALSE
  )
  trait2 <- data.frame(
    SNP = c("1:1000", "1:2000", "2:6000"),
    EA  = c("C", "A", "T"),
    NEA = c("T", "C", "A"),
    BETA = c(0.4, 0.2, 0.1), SE = 0.1, P_VALUE = c(1e-4, 0.05, 0.3),
    stringsAsFactors = FALSE
  )
  trait3 <- data.frame(
    VAR_ID = "1_1000_A_G", BETA = 0.3, SE = 0.1, P_VALUE = 0.003,
    stringsAsFactors = FALSE
  )

  files <- c(
    trait1 = write_tsv_fixture(trait1, dir, "trait1.txt"),
    trait2 = write_tsv_fixture(trait2, dir, "trait2.txt"),
    trait3 = write_tsv_fixture(trait3, dir, "trait3.txt")
  )

  df_input <- data.frame(
    SNP = c("1:1000", "1:2000", "2:1000", "2:5000", "2:6000"),
    VAR_ID = c("1_1000_A_G", "1_2000_C_T", "2_1000_A_G", "2_5000_A_C", "2_6000_A_T"),
    REF = c("A", "C", "A", "A", "A"),
    ALT = c("G", "T", "G", "C", "T"),
    Risk_Allele = c("G", "C", "G", "C", "T"),
    Nonrisk_Allele = c("A", "T", "A", "A", "A"),
    PVALUE_main = c(1e-20, 1e-3, 1e-30, 1e-4, 1e-4),
    BETA_main = c(0.2, -0.1, 0.3, 0.1, 0.1),
    SE_main = 0.01,
    stringsAsFactors = FALSE
  )

  list(
    dir = dir, files = files, df_input = df_input,
    sizes = c(trait1 = 10000, trait2 = 20000, trait3 = 30000)
  )
}
