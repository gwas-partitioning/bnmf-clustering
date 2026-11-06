# Variant selection: get_sig_snps(), get_biggest_gwas(), clumping/pruning
# ordering, and the proxy-need classifier.

fx <- make_discovery_fixtures()
vars_sig <- quietly(get_sig_snps(gwas = fx$gwas, PVCUTOFF = 1e-6))
vars_main <- quietly(get_biggest_gwas(fx$paths[["primary"]], vars_sig,
                                      main_pval_cutoff = 0.05))
dropped <- attr(vars_main, "dropped")

reason_for <- function(chr_pos) {
  dropped$reason[paste(dropped$CHR, dropped$POS, sep = ":") == chr_pos]
}

# ---- canonical key -----------------------------------------------------------

test_that("canonical key ignores which allele is listed first", {
  expect_identical(canonical_variant_key(1, 1000, "A", "G"),
                   canonical_variant_key("1", "1000", "G", "A"))
  expect_false(identical(canonical_variant_key(1, 1000, "A", "G"),
                         canonical_variant_key(1, 1000, "A", "C")))
})

# ---- get_sig_snps ------------------------------------------------------------

test_that("get_sig_snps names the discovery P explicitly and never returns a bare PVALUE", {
  expect_true("PVALUE_discovery" %in% names(vars_sig))
  expect_false("PVALUE" %in% names(vars_sig))
})

test_that("files with opposite allele order land on the same variant key", {
  key_primary <- vars_sig$VAR_KEY[vars_sig$GWAS == "primary_EUR" & vars_sig$POS == 1000 & vars_sig$CHR == 1]
  key_discB   <- vars_sig$VAR_KEY[vars_sig$GWAS == "discB_EAS"   & vars_sig$POS == 1000 & vars_sig$CHR == 1]
  expect_length(key_primary, 1L)
  expect_identical(key_primary, key_discB)
})

test_that("risk allele follows the stated effect allele, not the allele order", {
  # discB lists 1:1000 as EA = A with a negative effect, so the risk allele is G.
  row <- vars_sig[vars_sig$GWAS == "discB_EAS" & vars_sig$CHR == 1 & vars_sig$POS == 1000, ]
  expect_identical(row$Risk_Allele, "G")
  expect_identical(row$Nonrisk_Allele, "A")
})

test_that("an explicit effect-allele column wins over the ALT position in VAR_ID", {
  # discD: VAR_ID 3_2000_A_G, effect allele A, positive effect -> risk allele A.
  row <- vars_sig[vars_sig$GWAS == "discD_SAS", ]
  expect_identical(row$Risk_Allele, "A")
})

test_that("with no effect-allele column the effect is assumed to be ALT, and this is reported", {
  row <- vars_sig[vars_sig$GWAS == "discC_AFR", ]
  expect_identical(row$Risk_Allele, "G")
  qc <- attr(vars_sig, "input_qc")
  expect_match(qc$effect_allele_rule[qc$GWAS == "discC_AFR"], "ALT", fixed = TRUE)
  expect_match(qc$effect_allele_rule[qc$GWAS == "discD_SAS"], "Effect_Allele_PH", fixed = TRUE)
})

test_that("a manifest effect_allele column overrides the detected rule", {
  gwas_override <- fx$gwas
  gwas_override$effect_allele <- c(NA, NA, "REF", NA)
  overridden <- quietly(get_sig_snps(gwas = gwas_override, PVCUTOFF = 1e-6))
  # discC: 3_1000_A_G, positive effect, effect allele forced to REF -> risk A.
  expect_identical(overridden$Risk_Allele[overridden$GWAS == "discC_AFR"], "A")
})

# ---- get_biggest_gwas --------------------------------------------------------

test_that("get_biggest_gwas returns two named P columns and no bare PVALUE", {
  expect_true(all(c("PVALUE_discovery", "PVALUE_main") %in% names(vars_main)))
  expect_false(any(c("PVALUE", "PVALUE.Pop") %in% names(vars_main)))
})

test_that("a variant stays in when discovery P qualifies it and main P is only nominal", {
  row <- vars_main[vars_main$VAR_ID == "1_2000_C_T", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$PVALUE_discovery, 1e-8)
  expect_equal(row$PVALUE_main, 1e-3)
  expect_identical(row$Risk_Allele, "C")
})

test_that("kept variants carry the primary GWAS allele order in VAR_ID", {
  # discB wrote 1:1000 as A (effect) / G; the primary's order is A_G.
  expect_true("1_1000_A_G" %in% vars_main$VAR_ID)
  expect_setequal(vars_main$VAR_ID,
                  c("1_1000_A_G", "1_2000_C_T", "2_1000_A_G", "3_1000_A_G", "3_2000_A_G"))
  expect_identical(vars_main$ChrPos[vars_main$VAR_ID == "1_1000_A_G"], "1:1000")
  expect_identical(vars_main$ChrPos_LDlink[vars_main$VAR_ID == "1_1000_A_G"], "chr1:1000")
})

test_that("the strongest discovery record is the one retained per variant", {
  row <- vars_main[vars_main$VAR_ID == "1_1000_A_G", ]
  expect_equal(row$PVALUE_discovery, 1e-20)
  expect_identical(row$GWAS, "primary_EUR")
})

test_that("every removed variant is reported with one reason", {
  expect_identical(reason_for("2:2000"), "not_in_main_gwas")
  expect_identical(reason_for("2:4000"), "no_discovery_direction")
  expect_identical(reason_for("1:5000"), "no_main_direction")   # main BETA is NA
  expect_identical(reason_for("1:6000"), "no_main_direction")   # main BETA is 0
  expect_identical(reason_for("2:3000"), "main_pvalue_missing")
  expect_identical(reason_for("1:3000"), "main_pvalue_above_cutoff")
  expect_identical(reason_for("1:4000"), "risk_allele_mismatch")
  expect_equal(nrow(dropped) + nrow(vars_main), length(unique(vars_sig$VAR_KEY)))
})

test_that("honouring the effect-allele column avoids a false risk-allele mismatch", {
  # discD states effect allele A for 3_2000_A_G with a positive effect, and the
  # primary has effect allele G with a negative effect: both say A is the risk
  # allele. Reading discD as "effect = ALT" would have called this a mismatch.
  expect_true("3_2000_A_G" %in% vars_main$VAR_ID)
  expect_identical(vars_main$Risk_Allele[vars_main$VAR_ID == "3_2000_A_G"], "A")
})

test_that("a risk-allele mismatch is removed however significant the main GWAS is", {
  # 1:4000 has main P = 1e-4; the old rule kept mismatches below a Bonferroni band.
  strict <- quietly(get_biggest_gwas(fx$paths[["primary"]], vars_sig, main_pval_cutoff = 1))
  expect_false("1_4000_G_T" %in% strict$VAR_ID)
})

test_that("drop counts are summarised per discovery GWAS", {
  summary_tbl <- attr(vars_main, "drop_summary")
  expect_true(all(c("GWAS", "reason", "n") %in% names(summary_tbl)))
  n_mismatch_B <- summary_tbl$n[summary_tbl$GWAS == "discB_EAS" &
                                  summary_tbl$reason == "risk_allele_mismatch"]
  expect_equal(n_mismatch_B, 1L)
})

test_that("per-file agreement with the primary GWAS is reported", {
  concordance <- attr(vars_main, "concordance")
  expect_true(all(c("GWAS", "n_testable", "n_concordant", "n_flipped_order") %in% names(concordance)))
  rowB <- concordance[concordance$GWAS == "discB_EAS", ]
  # Testable = in the primary, both directions known, main P below the cutoff:
  # 1:1000 (concordant), 1:2000 (concordant), 1:4000 (discordant).
  expect_equal(rowB$n_testable, 3L)
  expect_equal(rowB$n_concordant, 2L)
})

test_that("BETA_main is the effect of ALT even when the primary's effect allele is REF", {
  dir <- tempfile("bnmf_ref_effect_")
  dir.create(dir)
  # Primary lists REF (A) as the effect allele with a negative effect, so ALT
  # (G) has a positive effect and is the risk allele.
  primary <- data.frame(VAR_ID = "5_1000_A_G", Effect_Allele_PH = "A",
                        BETA = -0.4, SE = 0.01, P_VALUE = 1e-12,
                        stringsAsFactors = FALSE)
  path <- write_tsv_fixture(primary, dir, "primary_ref_effect.txt")
  gwas <- data.frame(ID = "p", population = "EUR", full_path = path,
                     stringsAsFactors = FALSE)
  sig <- quietly(get_sig_snps(gwas, PVCUTOFF = 1e-6))
  out <- quietly(get_biggest_gwas(path, sig))
  expect_identical(out$ALT, "G")
  expect_identical(out$Risk_Allele, "G")
  expect_equal(out$BETA_main, 0.4)
})

test_that("a gzipped or space-delimited primary GWAS gives the same result", {
  primary_df <- data.table::fread(fx$paths[["primary"]], data.table = FALSE)

  gz_path <- file.path(fx$dir, "primary_copy.txt.gz")
  data.table::fwrite(primary_df, gz_path, sep = "\t", na = "NA", quote = FALSE)
  from_gz <- quietly(get_biggest_gwas(gz_path, vars_sig))

  space_path <- file.path(fx$dir, "primary_space.txt")
  data.table::fwrite(primary_df, space_path, sep = " ", na = "NA", quote = FALSE)
  from_space <- quietly(get_biggest_gwas(space_path, vars_sig))

  strip <- function(x) { attributes(x) <- attributes(x)[c("names", "row.names", "class")]; x }
  expect_equal(strip(from_gz), strip(vars_main))
  expect_equal(strip(from_space), strip(vars_main))
})

test_that("a primary GWAS located by CHR and POS columns is read without a VAR_ID", {
  primary_df <- data.table::fread(fx$paths[["primary"]], data.table = FALSE)
  parts <- data.table::tstrsplit(primary_df$VAR_ID, "_", fixed = TRUE)
  no_varid <- data.frame(
    CHR = paste0("chr", parts[[1]]), POS = parts[[2]],
    NEA = parts[[3]], EA = primary_df$Effect_Allele_PH,
    BETA = primary_df$BETA, SE = primary_df$SE, P_VALUE = primary_df$P_VALUE,
    stringsAsFactors = FALSE
  )
  path <- write_tsv_fixture(no_varid, fx$dir, "primary_chr_pos.txt")
  out <- quietly(get_biggest_gwas(path, vars_sig))
  expect_setequal(out$VAR_ID, vars_main$VAR_ID)
  expect_equal(out$PVALUE_main[match(vars_main$VAR_ID, out$VAR_ID)], vars_main$PVALUE_main)
})

# ---- ordering in clumping and LD pruning -------------------------------------

clump_input <- data.frame(
  VAR_ID = c("1_1000_A_G", "1_1500_C_T"),
  CHR = 1L, POS = c(1000L, 1500L),
  PVALUE_main = c(1e-9, 1e-6),
  PVALUE_discovery = c(1e-7, 1e-12),
  stringsAsFactors = FALSE
)

test_that("snp_clump refuses to guess which P-value orders the variants", {
  expect_error(quietly(snp_clump(clump_input, id = "VAR_ID", window = 1e5)),
               "order_by")
})

test_that("snp_clump keeps the lead variant under the column it is told to use", {
  expect_identical(
    quietly(snp_clump(clump_input, id = "VAR_ID", window = 1e5, order_by = "PVALUE_main")),
    "1_1000_A_G"
  )
  expect_identical(
    quietly(snp_clump(clump_input, id = "VAR_ID", window = 1e5, order_by = "PVALUE_discovery")),
    "1_1500_C_T"
  )
})

test_that("ld_pruning_SNP.clip requires order_by before contacting LDlink", {
  expect_error(
    quietly(ld_pruning_SNP.clip(clump_input, pop = "EUR", token = "unused",
                                output_dir = tempfile("ld_"))),
    "order_by"
  )
})

# ---- proxy-need classifier ---------------------------------------------------

proxy_input <- data.frame(
  VAR_ID = c("1_1000_A_G", "1_2000_C_T"),
  REF = c("A", "C"), ALT = c("G", "T"),
  RS_Number = c("rs1", "rs2"),
  PVALUE_main = c(1e-9, 1e-6),
  stringsAsFactors = FALSE
)

test_that("a sentinel with no non-missingness value is an error, not a proxy reason", {
  expect_error(
    quietly(find_variants_needing_proxies(proxy_input, c("1_1000_A_G" = 1))),
    "non-missingness"
  )
})

test_that("a sentinel below the non-missingness threshold is flagged high_missingness", {
  flagged <- quietly(find_variants_needing_proxies(
    proxy_input, c("1_1000_A_G" = 1, "1_2000_C_T" = 0.5)
  ))
  expect_identical(flagged$VAR_ID, "1_2000_C_T")
  expect_identical(flagged$reason, "high_missingness")
})
