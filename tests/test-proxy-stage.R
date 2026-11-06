# Proxy stage: find_proxy_candidates(), qualify_proxy_candidates(),
# select_proxies(). LDlink is replaced by a stand-in that writes the same file
# format as LDlinkR::LDproxy_batch(append = TRUE), so no network is needed.

px_dir <- tempfile("bnmf_proxy_")
dir.create(px_dir)

# ---- primary GWAS ------------------------------------------------------------
px_primary <- data.frame(
  VAR_ID = c("7_1000_A_G", "7_5000_C_G", "8_1000_A_C", "7_9000_A_G",   # sentinels
             "7_1100_C_T", "7_1200_A_G", "7_1300_C_T", "7_1400_A_T",
             "7_1800_C_T", "7_1900_T_C", "7_5100_A_G"),
  EA = c("G", "G", "C", "G", "T", "G", "T", "T", "T", "C", "G"),
  BETA = c(0.3, 0.2, 0.2, 0.2, 0.25, 0.01, -0.2, 0.2, 0.2, -0.2, 0.2),
  SE = 0.02,
  P_VALUE = c(1e-9, 1e-12, 1e-8, 1e-10, 1e-6, 0.5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-7),
  stringsAsFactors = FALSE
)
px_primary_path <- write_tsv_fixture(px_primary, px_dir, "primary.txt")

# ---- sentinels ---------------------------------------------------------------
px_sentinels <- data.frame(
  VAR_ID = c("7_1000_A_G", "7_5000_C_G", "8_1000_A_C", "7_9000_A_G"),
  ChrPos = c("7:1000", "7:5000", "8:1000", "7:9000"),
  Risk_Allele = c("G", "G", "C", "G"),
  Nonrisk_Allele = c("A", "C", "A", "A"),
  PVALUE_main = c(1e-9, 1e-12, 1e-8, 1e-10),
  stringsAsFactors = FALSE
)
# 7_9000_A_G is a healthy sentinel; the other three were flagged.
px_need <- data.frame(
  VAR_ID = c("7_1000_A_G", "7_5000_C_G", "8_1000_A_C"),
  reason = c("high_missingness", "ambiguous", "high_missingness"),
  PVALUE_main = c(1e-9, 1e-12, 1e-8),
  stringsAsFactors = FALSE
)

# ---- stand-in for LDlink -----------------------------------------------------
ld_row <- function(rs, coord, alleles, r2, corr) {
  data.frame(RS_Number = rs, Coord = coord, Alleles = alleles, MAF = 0.3,
             Distance = 100, Dprime = 1, R2 = r2, Correlated_Alleles = corr,
             FORGEdb = 5, RegulomeDB = "5", Function = NA, stringsAsFactors = FALSE)
}
ld_answers <- list(
  "chr7:1000" = rbind(
    ld_row("rs1",  "chr7:1000", "(A/G)", 1.00, "A=A,G=G"),   # the query itself
    ld_row("rs11", "chr7:1100", "(C/T)", 0.95, "A=C,G=T"),   # good proxy
    ld_row("rs12", "chr7:1200", "(A/G)", 0.99, "A=A,G=G"),   # primary P too high
    ld_row("rs13", "chr7:1300", "(C/T)", 0.90, "A=C,G=T"),   # primary risk allele opposes the sentinel
    ld_row("rs14", "chr7:1400", "(A/T)", 0.97, "A=A,G=T"),   # palindromic
    ld_row("rs15", "chr7:1500", "(C/T)", 0.85, "A=C,G=T"),   # not in the primary GWAS
    ld_row("rs16", "chr7:1600", "(G/T)", 0.70, "A=G,G=T"),   # r2 too low
    ld_row("rs17", "chr7:9000", "(A/G)", 0.90, "A=A,G=G"),   # another sentinel
    ld_row("rs18", "chr7:1800", "(C/T)", 0.96, "A=C,G=T"),   # too little trait data
    ld_row("rs19", "chr7:1900", "(C/T)", 0.82, "A=C,G=T"),   # good, but weaker LD than rs11
    ld_row(".",    "chr7:2000", "(-/A)", 0.90, "A=-,G=A")    # indel
  ),
  "chr7:5000" = rbind(
    ld_row("rs5",  "chr7:5000", "(C/G)", 1.00, "C=C,G=G"),
    ld_row("rs21", "chr7:5100", "(A/G)", 0.90, "C=A,G=G")
  )
  # chr8:1000 is deliberately absent: LDlink "fails" for it.
)
ld_calls <- list()
fake_ldproxy <- function(snps, population, token, out_file, genome_build = "grch37") {
  ld_calls[[length(ld_calls) + 1L]] <<- snps
  for (snp in snps) {
    if (is.null(ld_answers[[snp]])) next
    df <- cbind(query_snp = snp, ld_answers[[snp]])
    suppressWarnings(write.table(df, file = out_file, append = TRUE, quote = FALSE,
                                 row.names = TRUE, col.names = !file.exists(out_file),
                                 sep = "\t"))
  }
  invisible(NULL)
}

ld_dir <- file.path(px_dir, "ldlink")
candidates <- quietly(find_proxy_candidates(
  px_need, token = "unused", population = "EUR", output_dir = ld_dir,
  r2_min = 0.8, query_fun = fake_ldproxy
))

# ---- find_proxy_candidates ---------------------------------------------------

test_that("candidates are the LD partners at or above the r2 threshold, as SNPs, without the query itself", {
  expect_setequal(
    candidates$proxy_ChrPos[candidates$VAR_ID == "7_1000_A_G"],
    c("7:1100", "7:1200", "7:1300", "7:1400", "7:1500", "7:9000", "7:1800", "7:1900")
  )
  expect_identical(candidates$proxy_ChrPos[candidates$VAR_ID == "7_5000_C_G"], "7:5100")
  expect_false("7:1600" %in% candidates$proxy_ChrPos)   # r2 0.70
  expect_false("7:2000" %in% candidates$proxy_ChrPos)   # indel
  expect_false("7:1000" %in% candidates$proxy_ChrPos)   # the query itself
})

test_that("a sentinel LDlink did not answer for is reported, after one retry", {
  expect_identical(attr(candidates, "unanswered"), "8_1000_A_C")
  # First call asks for all three; the retry asks only for the unanswered one.
  expect_setequal(ld_calls[[1]], c("chr7:1000", "chr7:5000", "chr8:1000"))
  expect_identical(ld_calls[[2]], "chr8:1000")
})

test_that("a rerun reuses the answers on disk and only asks again for unanswered sentinels", {
  calls_before <- length(ld_calls)
  again <- quietly(find_proxy_candidates(
    px_need, token = "unused", population = "EUR", output_dir = ld_dir,
    r2_min = 0.8, query_fun = fake_ldproxy
  ))
  new_calls <- unlist(ld_calls[seq_along(ld_calls) > calls_before])
  expect_true(length(new_calls) > 0L)
  expect_true(all(new_calls == "chr8:1000"))
  expect_equal(nrow(again), nrow(candidates))
})

test_that("LDlink answers made with different settings are not reused", {
  expect_error(
    quietly(find_proxy_candidates(
      px_need, token = "unused", population = "EAS", output_dir = ld_dir,
      r2_min = 0.8, query_fun = fake_ldproxy
    )),
    "different settings"
  )
})

# ---- qualify_proxy_candidates ------------------------------------------------

qualified <- quietly(qualify_proxy_candidates(
  candidates, main_ss_filepath = px_primary_path, sentinels = px_sentinels,
  main_pval_cutoff = 0.05
))
status_of <- function(chr_pos) qualified$status[qualified$proxy_ChrPos == chr_pos]

test_that("each candidate is judged against the primary GWAS with one outcome", {
  expect_identical(status_of("7:1100"), "qualified")
  expect_identical(status_of("7:1200"), "main_pvalue_above_cutoff")
  expect_identical(status_of("7:1300"), "risk_allele_opposes_sentinel")
  expect_identical(status_of("7:1400"), "strand_ambiguous")
  expect_identical(status_of("7:1500"), "not_in_main_gwas")
  expect_identical(status_of("7:9000"), "is_a_sentinel")
  expect_identical(status_of("7:5100"), "qualified")
})

test_that("a candidate is matched to the primary GWAS whatever allele order the primary uses", {
  # LDlink lists chr7:1900 as (C/T); the primary lists it as T_C.
  row <- qualified[qualified$proxy_ChrPos == "7:1900", ]
  expect_identical(row$status, "qualified")
  expect_identical(row$proxy_VAR_ID, "7_1900_T_C")
  expect_identical(row$Risk_Allele, "T")
  expect_identical(row$Nonrisk_Allele, "C")
})

test_that("qualified candidates carry what the trait fetch needs", {
  ok <- qualified[qualified$status == "qualified", ]
  expect_true(all(c("proxy_VAR_ID", "proxy_ChrPos", "REF", "ALT", "Risk_Allele",
                    "Nonrisk_Allele", "PVALUE_main", "BETA_main", "SE_main") %in% names(ok)))
  expect_false(anyNA(ok$Risk_Allele))
})

# ---- select_proxies ----------------------------------------------------------

px_traits <- c("t1", "t2", "t3", "t4")
z_candidates <- matrix(1, nrow = 4, ncol = 4,
                       dimnames = list(c("7:1100", "7:1800", "7:1900", "7:5100"), px_traits))
z_candidates["7:1800", 2:4] <- NA   # 25% coverage

selection <- quietly(select_proxies(px_need, qualified, z_candidates, min_nonmissing = 0.8))

test_that("the proxy chosen is the qualified candidate in strongest LD with enough trait data", {
  chosen <- selection$chosen
  # chr7:1800 has higher r2 (0.96) than chr7:1100 (0.95) but too little trait data.
  expect_identical(chosen$proxy_ChrPos[chosen$original_SNP == "7_1000_A_G"], "7:1100")
  expect_identical(chosen$proxy_ChrPos[chosen$original_SNP == "7_5000_C_G"], "7:5100")
  expect_equal(nrow(chosen), 2L)
})

test_that("every candidate ends with one recorded outcome", {
  audit <- selection$candidates
  outcome <- function(chr_pos) audit$status[audit$proxy_ChrPos == chr_pos]
  expect_identical(outcome("7:1100"), "chosen")
  expect_identical(outcome("7:1800"), "insufficient_trait_data")
  expect_identical(outcome("7:1900"), "qualified_not_chosen")
  expect_identical(outcome("7:1300"), "risk_allele_opposes_sentinel")
})

test_that("a flagged sentinel with no proxy is reported with the reason and is not rescued", {
  no_proxy <- selection$no_proxy
  expect_identical(no_proxy$VAR_ID, "8_1000_A_C")
  expect_identical(no_proxy$reason, "ldlink_no_answer")
  expect_false("8_1000_A_C" %in% selection$chosen$original_SNP)
})

test_that("one variant is never the proxy for two sentinels; the stronger sentinel chooses first", {
  need <- data.frame(VAR_ID = c("9_100_A_G", "9_900_A_G"), reason = "high_missingness",
                     PVALUE_main = c(1e-9, 1e-20), stringsAsFactors = FALSE)
  qual <- data.frame(
    VAR_ID = c("9_100_A_G", "9_100_A_G", "9_900_A_G"),
    proxy_ChrPos = c("9:500", "9:200", "9:500"),
    proxy_VAR_ID = c("9_500_C_T", "9_200_C_T", "9_500_C_T"),
    proxy_rsID = c("rsX", "rsY", "rsX"),
    r2 = c(0.99, 0.85, 0.95),
    PVALUE_main = 1e-6, status = "qualified", stringsAsFactors = FALSE
  )
  z <- matrix(1, 2, 2, dimnames = list(c("9:500", "9:200"), c("a", "b")))
  picked <- quietly(select_proxies(need, qual, z, min_nonmissing = 0.8))$chosen
  expect_identical(picked$proxy_ChrPos[picked$original_SNP == "9_900_A_G"], "9:500")
  expect_identical(picked$proxy_ChrPos[picked$original_SNP == "9_100_A_G"], "9:200")

  # With no second candidate, the weaker sentinel gets no proxy and the reason says why.
  only_shared <- quietly(select_proxies(need, qual[qual$proxy_ChrPos == "9:500", ], z,
                                        min_nonmissing = 0.8))
  expect_identical(only_shared$no_proxy$VAR_ID, "9_100_A_G")
  expect_identical(only_shared$no_proxy$reason, "candidate_used_by_another_sentinel")
})

# ---- legacy choose_proxies() -------------------------------------------------

test_that("choose_proxies no longer rescues a flagged sentinel that has no proxy", {
  legacy_dir <- file.path(px_dir, "legacy")
  dir.create(legacy_dir)
  out_file <- file.path(legacy_dir, "combined_query_snp_list_grch37.txt")
  fake_ldproxy("chr7:1000", "EUR", "unused", out_file)   # answers on disk ...
  only_self <- read.table(out_file, sep = "\t", header = TRUE, row.names = NULL)
  only_self <- only_self[only_self$Coord == "chr7:1000", -1]
  unlink(out_file)                                       # ... reduced to the query itself
  write.table(only_self, out_file, quote = FALSE, row.names = TRUE, sep = "\t")
  saveRDS(list(query_snps = "chr7:1000", population = "EUR", genome_build = "grch37"),
          file.path(legacy_dir, "combined_query_snp_list_grch37.meta.rds"))

  # 75% trait coverage: the old rescue rule (clean alleles, >= 50% data) kept it.
  zmat <- matrix(c(1, 1, 1, NA), nrow = 1, dimnames = list("7:1000", px_traits))
  res <- quietly(choose_proxies(
    need_proxies = data.frame(VAR_ID = "7_1000_A_G", rsID = "rs1", stringsAsFactors = FALSE),
    rsid_map_dir = legacy_dir,
    pruned_variants = data.frame(VAR_ID = c("7_1000_A_G", "7_9000_A_G"), stringsAsFactors = FALSE),
    zmat_fullset = zmat, token = "unused", population = "EUR", output_dir = legacy_dir
  ))
  expect_false("7_1000_A_G" %in% res[[1]])
  expect_identical(res[[1]], "7_9000_A_G")
})
