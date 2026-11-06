# Trait fetch: fetch_summary_stats() and the trait reader.

tf <- make_trait_fixtures()
fetched <- quietly(fetch_summary_stats(
  df_input = tf$df_input,
  trait_ss_files = tf$files,
  trait_ss_size = tf$sizes
))
z <- as.matrix(fetched$df_z)

test_that("fetch_summary_stats runs without a primary GWAS when the input already carries risk alleles", {
  expect_true(all(c("df_z", "df_N", "df_gwas", "trait_qc") %in% names(fetched)))
  expect_setequal(fetched$df_gwas$SNP, tf$df_input$SNP)
})

test_that("every input variant gets a z-matrix row, including one with no trait data", {
  expect_setequal(rownames(z), tf$df_input$SNP)
  expect_true(all(is.na(z["2:5000", ])))
  expect_identical(rownames(fetched$df_N), rownames(fetched$df_z))
})

test_that("z-scores are oriented to the risk allele", {
  expect_equal(z["1:1000", "trait1"], 5)    # effect allele G is the risk allele
  expect_equal(z["1:2000", "trait1"], -2)   # effect allele T, risk allele C
  expect_equal(z["2:1000", "trait1"], -3)
})

test_that("a trait reported on the opposite strand is matched by complement", {
  # trait2 lists 1:1000 as C/T, the complement of G/A; effect allele C == G.
  expect_equal(z["1:1000", "trait2"], 4)
  qc <- fetched$trait_qc
  expect_equal(qc$n_complemented[qc$trait == "trait2"], 1L)
})

test_that("alleles that match on neither strand are rejected and counted", {
  # trait2 lists 1:2000 as A/C; the variant is C/T.
  expect_true(is.na(z["1:2000", "trait2"]))
  qc <- fetched$trait_qc
  expect_equal(qc$n_allele_mismatch[qc$trait == "trait2"], 1L)
})

test_that("a palindromic variant is matched as reported", {
  expect_equal(z["2:6000", "trait2"], 1)
})

test_that("a trait file with no effect-allele column is read as effect = ALT, and this is reported", {
  expect_equal(z["1:1000", "trait3"], 3)
  qc <- fetched$trait_qc
  expect_match(qc$effect_allele_rule[qc$trait == "trait3"], "ALT", fixed = TRUE)
  expect_match(qc$effect_allele_rule[qc$trait == "trait1"], "Effect_Allele_PH", fixed = TRUE)
})

test_that("a trait effect-allele override changes the orientation", {
  overridden <- quietly(fetch_summary_stats(
    df_input = tf$df_input,
    trait_ss_files = tf$files,
    trait_ss_size = tf$sizes,
    trait_effect_allele = c(trait3 = "REF")
  ))
  expect_equal(as.matrix(overridden$df_z)["1:1000", "trait3"], -3)
})

test_that("unclean input is refused instead of being filtered silently", {
  bad_allele <- tf$df_input
  bad_allele$Risk_Allele[2] <- NA
  expect_error(
    quietly(fetch_summary_stats(bad_allele, trait_ss_files = tf$files, trait_ss_size = tf$sizes)),
    "Risk_Allele"
  )
  duplicated_position <- rbind(tf$df_input, tf$df_input[1, ])
  expect_error(
    quietly(fetch_summary_stats(duplicated_position, trait_ss_files = tf$files,
                                trait_ss_size = tf$sizes)),
    "duplicate"
  )
})

test_that("the legacy primary-GWAS path removes every risk-allele mismatch and says why", {
  legacy_input <- data.frame(
    SNP = c("1:1000", "1:2000", "2:1000"),
    Risk_Allele = c("G", "T", "G"),   # 1:2000 disagrees with the primary below
    stringsAsFactors = FALSE
  )
  primary <- data.frame(
    SNP = c("1:1000", "1:2000", "2:1000"),
    REF = c("A", "C", "A"), ALT = c("G", "T", "G"),
    BETA = c(0.2, -0.1, NA),
    P_VALUE = c(1e-20, 1e-30, 1e-10),
    stringsAsFactors = FALSE
  )
  legacy <- quietly(fetch_summary_stats(
    df_input = legacy_input, gwas_ss_file = primary,
    trait_ss_files = tf$files, trait_ss_size = tf$sizes, pval_cutoff = 0.05
  ))
  expect_identical(legacy$df_gwas$SNP, "1:1000")
  expect_identical(legacy$dropped$reason[legacy$dropped$SNP == "1:2000"], "risk_allele_mismatch")
  expect_identical(legacy$dropped$reason[legacy$dropped$SNP == "2:1000"], "no_main_direction")
})

# ---- parallel trait fetch -----------------------------------------------------

test_that("fetching traits with several workers gives the same result as one worker", {
  parallel_fetch <- quietly(fetch_summary_stats(
    df_input = tf$df_input,
    trait_ss_files = tf$files,
    trait_ss_size = tf$sizes,
    n_workers = 2
  ))
  expect_equal(parallel_fetch$df_z, fetched$df_z)
  expect_equal(parallel_fetch$df_N, fetched$df_N)
  expect_equal(parallel_fetch$trait_qc, fetched$trait_qc)
})

test_that("parallel workers write the same checkpoints a single worker can reuse", {
  cp_dir <- tempfile("bnmf_cp_")
  first <- quietly(fetch_summary_stats(
    df_input = tf$df_input, trait_ss_files = tf$files, trait_ss_size = tf$sizes,
    checkpoint_dir = cp_dir, n_workers = 2
  ))
  expect_setequal(list.files(cp_dir), paste0(names(tf$files), ".rds"))
  # A second, single-worker call must find every checkpoint valid.
  messages <- character(0)
  second <- withCallingHandlers(
    {
      utils::capture.output(res <- fetch_summary_stats(
        df_input = tf$df_input, trait_ss_files = tf$files, trait_ss_size = tf$sizes,
        checkpoint_dir = cp_dir, n_workers = 1
      ))
      res
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_equal(sum(grepl("[checkpoint]", messages, fixed = TRUE)), length(tf$files))
  expect_equal(second$df_z, first$df_z)
})

test_that("a trait that fails in a parallel worker is reported by name", {
  bad_dir <- tempfile("bnmf_bad_")
  dir.create(bad_dir)
  bad_path <- write_tsv_fixture(
    data.frame(VAR_ID = "1_1000_A_G", Effect_Allele_PH = "G", P_VALUE = 0.01),  # no BETA / SE
    bad_dir, "badtrait.txt"
  )
  files_with_bad <- c(tf$files, badtrait = bad_path)
  expect_error(
    quietly(fetch_summary_stats(
      df_input = tf$df_input, trait_ss_files = files_with_bad,
      trait_ss_size = c(tf$sizes, badtrait = 1000), n_workers = 2
    )),
    "badtrait"
  )
})
