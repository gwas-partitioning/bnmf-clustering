# Z-score alignment audit: scripts/audit_zscore_alignment.R
# No network access: Knowledge Portal responses are canned lists shaped like
# the parsed JSON of its BioIndex API.

suppressWarnings(suppressPackageStartupMessages({
  source(file.path(scripts_dir, "qc_report.R"))
  source(file.path(scripts_dir, "audit_zscore_alignment.R"))
}))

kp_row <- function(phenotype, beta, se = 0.01, p = 1e-6, n = 1e5) {
  list(phenotype = phenotype, beta = beta, stdErr = se, pValue = p, n = n, ancestry = "Mixed")
}

fake_portal_api <- function(responses) {
  # responses: named list, var_id -> list of kp_row(); unknown IDs return none.
  env <- new.env()
  env$calls <- character(0)
  env$get <- function(url) {
    env$calls <- c(env$calls, url)
    var_id <- sub("^.*q=(?:[A-Za-z]+,)?", "", url, perl = TRUE)
    rows <- responses[[var_id]]
    list(count = length(rows), data = if (length(rows)) rows else list())
  }
  env
}

test_that("kp_variant_url uses the ancestry-specific index when asked", {
  expect_match(kp_variant_url("1:100:A:G"), "phewas-associations?q=1:100:A:G", fixed = TRUE)
  expect_match(kp_variant_url("1:100:A:G", "EU"),
               "ancestry-phewas-associations?q=EU,1:100:A:G", fixed = TRUE)
})

test_that("kp_parse_rows reads beta, se, P and N, tolerating nulls", {
  d <- kp_parse_rows(list(kp_row("HDL", -0.02, p = 3e-9), list(phenotype = "BMI", beta = NULL)))
  expect_equal(d$phenotype, c("HDL", "BMI"))
  expect_equal(d$beta, c(-0.02, NA))
  expect_equal(d$se, c(0.01, NA))
  expect_equal(d$p[1], 3e-9)
})

test_that("kp_fetch_variant tries both allele orders and caches the result", {
  api <- fake_portal_api(list("1:100:G:A" = list(kp_row("HDL", 0.05))))
  cache <- file.path(tempfile(), "cache")
  d <- kp_fetch_variant("1", "100", "A", "G", cache_dir = cache, get_json = api$get)
  expect_equal(d$effect_allele, "A")
  expect_equal(d$other_allele, "G")
  expect_equal(d$kp_var_id, "1:100:G:A")
  expect_length(api$calls, 2L)
  again <- kp_fetch_variant("1", "100", "A", "G", cache_dir = cache, get_json = api$get)
  expect_identical(again, d)
  expect_length(api$calls, 2L)
  none <- kp_fetch_variant("1", "999", "A", "G", cache_dir = cache, get_json = api$get)
  expect_equal(nrow(none), 0L)
  kp_fetch_variant("1", "999", "A", "G", cache_dir = cache, get_json = api$get)
  expect_length(api$calls, 4L)
  expect_error(kp_fetch_variant("1", "5", "A", "G", get_json = function(url) list(detail = "x")),
               "Unexpected Knowledge Portal response")
})

audit_fixture <- function() {
  # Eight matrix variants (rs7 is A/T, which needs no special handling here).
  # The Portal reports alternate variants for the risk and the non-risk
  # allele; every beta implies +5 for the T2D risk allele.
  alignment <- data.frame(
    ChrPos = paste0("1:", 1:8 * 100),
    rsID = paste0("rs", 1:8),
    REF = c("A", "C", "A", "G", "A", "C", "A", "G"),
    ALT = c("G", "T", "C", "T", "G", "A", "T", "A"),
    Risk_Allele = c("G", "C", "A", "T", "A", "A", "T", "G"),
    stringsAsFactors = FALSE
  )
  nonrisk <- ifelse(alignment$Risk_Allele == alignment$REF, alignment$ALT, alignment$REF)
  effect <- ifelse(seq_len(8) %% 2 == 1, alignment$Risk_Allele, nonrisk)
  other <- ifelse(effect == alignment$REF, alignment$ALT, alignment$REF)
  pos <- sub("^1:", "", alignment$ChrPos)
  portal <- data.frame(
    phenotype = "TRAITX",
    beta = ifelse(effect == alignment$Risk_Allele, 0.05, -0.05),
    se = 0.01, p = 1e-7, n = 1e5,
    effect_allele = effect, other_allele = other,
    kp_var_id = paste("1", pos, other, effect, sep = ":"),
    ChrPos = alignment$ChrPos, stringsAsFactors = FALSE
  )
  z <- cbind(good = rep(5, 8), flipped = rep(-5, 8),
             scrambled = c(5, -5, 5, -5, 5, -5, 5, -5),
             unknown = rep(5, 8), blank = rep(5, 8))
  rownames(z) <- alignment$ChrPos
  list(alignment = alignment, portal = portal, z = z)
}

test_that("kp_fetch_matrix gathers every variant and records failures", {
  fx <- audit_fixture()
  p <- fx$portal
  responses <- setNames(lapply(seq_len(nrow(p)), function(i) {
    list(kp_row("TRAITX", p$beta[i]), kp_row("OTHER", 0.01))
  }), p$kp_var_id)
  responses[["1:800:G:A"]] <- NULL; responses[["1:800:A:G"]] <- NULL   # not in the Portal
  api <- fake_portal_api(responses)
  get <- function(url) if (grepl("1:300:", url, fixed = TRUE)) stop("HTTP 503") else api$get(url)
  a <- prepare_alignment(fx$alignment, fx$z)
  for (workers in c(1L, 2L)) {
    d <- suppressWarnings(quietly(kp_fetch_matrix(a, get_json = get, workers = workers)))
    expect_equal(attr(d, "failed"), "1:300")
    expect_equal(attr(d, "not_found"), "1:800")
    expect_equal(sort(unique(d$ChrPos)), paste0("1:", c(1, 2, 4, 5, 6, 7) * 100))
    expect_equal(nrow(d), 12L)
  }
  expect_warning(kp_fetch_matrix(a, get_json = get, workers = 1L), "could not be fetched")
})

test_that("suggest_kp_phenotypes matches codes and descriptions", {
  ph <- data.frame(
    name = c("HDL", "WHRadjBMI", "WaistC", "FI", "FIadjBMI", "BMI"),
    description = c("HDL cholesterol", "Waist-hip ratio adj BMI", "Waist circumference",
                    "Fasting insulin", "Fasting insulin adj BMI", "Body mass index (BMI)"),
    group = "", stringsAsFactors = FALSE
  )
  expect_equal(suggest_kp_phenotypes(ph, "HDL")$name[1], "HDL")
  expect_equal(suggest_kp_phenotypes(ph, "WHR_adjBMI")$name[1], "WHRadjBMI")
  expect_equal(suggest_kp_phenotypes(ph, "Waist circumference male")$name, "WaistC")
  expect_equal(suggest_kp_phenotypes(ph, "Fasting insulin")$name, c("FI", "FIadjBMI"))
  expect_equal(suggest_kp_phenotypes(ph, "Body mass index")$name[1], "BMI")
  expect_equal(nrow(suggest_kp_phenotypes(ph, "XINSG30")), 0L)
})

test_that("prepare_alignment orders rows like the matrix", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment[8:1, ], fx$z)
  expect_equal(a$ChrPos, rownames(fx$z))
  z_chr <- fx$z
  rownames(z_chr) <- paste0("chr", rownames(z_chr))
  expect_equal(prepare_alignment(fx$alignment, z_chr)$ChrPos, rownames(z_chr))
  no_rsid <- fx$alignment[setdiff(names(fx$alignment), "rsID")]
  expect_true(all(is.na(prepare_alignment(no_rsid, fx$z)$rsID)))
})

test_that("prepare_alignment rejects missing rows and bad risk alleles", {
  fx <- audit_fixture()
  expect_error(prepare_alignment(fx$alignment[-1, ], fx$z), "no alignment row")
  bad <- fx$alignment
  bad$Risk_Allele[1] <- "C"
  expect_error(prepare_alignment(bad, fx$z), "Risk_Allele")
})

test_that("compare_trait expresses Portal z-scores for the T2D risk allele", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment, fx$z)
  d <- compare_trait(fx$portal, a, fx$z[, "good"], "TRAITX")
  expect_equal(d$ChrPos, rownames(fx$z))
  expect_equal(d$portal_z, rep(5, 8))
  expect_true(all(d$status == "concordant"))
  z <- fx$z[, "good"]
  z[2] <- 1
  z[3] <- NA
  d <- compare_trait(fx$portal, a, z, "TRAITX")
  expect_equal(d$status[2:3], c("weak", "missing_in_matrix"))
  expect_equal(d$informative[1:3], c(TRUE, FALSE, FALSE))
  weak_portal <- fx$portal
  weak_portal$se[1] <- 1
  expect_equal(compare_trait(weak_portal, a, fx$z[, "good"], "TRAITX")$status[1], "weak")
  expect_null(compare_trait(fx$portal, a, fx$z[, "good"], "NOPE"))
})

test_that("trait_verdict separates aligned, flipped and scrambled traits", {
  fx <- audit_fixture()
  a <- prepare_alignment(fx$alignment, fx$z)
  verdict <- function(trait) trait_verdict(compare_trait(fx$portal, a, fx$z[, trait], "TRAITX"))
  good <- verdict("good")
  flipped <- verdict("flipped")
  scrambled <- verdict("scrambled")
  expect_equal(c(good$verdict, flipped$verdict, scrambled$verdict),
               c("PASS", "FLIPPED", "REVIEW"))
  expect_equal(good$n_informative, 8L)
  expect_equal(flipped$concordance, 0)
  expect_equal(scrambled$n_concordant, 4L)
  expect_lt(flipped$binom_p, 0.05)
  calls <- compare_trait(fx$portal, a, c(1:8) * 0.5 + 2, "TRAITX")
  calls$matrix_z <- seq(2, 9)
  calls$portal_z <- seq(9, 2)
  expect_lt(trait_verdict(calls)$z_correlation, 0)
  expect_equal(trait_verdict(NULL)$verdict, "INSUFFICIENT")
  expect_equal(trait_verdict(calls[1:4, ])$verdict, "INSUFFICIENT")
})

test_that("concordance_plot_svg escapes trait names and skips traits without data", {
  s <- data.frame(trait = c("A<B", "none"), verdict = c("PASS", "INSUFFICIENT"),
                  n_informative = c(10L, 0L), n_concordant = c(9L, 0L),
                  concordance = c(0.9, NA))
  svg <- concordance_plot_svg(s)
  expect_match(svg, "A&lt;B", fixed = TRUE)
  expect_false(grepl(">none<", svg, fixed = TRUE))
  expect_match(concordance_plot_svg(s[2, ]), "No trait", fixed = TRUE)
})

write_lines_tsv <- function(lines) {
  f <- tempfile(fileext = ".tsv")
  writeLines(lines, f)
  f
}

test_that("read_audit_manifest validates the trait list", {
  traits <- c("good", "flipped")
  m <- read_audit_manifest(write_lines_tsv(c("trait\tkp_phenotype\tsuggestions",
                                             "good\tHDL\tx", "flipped\t\t")), traits)
  expect_equal(m$kp_phenotype, c("HDL", ""))
  expect_named(m, c("trait", "kp_phenotype"))
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tkp_phenotype", "nope\tHDL")), traits),
               "not in the matrix")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tkp_phenotype", "good\tHDL;LDL")),
                                   traits), "one Portal phenotype")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tefo_id", "good\tEFO_1")), traits),
               "kp_phenotype")
  expect_error(read_audit_manifest(write_lines_tsv(c("trait\tkp_phenotype", "good\tHDL",
                                                     "good\tLDL")), traits), "Duplicate")
})

test_that("draft_audit_manifest writes suggestions and never overwrites", {
  out <- tempfile(fileext = ".tsv")
  ph <- data.frame(name = c("WaistC", "WaistCadjBMI"),
                   description = c("Waist circumference", "Waist circumference adj BMI"),
                   group = "", stringsAsFactors = FALSE)
  quietly(draft_audit_manifest(c("Waist circumference_male", "XINSG30"), out, ph))
  m <- utils::read.delim(out, colClasses = "character", na.strings = character(0))
  expect_equal(names(m), c("trait", "kp_phenotype", "suggestions"))
  expect_equal(m$kp_phenotype, c("WaistC", ""))
  expect_equal(m$suggestions[1], "WaistC (Waist circumference)")
  expect_error(draft_audit_manifest("x", out, ph), "already exists")
  expect_equal(nrow(read_audit_manifest(out, c("Waist circumference_male", "XINSG30"))), 2L)
})

test_that("audit_zscore_alignment labels every trait and writes its outputs", {
  fx <- audit_fixture()
  manifest <- data.frame(
    trait = c("good", "flipped", "scrambled", "unknown", "blank"),
    kp_phenotype = c("TRAITX", "TRAITX", "TRAITX", "NOPE", ""),
    stringsAsFactors = FALSE
  )
  prefix <- file.path(tempfile(), "audit")
  res <- quietly(audit_zscore_alignment(fx$z, fx$alignment, manifest, fx$portal, prefix,
                                        phenotype_names = c("TRAITX", "OTHER"),
                                        settings = list(ancestry = "EU", failed = "1:999")))
  s <- res$summary
  expect_equal(s$verdict[match(manifest$trait, s$trait)],
               c("PASS", "FLIPPED", "REVIEW", "INSUFFICIENT", "NOT_CHECKED"))
  expect_match(s$note[s$trait == "unknown"], "not a Portal phenotype")
  for (suffix in c(".csv", "_traits.csv", ".html")) {
    expect_true(file.exists(paste0(prefix, suffix)))
  }
  comp <- utils::read.csv(paste0(prefix, ".csv"))
  expect_equal(sort(unique(comp$trait)), c("flipped", "good", "scrambled"))
  expect_equal(nrow(comp), 24L)
  html <- paste(readLines(paste0(prefix, ".html")), collapse = "\n")
  for (t in manifest$trait) expect_match(html, paste0(">", t, "<"), fixed = TRUE)
  expect_match(html, "<details open>", fixed = TRUE)
  expect_match(html, "EU analysis", fixed = TRUE)
  expect_match(html, "could not be fetched", fixed = TRUE)
  expect_false(grepl("<script", html, fixed = TRUE))
})
