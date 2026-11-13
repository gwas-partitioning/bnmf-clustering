# =============================================================================
# audit_zscore_alignment.R
# =============================================================================
# Checks that each trait column of the pre-imputed z-score matrix
# (zmat_preImputed.rds) is oriented to the T2D risk allele, by comparing signs
# with genome-wide-significant associations in the NHGRI-EBI GWAS Catalog.
#
# The check is independent of your own summary-statistics files, so it can
# catch an allele mistake made while preprocessing them (for example, the
# wrong effect-allele column), which flips or scrambles a whole trait.
#
# The Catalog data come from its bulk "All associations with added ontology
# annotations" file (about 70 MB zipped, updated weekly):
#   https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip
# Pass the .zip (needs the unzip command) or the unzipped .tsv. The REST API
# is too slow for whole-trait queries (about 45 s per 500 records).
#
# Usage:
#   1. Draft the audit manifest (EFO terms suggested from each trait name):
#        Rscript scripts/audit_zscore_alignment.R --draft --catalog CATALOG.zip \
#          results/zmat_preImputed.rds audit_traits.tsv
#   2. Review audit_traits.tsv: correct efo_id (see the suggestions column;
#      separate several IDs with ";"), blank it to skip a trait, and
#      optionally add the GCST accession of the study the clustering used
#      and a reported_trait regular expression.
#   3. Run the audit:
#        Rscript scripts/audit_zscore_alignment.R --catalog CATALOG.zip \
#          results/zmat_preImputed.rds results/alignment_GWAS_summStats.csv \
#          audit_traits.tsv results/zscore_audit [--min-abs-z 2]
#      Writes results/zscore_audit.csv (one row per SNP-trait comparison),
#      results/zscore_audit_traits.csv (one row per trait), and
#      results/zscore_audit.html.
#
# Verdicts: PASS (>= 80% of informative SNPs agree in sign), FLIPPED (<= 20%),
# REVIEW (in between), INSUFFICIENT (< 5 informative SNPs), NOT_CHECKED (no
# EFO ID).
#
# From R, source scripts/qc_report.R before this file.
#
# Dependencies: data.table.
# =============================================================================

AUDIT_GWS_P <- 5e-8
AUDIT_BASES <- c("A", "C", "G", "T")
AUDIT_TERM_PATTERN <- "^[A-Za-z]+_[A-Za-z0-9]+$"   # EFO_0004340, OBA_VT0005416, MONDO_0005148
CATALOG_COLUMNS <- c(
  risk = "STRONGEST SNP-RISK ALLELE", snps = "SNPS", current = "SNP_ID_CURRENT",
  p = "P-VALUE", value = "OR or BETA", ci = "95% CI (TEXT)", mapped = "MAPPED_TRAIT",
  uri = "MAPPED_TRAIT_URI", accession = "STUDY ACCESSION", reported = "DISEASE/TRAIT"
)

audit_complement <- function(x) chartr("ACGT", "TGCA", x)

catalog_direction <- function(ci_text, value) {
  # The Catalog stores effect sizes unsigned. A beta's direction is written in
  # the CI text ("[0.01-0.02] unit increase"); without it, a value above 1 is
  # an odds ratio for the risk allele. Anything else (including a bare value
  # below 1, which could be an unsigned beta) gives no direction.
  ci_text <- tolower(ifelse(is.na(ci_text), "", ci_text))
  value <- suppressWarnings(as.numeric(value))
  ifelse(grepl("increase", ci_text, fixed = TRUE), 1,
  ifelse(grepl("decrease", ci_text, fixed = TRUE), -1,
  ifelse(is.finite(value) & value > 1, 1, NA_real_)))
}

parse_catalog_rows <- function(raw) {
  # Normalise rows of the bulk associations file (all columns character).
  missing_cols <- setdiff(CATALOG_COLUMNS, names(raw))
  if (length(missing_cols)) {
    stop("The GWAS Catalog file lacks column(s): ", paste(missing_cols, collapse = ", "),
         ". Use the 'All associations' file with added ontology annotations.")
  }
  col <- function(key) {
    x <- as.character(raw[[CATALOG_COLUMNS[[key]]]])
    x[is.na(x)] <- ""
    trimws(x)
  }
  snps <- col("snps")
  risk <- gsub("\\s+", "", col("risk"))
  single <- grepl("^rs[0-9]+$", snps)
  # "rs123-A": the effect allele follows the last "-" and must name this SNP.
  effect_allele <- toupper(sub("^.*-", "", risk))
  effect_allele[!grepl("-", risk, fixed = TRUE) | sub("-[^-]*$", "", risk) != snps |
                  !nzchar(effect_allele)] <- NA_character_
  current <- col("current")
  # "http://www.ebi.ac.uk/efo/EFO_1, http://purl.obolibrary.org/obo/MONDO_2"
  # becomes "EFO_1,MONDO_2".
  term_ids <- gsub("\\s", "", gsub("[^,]*/", "", col("uri")))
  data.frame(
    catalog_row = as.character(seq_along(snps)),
    accession = col("accession"),
    rsID = ifelse(single, snps, NA_character_),
    rsID_current = ifelse(single & grepl("^[0-9]+$", current), paste0("rs", current),
                          NA_character_),
    effect_allele = effect_allele,
    p_value = suppressWarnings(as.numeric(col("p"))),
    ci_text = col("ci"),
    effect_value = suppressWarnings(as.numeric(col("value"))),
    reported_trait = col("reported"),
    mapped_trait = col("mapped"),
    term_ids = term_ids,
    single_snp = single,
    stringsAsFactors = FALSE
  )
}

read_catalog_bulk <- function(path) {
  if (!file.exists(path)) stop("GWAS Catalog file not found: ", path)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install the data.table package to read the GWAS Catalog file.")
  }
  args <- list(sep = "\t", quote = "", colClasses = "character", na.strings = "",
               select = unname(CATALOG_COLUMNS), showProgress = FALSE)
  source_arg <- if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    if (!nzchar(Sys.which("unzip"))) stop("Reading a .zip needs the unzip command; unzip it and pass the .tsv.")
    list(cmd = paste("unzip -p", shQuote(path)))
  } else {
    list(file = path)
  }
  raw <- suppressWarnings(do.call(data.table::fread, c(source_arg, args)))
  out <- parse_catalog_rows(as.data.frame(raw))
  message(sprintf("Read %s GWAS Catalog associations from %s",
                  format(nrow(out), big.mark = ","), basename(path)))
  out
}

catalog_for_terms <- function(catalog, term_ids) {
  # Rows annotated with any of term_ids (whole IDs: EFO_1 does not match EFO_10).
  keys <- paste0(",", catalog$term_ids, ",")
  hit <- Reduce(`|`, lapply(term_ids, function(id) {
    grepl(paste0(",", id, ","), keys, fixed = TRUE)
  }), FALSE)
  catalog[hit, , drop = FALSE]
}

catalog_terms <- function(catalog) {
  # One row per ontology term in the file: efo_id, efo_trait, n (associations).
  has <- nzchar(catalog$term_ids)
  tab <- table(paste(catalog$term_ids[has], catalog$mapped_trait[has], sep = "\t"))
  parts <- strsplit(names(tab), "\t", fixed = TRUE)
  ids <- strsplit(vapply(parts, `[`, character(1), 1L), ",", fixed = TRUE)
  labels <- strsplit(vapply(parts, function(p) if (length(p) > 1L) p[2] else "", character(1)),
                     ", ", fixed = TRUE)
  # A label containing ", " cannot be paired with its ID; skip those rows.
  ok <- lengths(ids) == lengths(labels)
  long <- data.frame(efo_id = unlist(ids[ok]), efo_trait = unlist(labels[ok]),
                     n = rep(as.integer(tab)[ok], lengths(ids[ok])), stringsAsFactors = FALSE)
  if (!nrow(long)) return(data.frame(efo_id = character(0), efo_trait = character(0), n = integer(0)))
  key <- paste(long$efo_id, long$efo_trait, sep = "\t")
  n <- rowsum(long$n, key, reorder = FALSE)
  first <- long[match(rownames(n), key), c("efo_id", "efo_trait")]
  first$n <- as.integer(n[, 1])
  rownames(first) <- NULL
  first
}

# Common trait-name abbreviations, expanded before matching Catalog labels.
AUDIT_TRAIT_SYNONYMS <- c(
  "bmi" = "body mass index", "hdl" = "high density lipoprotein cholesterol",
  "ldl" = "low density lipoprotein cholesterol", "tg" = "triglyceride",
  "logtg" = "triglyceride", "sbp" = "systolic blood pressure",
  "dbp" = "diastolic blood pressure", "whr" = "waist hip ratio",
  "shbg" = "sex hormone binding globulin", "alt" = "alanine aminotransferase",
  "ast" = "aspartate aminotransferase", "alp" = "alkaline phosphatase",
  "ggt" = "gamma glutamyl transferase", "gamma glutamyltransferase" = "gamma glutamyl transferase",
  "crp" = "c reactive protein", "white blood cell" = "leukocyte",
  "red blood cell" = "erythrocyte", "oestradiol" = "estradiol", "haemoglobin" = "hemoglobin",
  "adjbmi" = "bmi adjusted"
)

audit_norm_label <- function(x) trimws(gsub("[^a-z0-9]+", " ", tolower(x)))

suggest_efo_terms <- function(terms, query, n = 3L) {
  # Up to n terms whose label matches query, best first:
  #   0 label and query have the same words
  #   1 label has every word of the query (and more)
  #   2 query has every word of the label
  # ignoring a trailing "measurement" and word order, and trying the query
  # both as given and with abbreviations expanded (BMI, HDL, adjBMI, ...).
  # More associations first within a tier.
  empty <- terms[0, c("efo_id", "efo_trait", "n")]
  q <- audit_norm_label(query)
  if (!nzchar(q) || !nrow(terms)) return(empty)
  expanded <- q
  for (abbr in names(AUDIT_TRAIT_SYNONYMS)) {
    expanded <- gsub(paste0("(^| )", abbr, "( |$)"),
                     paste0("\\1", AUDIT_TRAIT_SYNONYMS[[abbr]], "\\2"), expanded)
  }
  core <- sub(" measurement$", "", audit_norm_label(terms$efo_trait))
  core_words <- strsplit(core, " ", fixed = TRUE)
  tier_for <- function(q) {
    qw <- unique(strsplit(q, " ", fixed = TRUE)[[1]])
    # Only terms sharing a word with the query can match.
    tier <- rep(NA_integer_, length(core))
    cand <- which(vapply(core_words, function(w) any(w %in% qw), logical(1)))
    has_query <- vapply(core_words[cand], function(w) all(qw %in% w), logical(1))
    in_query <- vapply(core_words[cand], function(w) all(w %in% qw), logical(1)) &
      nchar(core[cand]) >= 4L
    tier[cand] <- ifelse(has_query & in_query, 0L, ifelse(has_query, 1L,
                         ifelse(in_query, 2L, NA_integer_)))
    tier
  }
  tier <- pmin(tier_for(q), tier_for(expanded), na.rm = TRUE)
  keep <- which(!is.na(tier))
  if (!length(keep)) return(empty)
  keep <- keep[order(tier[keep], -terms$n[keep], nchar(core[keep]))]
  out <- utils::head(terms[keep, c("efo_id", "efo_trait", "n")], n)
  rownames(out) <- NULL
  out
}

prepare_alignment <- function(alignment, z) {
  # One alignment row per matrix row, in matrix order.
  needed <- c("ChrPos", "rsID", "REF", "ALT", "Risk_Allele")
  missing_cols <- setdiff(needed, names(alignment))
  if (length(missing_cols)) {
    stop("The alignment file lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  # Runs from earlier pipeline versions name matrix rows "chr1:123" while the
  # alignment file has "1:123"; match on the position without the prefix.
  strip_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
  idx <- match(strip_chr(rownames(z)), strip_chr(alignment$ChrPos))
  if (anyNA(idx)) {
    stop(sum(is.na(idx)), " matrix row(s) have no alignment row (first: ",
         rownames(z)[is.na(idx)][1], ").")
  }
  a <- as.data.frame(alignment, stringsAsFactors = FALSE)[idx, needed]
  a[] <- lapply(a, as.character)
  a$ChrPos <- rownames(z)
  a[c("REF", "ALT", "Risk_Allele")] <- lapply(a[c("REF", "ALT", "Risk_Allele")], toupper)
  bad_risk <- is.na(a$Risk_Allele) | !(a$Risk_Allele == a$REF | a$Risk_Allele == a$ALT)
  bad_risk[is.na(bad_risk)] <- TRUE
  if (any(bad_risk)) {
    stop(sum(bad_risk), " variant(s) have a Risk_Allele that is neither REF nor ALT (first: ",
         a$ChrPos[bad_risk][1], ").")
  }
  a$rsID[is.na(a$rsID) | !grepl("^rs[0-9]+$", a$rsID)] <- NA_character_
  if (all(is.na(a$rsID))) {
    stop("The alignment file has no rsIDs; the audit matches the Catalog by rsID. ",
         "rsIDs come from LDlink, so they are blank when LD pruning is off.")
  }
  a$palindromic <- paste0(a$REF, a$ALT) %in% c("AT", "TA", "CG", "GC")
  rownames(a) <- NULL
  a
}

align_catalog_to_matrix <- function(catalog, alignment, accession = "",
                                    reported_trait = "", p_threshold = AUDIT_GWS_P) {
  # Keep usable genome-wide-significant records, match them to matrix rows by
  # rsID, and express each record's direction for the T2D risk allele.
  catalog$catalog_direction <- catalog_direction(catalog$ci_text, catalog$effect_value)
  keep <- catalog$single_snp & !is.na(catalog$rsID) & is.finite(catalog$p_value) &
    catalog$p_value <= p_threshold & !is.na(catalog$catalog_direction)
  if (nzchar(reported_trait)) {
    keep <- keep & grepl(reported_trait, catalog$reported_trait, ignore.case = TRUE,
                         perl = TRUE)
  }
  kept <- catalog[keep, , drop = FALSE]
  # Also match on the current rsID when dbSNP has merged the reported one.
  merged <- kept[!is.na(kept$rsID_current) & kept$rsID_current != kept$rsID, , drop = FALSE]
  merged$rsID <- merged$rsID_current
  kept <- rbind(kept, merged)

  usable <- alignment[!is.na(alignment$rsID) & !alignment$palindromic, , drop = FALSE]
  d <- merge(usable, kept, by = "rsID")
  d <- d[!duplicated(paste(d$catalog_row, d$ChrPos)), , drop = FALSE]

  ea <- d$effect_allele
  base <- !is.na(ea) & ea %in% AUDIT_BASES
  direct <- base & (ea == d$REF | ea == d$ALT)
  flipped_strand <- base & !direct &
    (audit_complement(ea) == d$REF | audit_complement(ea) == d$ALT)
  matched <- direct | flipped_strand
  d <- d[matched, , drop = FALSE]
  d$strand_complemented <- flipped_strand[matched]
  d$aligned_allele <- ifelse(d$strand_complemented, audit_complement(d$effect_allele),
                             d$effect_allele)
  d$expected_sign <- d$catalog_direction * ifelse(d$aligned_allele == d$Risk_Allele, 1, -1)
  d$same_study <- nzchar(accession) & !is.na(d$accession) & d$accession == accession
  d <- d[order(match(d$ChrPos, alignment$ChrPos), d$p_value), , drop = FALSE]
  rownames(d) <- NULL
  d
}

call_snp_trait <- function(aligned, z_values, min_abs_z = 2) {
  # One call per matrix variant: the smallest-P record from the clustering
  # study when there is one, otherwise from any study. Records in that pool
  # that disagree in direction make the call a catalog_conflict.
  if (!nrow(aligned)) return(NULL)
  rows <- lapply(split(seq_len(nrow(aligned)), aligned$ChrPos), function(i) {
    g <- aligned[i, , drop = FALSE]
    pool <- if (any(g$same_study)) g[g$same_study, , drop = FALSE] else g
    best <- pool[which.min(pool$p_value), , drop = FALSE]
    best$n_catalog_records <- nrow(g)
    best$catalog_conflict <- length(unique(pool$expected_sign)) > 1L
    best
  })
  calls <- do.call(rbind, rows)
  calls <- calls[order(match(calls$ChrPos, names(z_values))), , drop = FALSE]
  calls$matrix_z <- unname(z_values[calls$ChrPos])
  calls$status <- ifelse(calls$catalog_conflict, "catalog_conflict",
                  ifelse(!is.finite(calls$matrix_z), "missing_in_matrix",
                  ifelse(abs(calls$matrix_z) < min_abs_z, "weak_matrix_z",
                  ifelse(sign(calls$matrix_z) == calls$expected_sign,
                         "concordant", "discordant"))))
  calls$informative <- calls$status %in% c("concordant", "discordant")
  calls$concordant <- ifelse(calls$informative, calls$status == "concordant", NA)
  rownames(calls) <- NULL
  calls
}

AUDIT_VERDICT_ORDER <- c("FLIPPED", "REVIEW", "INSUFFICIENT", "PASS", "NOT_CHECKED")

trait_verdict <- function(calls, accession = "", min_n = 5L, pass_at = 0.8,
                          flip_at = 0.2) {
  informative <- if (is.null(calls)) logical(0) else calls$informative
  n <- sum(informative)
  k <- if (n) sum(calls$concordant[informative]) else 0L
  concordance <- if (n) k / n else NA_real_
  same <- if (n) informative & calls$same_study else logical(0)
  n_same <- sum(same)
  verdict <- if (n < min_n) {
    "INSUFFICIENT"
  } else if (concordance >= pass_at) {
    "PASS"
  } else if (concordance <= flip_at) {
    "FLIPPED"
  } else {
    "REVIEW"
  }
  data.frame(
    verdict = verdict,
    n_catalog_snps = if (is.null(calls)) 0L else nrow(calls),
    n_informative = n,
    n_concordant = k,
    concordance = concordance,
    binom_p = if (n) stats::binom.test(k, n)$p.value else NA_real_,
    n_same_study = n_same,
    same_study_concordance = if (n_same) mean(calls$concordant[same]) else NA_real_,
    study_differs = if (nzchar(accession)) n_same == 0L else NA,
    stringsAsFactors = FALSE
  )
}

AUDIT_COLORS <- c(PASS = "#2f7d4f", REVIEW = "#b26b00", FLIPPED = "#b3261e")

audit_setting <- function(settings, name, default) {
  if (is.null(settings[[name]])) default else settings[[name]]
}

concordance_plot_svg <- function(summary) {
  # Dot plot: one row per trait with informative SNPs, x = share concordant,
  # dot area grows with the number of informative SNPs.
  s <- summary[summary$n_informative > 0 & is.finite(summary$concordance), , drop = FALSE]
  if (!nrow(s)) return("<p class=\"muted\">No trait had informative SNPs to plot.</p>")
  s <- s[order(s$concordance, s$trait), , drop = FALSE]
  left <- 240; plot_w <- 480; right <- 60; top <- 16; row_h <- 18
  width <- left + plot_w + right
  bottom <- top + nrow(s) * row_h
  height <- bottom + 40
  f <- function(v) formatC(v, format = "f", digits = 1)
  x <- function(v) left + v * plot_w
  y <- top + (seq_len(nrow(s)) - 0.5) * row_h
  r <- pmin(8, 2.5 + sqrt(s$n_informative) / 1.5)
  col <- unname(AUDIT_COLORS[s$verdict])
  col[is.na(col)] <- "#7a8794"
  ticks <- c(0, 0.2, 0.5, 0.8, 1)
  paste(c(
    sprintf(paste0('<svg viewBox="0 0 %d %s" width="100%%" role="img" ',
                   'aria-label="Sign concordance by trait" ',
                   'style="max-width:%dpx;font-size:11px;font-family:inherit">'),
            width, f(height), width),
    sprintf('<rect x="%s" y="%d" width="%s" height="%s" fill="#b3261e" fill-opacity="0.08"/>',
            f(x(0)), top, f(0.2 * plot_w), f(bottom - top)),
    sprintf('<rect x="%s" y="%d" width="%s" height="%s" fill="#b26b00" fill-opacity="0.07"/>',
            f(x(0.2)), top, f(0.6 * plot_w), f(bottom - top)),
    sprintf(paste0('<line x1="%1$s" x2="%1$s" y1="%2$d" y2="%3$s" stroke="#c9ced8"/>',
                   '<text x="%1$s" y="%4$s" text-anchor="middle" fill="#6a7486">%5$s</text>'),
            f(x(ticks)), top, f(bottom), f(bottom + 14), format(ticks)),
    sprintf(paste0('<text x="%s" y="%s" text-anchor="middle" fill="#3b4658">',
                   'Share of informative SNPs whose matrix sign matches the Catalog</text>'),
            f(x(0.5)), f(bottom + 32)),
    sprintf('<text x="%d" y="%s" text-anchor="end" dominant-baseline="middle" fill="#1c2330">%s</text>',
            left - 8, f(y), qc_html_escape(s$trait)),
    sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s"><title>%s: %d of %d concordant (%s)</title></circle>',
            f(x(s$concordance)), f(y), f(r), col, qc_html_escape(s$trait),
            as.integer(s$n_concordant), as.integer(s$n_informative), s$verdict),
    sprintf('<text x="%s" y="%s" dominant-baseline="middle" fill="#6a7486">n=%d</text>',
            f(x(1) + 12), f(y), as.integer(s$n_informative)),
    "</svg>"
  ), collapse = "\n")
}

write_audit_report <- function(summary, comparisons, file, settings = list()) {
  if (!exists("qc_html_table", mode = "function")) {
    stop("Source scripts/qc_report.R before writing the audit report.")
  }
  min_abs_z <- audit_setting(settings, "min_abs_z", 2)
  p_threshold <- audit_setting(settings, "p_threshold", AUDIT_GWS_P)
  summary <- summary[order(match(summary$verdict, AUDIT_VERDICT_ORDER), summary$trait), ,
                     drop = FALSE]
  count <- table(factor(summary$verdict, levels = AUDIT_VERDICT_ORDER))
  traits_with <- function(v) qc_name_list(summary$trait[summary$verdict == v])

  flags <- character(0)
  if (count[["FLIPPED"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) look flipped: the matrix sign disagrees with the Catalog at most informative SNPs (%s). Check which column was used as the effect allele when these GWAS were formatted.",
    count[["FLIPPED"]], traits_with("FLIPPED")))
  if (count[["REVIEW"]] > 0) flags <- c(flags, sprintf(
    "%d trait(s) have mixed agreement (%s). Possible causes: a partial allele mix-up, or an EFO term that pools a differently defined phenotype (for example adjusted and unadjusted). Their SNP tables are open below.",
    count[["REVIEW"]], traits_with("REVIEW")))
  unknown <- summary$trait[grepl("not in the Catalog file", summary$note, fixed = TRUE)]
  if (length(unknown)) flags <- c(flags, sprintf(
    "Some EFO IDs were not found in the Catalog file, for: %s. Check them against the Catalog's trait pages.",
    qc_name_list(unknown)))
  differs <- summary$trait[summary$verdict %in% c("PASS", "REVIEW", "FLIPPED") &
                             summary$study_differs %in% TRUE]
  if (length(differs)) flags <- c(flags, sprintf(
    "No informative Catalog record came from the clustering study for: %s. Their verdicts rest on other studies of the trait.",
    qc_name_list(differs)))

  check_block <- if (length(flags)) {
    paste0("<div class=\"checks\"><h2>Things to check</h2><ul>",
           paste0("<li>", qc_html_escape(flags), "</li>", collapse = ""), "</ul></div>")
  } else {
    "<div class=\"checks ok\"><h2>Things to check</h2><p>No flags were raised.</p></div>"
  }

  summary_cols <- c("trait", "verdict", "n_informative", "n_concordant", "concordance",
                    "binom_p", "n_catalog_snps", "n_same_study", "same_study_concordance",
                    "study_differs", "efo_id", "accession", "reported_trait", "note")
  detail_cols <- c("rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "effect_allele",
                   "strand_complemented", "catalog_direction", "expected_sign", "matrix_z",
                   "status", "p_value", "accession", "same_study", "catalog_row",
                   "reported_trait")
  details <- vapply(seq_len(nrow(summary)), function(i) {
    s <- summary[i, ]
    rows <- if (is.null(comparisons)) NULL else comparisons[comparisons$trait == s$trait, , drop = FALSE]
    if (is.null(rows) || !nrow(rows)) return("")
    sprintf("<details%s><summary><strong>%s</strong> &middot; %s (%d of %d informative SNPs concordant)</summary>%s</details>",
            if (s$verdict %in% c("FLIPPED", "REVIEW")) " open" else "",
            qc_html_escape(s$trait), s$verdict, as.integer(s$n_concordant),
            as.integer(s$n_informative), qc_html_table(rows[detail_cols], max_rows = 500L))
  }, character(1))

  meta <- sprintf("%s variants &times; %s traits &middot; Catalog file %s &middot; generated %s",
                  audit_setting(settings, "n_variants", "?"),
                  audit_setting(settings, "n_traits", nrow(summary)),
                  qc_html_escape(audit_setting(settings, "catalog_label", "not recorded")),
                  format(Sys.time(), "%Y-%m-%d %H:%M"))
  counts_line <- paste(sprintf("%s: %d", names(count), as.integer(count)), collapse = " &middot; ")
  methods <- sprintf(paste(
    "<p>For each trait, every GWAS Catalog association annotated with the trait's ontology term(s) was taken from any study,",
    "keeping single-SNP records with P &le; %s and a direction: a beta whose CI text says increase or decrease, or an",
    "unsigned value above 1, read as an odds ratio for the reported allele. Records were matched to matrix variants by rsID",
    "(or by the current rsID when dbSNP has merged the reported one). A/T and C/G variants were excluded; an effect allele",
    "matching only the complement of REF/ALT was read on the opposite strand. Each record's direction was expressed for the",
    "T2D risk allele.</p>",
    "<p>Each variant gets one call: the smallest-P record from the clustering study when its accession is given and present,",
    "otherwise from any study. When those records disagree in direction the call is a Catalog conflict and is not counted.",
    "A call is informative when the matrix |z| &ge; %s. Verdicts need at least 5 informative variants: PASS when at least 80%%",
    "agree in sign, FLIPPED when at most 20%% do, REVIEW in between. The binomial P tests agreement against 50%%.</p>",
    "<p>Limitations: an ontology term can pool related phenotypes whose directions differ (adjusted versus unadjusted,",
    "sex-specific); the reported_trait filter narrows a term. Catalog records are curated top hits, so traits with few",
    "significant loci among these variants may be INSUFFICIENT. An unsigned beta above 1 with no direction text would be",
    "misread as an odds ratio. Only the sign is checked, not the magnitude.</p>"),
    format(p_threshold), format(min_abs_z))

  css <- paste(
    "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f6f7f9;color:#1c2330;}",
    "main{max-width:1100px;margin:0 auto;padding:24px 16px 64px;}",
    "h1{font-size:24px;margin:8px 0 4px;} h2{font-size:18px;margin:32px 0 8px;}",
    "p{line-height:1.5;margin:6px 0;} .muted{color:#6a7486;font-size:13px;}",
    ".meta{color:#6a7486;font-size:13px;margin-bottom:16px;}",
    ".checks{background:#fff7e6;border:1px solid #f0c36d;border-radius:8px;padding:4px 16px 12px;}",
    ".checks.ok{background:#eef8f0;border-color:#9fd3aa;} .checks h2{margin-top:12px;} .checks li{margin:6px 0;line-height:1.45;}",
    ".tablewrap{overflow-x:auto;background:#fff;border:1px solid #dfe3ea;border-radius:8px;}",
    "table{border-collapse:collapse;width:100%;font-size:13px;}",
    "th,td{padding:6px 10px;border-bottom:1px solid #eceff4;text-align:left;white-space:nowrap;}",
    "th{background:#f0f2f6;font-weight:600;} tr:last-child td{border-bottom:none;}",
    ".num{text-align:right;font-variant-numeric:tabular-nums;}",
    "details{margin:8px 0;} summary{cursor:pointer;padding:4px 0;}",
    ".plot{background:#fff;border:1px solid #dfe3ea;border-radius:8px;padding:12px;}",
    sep = "\n"
  )
  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Z-score alignment audit</title>",
    paste0("<style>", css, "</style></head><body><main>"),
    "<h1>Z-score alignment audit</h1>",
    paste0("<div class=\"meta\">", meta, "</div>"),
    check_block,
    "<h2>Results by trait</h2>",
    paste0("<p class=\"muted\">", counts_line, "</p>"),
    qc_html_table(summary[summary_cols], max_rows = 1000L),
    "<h2>Sign concordance by trait</h2>",
    paste0("<div class=\"plot\">", concordance_plot_svg(summary), "</div>"),
    "<p class=\"muted\">Red band: at most 20% agree (FLIPPED). Amber band: 20&ndash;80% (REVIEW).</p>",
    "<h2>Variant calls by trait</h2>",
    details[nzchar(details)],
    "<h2>Method</h2>",
    methods,
    "</main></body></html>"
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, file, useBytes = TRUE)
  message(sprintf("Audit report written to %s (%d item(s) to check).", file, length(flags)))
  invisible(flags)
}

AUDIT_COMPARISON_COLS <- c(
  "trait", "rsID", "ChrPos", "REF", "ALT", "Risk_Allele", "catalog_row", "accession",
  "same_study", "reported_trait", "p_value", "effect_allele", "strand_complemented",
  "catalog_direction", "expected_sign", "matrix_z", "n_catalog_records", "informative",
  "concordant", "status"
)

audit_read_matrix <- function(path) {
  z <- readRDS(path)
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z) || !is.numeric(z) || is.null(rownames(z)) || is.null(colnames(z)) ||
      anyDuplicated(rownames(z)) || anyDuplicated(colnames(z))) {
    stop(path, " is not a numeric matrix with unique row and column names.")
  }
  z
}

audit_split_ids <- function(x) {
  ids <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  ids[nzchar(ids)]
}

read_audit_manifest <- function(path, traits) {
  m <- utils::read.delim(path, colClasses = "character", na.strings = character(0),
                         quote = "", check.names = FALSE)
  missing_cols <- setdiff(c("trait", "efo_id"), names(m))
  if (length(missing_cols)) {
    stop(path, " lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  for (col in c("accession", "reported_trait")) if (!col %in% names(m)) m[[col]] <- ""
  m <- m[c("trait", "efo_id", "accession", "reported_trait")]
  m[] <- lapply(m, function(x) trimws(ifelse(is.na(x), "", x)))

  if (anyDuplicated(m$trait)) {
    stop("Duplicate traits in ", path, ": ", qc_name_list(unique(m$trait[duplicated(m$trait)])))
  }
  unknown <- setdiff(m$trait, traits)
  if (length(unknown)) stop("Traits not in the matrix: ", paste(unknown, collapse = ", "))
  ids <- unlist(lapply(m$efo_id, audit_split_ids))
  bad_ids <- ids[!grepl(AUDIT_TERM_PATTERN, ids)]
  if (length(bad_ids)) stop("Invalid EFO ID(s): ", paste(bad_ids, collapse = ", "))
  bad_acc <- m$accession[nzchar(m$accession) & !grepl("^GCST[0-9]+$", m$accession)]
  if (length(bad_acc)) {
    stop("Invalid accession(s) (expected GCST followed by digits): ", paste(bad_acc, collapse = ", "))
  }
  for (re in m$reported_trait[nzchar(m$reported_trait)]) {
    ok <- tryCatch({ suppressWarnings(grepl(re, "", perl = TRUE)); TRUE },
                   error = function(e) FALSE)
    if (!ok) stop("Invalid reported_trait regular expression: ", re)
  }
  m
}

draft_audit_manifest <- function(traits, out_file, catalog) {
  # Suggest ontology terms from each trait name. The result needs human review.
  if (file.exists(out_file)) {
    stop(out_file, " already exists; edit it, or choose another file name.")
  }
  terms <- catalog_terms(catalog)
  rows <- lapply(traits, function(trait) {
    hits <- suggest_efo_terms(terms, gsub("[_.]+", " ", trait))
    found <- nrow(hits) > 0L
    message(sprintf("%s: %s", trait,
                    if (found) sprintf("%s (%s)", hits$efo_id[1], hits$efo_trait[1]) else "no match"))
    data.frame(
      trait = trait,
      efo_id = if (found) hits$efo_id[1] else "",
      accession = "", reported_trait = "",
      suggestions = if (found) paste(sprintf("%s (%s, n=%d)", hits$efo_id, hits$efo_trait, hits$n),
                                     collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(out)
}

audit_zscore_alignment <- function(z, alignment, manifest, catalog, output_prefix = NULL,
                                   min_abs_z = 2, p_threshold = AUDIT_GWS_P,
                                   catalog_label = "") {
  a <- prepare_alignment(alignment, z)
  known_ids <- unique(unlist(strsplit(unique(catalog$term_ids), ",", fixed = TRUE)))
  # Only records for the matrix's variants can matter.
  matrix_rsids <- a$rsID[!is.na(a$rsID)]
  local <- catalog[catalog$rsID %in% matrix_rsids | catalog$rsID_current %in% matrix_rsids, ,
                   drop = FALSE]

  summaries <- list()
  comparisons <- list()
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, ]
    ids <- audit_split_ids(row$efo_id)
    note <- ""
    if (!length(ids)) {
      v <- trait_verdict(NULL)
      v$verdict <- "NOT_CHECKED"
      note <- "No EFO ID given."
    } else {
      missing_ids <- setdiff(ids, known_ids)
      if (length(missing_ids)) {
        note <- sprintf("EFO ID(s) not in the Catalog file: %s.", paste(missing_ids, collapse = ", "))
      }
      aligned <- align_catalog_to_matrix(catalog_for_terms(local, ids), a, row$accession,
                                         row$reported_trait, p_threshold)
      calls <- call_snp_trait(aligned, z[, row$trait], min_abs_z)
      v <- trait_verdict(calls, row$accession)
      if (!is.null(calls)) {
        comparisons[[length(comparisons) + 1L]] <- cbind(trait = row$trait, calls)
      }
    }
    message(sprintf("%s: %s (%d of %d informative SNPs concordant)", row$trait, v$verdict,
                    as.integer(v$n_concordant), as.integer(v$n_informative)))
    summaries[[i]] <- cbind(row, v, note = note, stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  comparisons <- if (length(comparisons)) do.call(rbind, comparisons)[AUDIT_COMPARISON_COLS]

  if (!is.null(output_prefix)) {
    dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(if (is.null(comparisons)) data.frame() else comparisons,
                     paste0(output_prefix, ".csv"), row.names = FALSE)
    utils::write.csv(summary, paste0(output_prefix, "_traits.csv"), row.names = FALSE)
    write_audit_report(summary, comparisons, paste0(output_prefix, ".html"),
                       settings = list(min_abs_z = min_abs_z, p_threshold = p_threshold,
                                       n_variants = nrow(z), n_traits = nrow(manifest),
                                       catalog_label = catalog_label))
  }
  invisible(list(summary = summary, comparisons = comparisons))
}

if (sys.nframe() == 0L) {
  # Rscript writes spaces in the script path as "~+~".
  script_path <- gsub("~+~", " ", sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), fixed = TRUE)
  source(file.path(dirname(normalizePath(script_path)), "qc_report.R"))

  usage <- paste(
    "Usage:",
    "  Rscript audit_zscore_alignment.R --draft --catalog CATALOG.zip zmat_preImputed.rds audit_traits.tsv",
    "  Rscript audit_zscore_alignment.R --catalog CATALOG.zip zmat_preImputed.rds \\",
    "    alignment_GWAS_summStats.csv audit_traits.tsv output_prefix [--min-abs-z 2]",
    "CATALOG.zip: gwas-catalog-associations_ontology-annotated-full.zip from",
    "  https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/ (or the unzipped .tsv).",
    "See the header of this script for the manifest columns.", sep = "\n")
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat(usage, "\n")
    quit(status = 0)
  }
  take_value <- function(flag) {
    at <- which(args == flag)
    if (!length(at)) return(NULL)
    if (at[1] == length(args)) { cat(usage, "\n"); quit(status = 1) }
    value <- args[at[1] + 1L]
    args <<- args[-c(at[1], at[1] + 1L)]
    value
  }
  catalog_path <- take_value("--catalog")
  min_abs_z_arg <- take_value("--min-abs-z")
  min_abs_z <- if (is.null(min_abs_z_arg)) 2 else suppressWarnings(as.numeric(min_abs_z_arg))
  draft <- "--draft" %in% args
  args <- setdiff(args, "--draft")
  if (is.null(catalog_path) || !is.finite(min_abs_z) || min_abs_z < 0 ||
      length(args) != if (draft) 2L else 4L) {
    cat(usage, "\n")
    quit(status = 1)
  }
  catalog_label <- sprintf("%s (modified %s)", basename(catalog_path),
                           format(file.info(catalog_path)$mtime, "%Y-%m-%d"))
  if (draft) {
    z <- audit_read_matrix(args[1])
    if (file.exists(args[2])) stop(args[2], " already exists; edit it, or choose another file name.")
    draft_audit_manifest(colnames(z), args[2], read_catalog_bulk(catalog_path))
    message("Wrote ", args[2], ". Review efo_id for every trait before running the audit.")
  } else {
    z <- audit_read_matrix(args[1])
    alignment <- utils::read.csv(args[2], colClasses = "character", check.names = FALSE)
    manifest <- read_audit_manifest(args[3], colnames(z))
    result <- audit_zscore_alignment(z, alignment, manifest, read_catalog_bulk(catalog_path),
                                     args[4], min_abs_z = min_abs_z,
                                     catalog_label = catalog_label)
    print(table(factor(result$summary$verdict, levels = AUDIT_VERDICT_ORDER)))
  }
}
