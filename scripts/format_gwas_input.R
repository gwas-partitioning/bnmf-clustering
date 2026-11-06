# =============================================================================
# format_gwas_input.R
# =============================================================================
# Converts raw GWAS summary statistics into the format the bNMF pipeline reads
# (choose_variants_2025.R, prep_bNMF_2025.R). One gzipped, tab-delimited output
# file per GWAS: <trait>_<study>_formatted.txt.gz
#
# QUICK START
#   1. Set the paths in section 1 (INPUT_DIR, OUTPUT_DIR, rsID map, chain file).
#   2. Add one row per raw GWAS to the MANIFEST (section 2).
#   3. Run:  Rscript format_gwas_input.R
#   4. Read the log: check the "[auto]" lines show the right columns were
#      picked (especially ea / nea / beta), and see how many variants were
#      dropped at each step. If a column is not detected, add it to
#      COL_OVERRIDES (section 3).
#   5. Point the pipeline manifest's full_path at the output files.
#
# WHICH MODE? (controlled by the rsID map in section 1)
#   The pipeline accepts the output of either mode for both discovery and trait
#   GWAS: it tracks the effect allele (EA) itself and does not need files to
#   share one REF/ALT order.
#   Without the map (both map options NULL, the default) -> output is
#       SNP (CHR:POS), EA, NEA, BETA, SE, P_VALUE, N, in the input's own allele
#       order. Keeps variants that are absent from the 1000 Genomes map.
#   With the map -> each variant is aligned to the map's REF/ALT. Output adds
#       VAR_ID (CHR_POS_REF_ALT), rsID, REF, ALT, with EA = ALT and BETA aligned
#       to ALT. Use this when you want reference-ordered IDs and rsIDs in the
#       files themselves; variants absent from the map are dropped.
#   The map is built once with generate_varid_to_rsid_map_file.R. A script run
#   uses one mode for all rows.
#
# WHAT IT DOES
#   Detects columns -> optional hg38->hg19 liftover -> keeps biallelic SNPs
#   (ACGT) -> (main mode) aligns to the map's REF/ALT, flipping BETA so it is
#   ALT-aligned -> writes the output. ORs are converted to log(OR) when
#   odds_ratio = TRUE. Existing outputs are skipped unless OVERWRITE <- TRUE.
#
# ASSUMPTIONS TO CHECK FOR A NEW FILE
#   - Effect-allele column really is the allele BETA refers to (A1/ALT/etc. are
#     guesses from the header name; verify against the study's README).
#   - BETA is on the trait scale (log-OR for binary traits, not OR).
#   - Files must include a non-effect allele (or a VAR_ID override, see below).
#   - N: use the file's per-variant N when present; otherwise set N in the
#     manifest. For binary traits use total or effective N, not cases only.
#
# Dependencies (CRAN): data.table, dplyr, readr
# Optional (hg38 input only): rtracklayer, GenomicRanges, IRanges, S4Vectors
# =============================================================================

library(data.table)
library(dplyr)
library(readr)


# =============================================================================
# 1. PATHS  (edit these)
# =============================================================================

INPUT_DIR  <- "/path/to/raw_gwas"         # folder with the raw GWAS files
OUTPUT_DIR <- "/path/to/formatted_gwas"   # where outputs are written
OVERWRITE  <- FALSE                        # TRUE = regenerate existing outputs

# rsID map (optional; leave BOTH as NULL to keep the input's allele order).
# To align to the 1000 Genomes REF/ALT, give ONE of the two; both come from
# generate_varid_to_rsid_map_file.R.
#   RSID_MAP_DIR : folder of headerless per-chromosome files chr1.txt, chr2.txt, ...
#                  (columns: hg19_posID, rsID, ref_allele, alt_allele). Fast.
#   RSID_MAP_FILE: single 2-column file VAR_ID (CHR_POS_REF_ALT), rsID (.gz ok).
RSID_MAP_DIR  <- NULL   # e.g. "/path/to/rsid_map_by_chr"
RSID_MAP_FILE <- NULL   # e.g. "/path/to/varid_rsid_map.txt.gz"

# Only needed if any manifest row has build = "hg38":
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
# (unzip first; rtracklayer reads the uncompressed .chain)
CHAIN_FILE <- NULL      # e.g. "/path/to/hg38ToHg19.over.chain"


# =============================================================================
# 2. MANIFEST  (one row per raw GWAS; replace the example rows)
# =============================================================================
#   file       : file name inside INPUT_DIR (.gz ok)
#   study      : study ID, used in the output file name
#   trait      : short trait label, used in the output file name and to link
#                COL_OVERRIDES entries
#   build      : "hg19" or "hg38" (hg38 is lifted over to hg19)
#   N          : sample size if the file has no N column; NA = use the file's N
#   odds_ratio : TRUE if the effect column holds odds ratios (converted to log)
MANIFEST <- data.frame(
  file       = c("example_trait_gwas.tsv.gz", "example_disease_gwas.txt.gz"),
  study      = c("StudyA",                    "StudyB"),
  trait      = c("HDL",                       "T2D"),
  build      = c("hg19",                      "hg38"),
  N          = c(NA_real_,                    NA_real_),
  odds_ratio = c(FALSE,                       FALSE),
  stringsAsFactors = FALSE
)


# =============================================================================
# 3. COL_OVERRIDES  (optional; only for files where auto-detection fails)
# =============================================================================
# Name = the manifest trait. Value = standard_name = "column name in your file".
# Standard names: chr, pos, ea, nea, beta, se, p, log10p, n, eaf
# Special names:
#   chrpos = one column holding "CHR:POS" (split into chr and pos)
#   varid  = column holding CHR_POS_REF_ALT; used to derive NEA if the file
#            has no non-effect-allele column
# Only list the columns that need fixing; the rest are still auto-detected.
COL_OVERRIDES <- list(
  # HDL = c(
  #   chr = "chromosome", pos = "base_pair_location",
  #   ea  = "effect_allele", nea = "other_allele",
  #   beta = "beta", se = "standard_error", p = "p_value", n = "n"
  # ),
  # T2D = c(chrpos = "ChrPos", varid = "VAR_ID", ea = "EA")
)


# =============================================================================
# AUTO-DETECTION CANDIDATES
# Common column name variants for each field.
# Add your file's column names here if auto-detection misses them.
# =============================================================================
CANDIDATES <- list(
  chr    = c("chromosome", "CHR", "CHROM", "#CHROM", "chr", "Chromosome", "chrom",
             "#chromosome"),
  pos    = c("base_pair_location", "POS", "BP", "GENPOS", "position", "pos",
             "Position", "base_pair_location_grch37", "POS_HG19", "bp"),
  ea     = c("effect_allele", "Effect_Allele", "EA", "A1", "Allele1", "allele1",
             "ALLELE1", "tested_allele", "alt", "ALT", "effect_allele_PH"),
  nea    = c("other_allele", "non_effect_allele", "NEA", "A2", "Allele2",
             "allele2", "ALLELE2", "reference_allele", "nea", "ref", "REF", "ALLELE0", "Allele0", "A0"),
  beta   = c("beta", "BETA", "Effect", "effect", "b", "EFFECT", "es",
             "effect_size"),
  se     = c("standard_error", "SE", "StdErr", "stderr", "se",
             "Standard_Error", "SE_EUR", "standard_error_EUR"),
  p      = c("p_value", "P_VALUE", "p.value", "pvalue", "P", "p", "PVALUE",
             "pval", "p_value_EUR", "p.value.exposure","P-value", "Pval", "Pvalue",
             "P_BOLT_LMM", "P_BOLT_LMM_INF"),
  log10p = c("LOG10P", "log10p", "MLOG10P", "mlog10p", "neg_log10_p",
             "neglog10p", "NEGLOG10P"),
  n      = c("n", "N", "N_PH", "sample_size", "obs_ct", "OBS_CT",
             "TotalSampleSize", "N_eff"),   # N_cases deliberately excluded (wrong N for binary traits)
  eaf    = c("effect_allele_frequency", "EAF", "AF", "A1_FREQ", "freq1",
             "eaf", "AF_EUR", "effect_allele_freq", "A1FREQ")
)


# =============================================================================
# FUNCTIONS
# =============================================================================

# Detect the first matching column name from a list of candidates.
detect_col <- function(headers, candidates, label, required = TRUE) {
  # Case-insensitive; candidate order sets priority, and the file's actual header is returned
  hit   <- match(tolower(candidates), tolower(headers))
  found <- headers[hit[!is.na(hit)]]
  if (length(found) > 0) {
    message(sprintf("    [auto] %-8s-> '%s'", label, found[1]))
    return(found[1])
  }
  if (required) {
    stop(sprintf(
      "Cannot detect '%s' column. Headers present: %s\n  Add an override to COL_OVERRIDES.",
      label, paste(headers, collapse = ", ")
    ))
  }
  NULL
}

# Resolve all column names for a given trait, applying user overrides first.
resolve_cols <- function(headers, overrides, trait) {
  ov <- if (!is.null(overrides[[trait]])) overrides[[trait]] else character(0)
  get <- function(key, required = TRUE) {
    if (key %in% names(ov) && !is.na(ov[[key]])) return(ov[[key]])
    detect_col(headers, CANDIDATES[[key]], key, required)
  }
  # chrpos = combined "CHR:POS" column; if provided, chr/pos auto-detection is skipped
  has_chrpos <- "chrpos" %in% names(ov) && !is.na(ov[["chrpos"]])
  list(
    chrpos = if (has_chrpos) ov[["chrpos"]] else NULL,
    varid  = if ("varid" %in% names(ov) && !is.na(ov[["varid"]])) ov[["varid"]] else NULL,
    chr    = if (has_chrpos) NULL else get("chr"),
    pos    = if (has_chrpos) NULL else get("pos"),
    ea     = get("ea"),
    nea    = get("nea",    required = FALSE),
    beta   = get("beta"),
    se     = get("se"),
    p      = get("p",      required = FALSE),
    log10p = get("log10p", required = FALSE),
    n      = get("n",      required = FALSE),
    eaf    = get("eaf",    required = FALSE)
  )
}

# Lift over a data.table from hg38 to hg19. Expects CHR, POS columns.
liftover_to_hg19 <- function(dt, chain_file) {
  for (pkg in c("rtracklayer", "GenomicRanges", "IRanges", "S4Vectors")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("'%s' is required for liftover. Install via BiocManager::install('%s').", pkg, pkg))
  }
  chain  <- rtracklayer::import.chain(chain_file)
  chroms <- paste0("chr", gsub("^chr", "", dt$CHR, ignore.case = TRUE))
  gr     <- GenomicRanges::GRanges(
    seqnames = chroms,
    ranges   = IRanges::IRanges(start = as.integer(dt$POS), end = as.integer(dt$POS))
  )
  lifted <- rtracklayer::liftOver(gr, chain)
  keep   <- which(S4Vectors::elementNROWS(lifted) == 1)
  message(sprintf("    Liftover: %i / %i variants mapped to hg19", length(keep), nrow(dt)))
  if (length(keep) == 0) stop("Liftover returned zero variants.")
  dt_out       <- dt[keep]
  gr_hg19      <- unlist(lifted[keep])
  dt_out$CHR   <- gsub("^chr", "", as.character(GenomicRanges::seqnames(gr_hg19)))
  dt_out$POS   <- as.integer(GenomicRanges::start(gr_hg19))
  dt_out
}

# Load one chromosome's worth of data from the 2-column VAR_ID->rsID file and
# convert it to the 4-column format expected by align_to_rsid_map().
# Filters by grepping "^{chr}_" to avoid loading the full (potentially huge) file.
load_rsid_map_chr_from_file <- function(rsid_map_file, chr) {
  is_gz  <- endsWith(rsid_map_file, ".gz")
  prefix <- paste0("^", chr, "_")
  cmd <- if (is_gz) {
    sprintf("gzip -cd %s | grep -E '%s'", shQuote(rsid_map_file), prefix)
  } else {
    sprintf("grep -E '%s' %s", prefix, shQuote(rsid_map_file))
  }
  dt <- tryCatch(
    fread(cmd = cmd, header = FALSE,
          col.names = c("VAR_ID", "map_rsID"),
          data.table = TRUE, showProgress = FALSE),
    error = function(e) data.table()
  )
  if (nrow(dt) == 0L) return(data.table())

  # Extract REF/ALT from VAR_ID (format: CHR_POS_REF_ALT)
  parts <- tstrsplit(dt$VAR_ID, "_", fixed = TRUE)
  if (length(parts) < 4L) return(data.table())
  dt[, `:=`(
    hg19_posID  = paste0("chr", parts[[1]], ":", parts[[2]]),
    ref_allele  = parts[[3]],
    alt_allele  = parts[[4]]
  )]
  dt[, .(hg19_posID, map_rsID, ref_allele, alt_allele)]
}

# Align variants to a canonical REF/ALT reference using an rsID map.
# Input dt must have columns: CHR, POS, EA, NEA, BETA, SE, P_VALUE, N, [EAF], [rsID]
# Returns a new data.table with VAR_ID (CHR_POS_REF_ALT), SNP (CHR:POS),
# EA (= ALT), NEA (= REF), BETA (ALT-aligned), and original BETA/SE/P/N columns.
# Supply either rsid_map_dir (directory of per-chr files) or rsid_map_file (2-col file).
align_to_rsid_map <- function(dt, rsid_map_dir = NULL, rsid_map_file = NULL) {
  if (is.null(rsid_map_dir) && is.null(rsid_map_file))
    stop("Provide either rsid_map_dir or rsid_map_file.")
  result_list <- list()

  for (chr in sort(unique(dt$CHR))) {
    # Load the per-chromosome map from whichever source is available
    if (!is.null(rsid_map_dir)) {
      chr_file <- file.path(rsid_map_dir, paste0("chr", chr, ".txt"))
      if (!file.exists(chr_file)) {
        message(sprintf("    chr%s: no map file, skipping", chr))
        next
      }
      map <- fread(chr_file, header = FALSE,
                   col.names = c("hg19_posID", "map_rsID", "ref_allele", "alt_allele"),
                   key = "hg19_posID")
    } else {
      map <- load_rsid_map_chr_from_file(rsid_map_file, chr)
      if (nrow(map) == 0L) {
        message(sprintf("    chr%s: no entries in map file, skipping", chr))
        next
      }
      setkey(map, hg19_posID)
    }

    chr_dt <- dt[CHR == chr]
    chr_dt[, hg19_posID := paste0("chr", CHR, ":", POS)]
    setkey(chr_dt, hg19_posID)

    merged <- map[chr_dt, nomatch = 0]   # inner join: keep only map-matched variants
    if (nrow(merged) == 0) next

    # Case-insensitive allele matching (unordered pair)
    merged[, `:=`(
      ea_up  = toupper(EA),
      nea_up = toupper(NEA),
      ref_up = toupper(ref_allele),
      alt_up = toupper(alt_allele)
    )]
    merged[, alleles_match :=
             (ea_up == ref_up & nea_up == alt_up) |
             (ea_up == alt_up & nea_up == ref_up)]
    # need_flip = TRUE when EA is REF (BETA is REF-aligned; flip to make it ALT-aligned)
    merged[, need_flip := ea_up == ref_up & nea_up == alt_up]

    matched <- merged[alleles_match == TRUE]
    if (nrow(matched) == 0) next

    # Construct output columns; BETA and EAF are flipped when EA was REF
    matched[, `:=`(
      VAR_ID = paste(CHR, POS, ref_allele, alt_allele, sep = "_"),
      SNP    = paste0(CHR, ":", POS),
      rsID   = fifelse(!is.na(rsID) & rsID != "" & rsID != ".", rsID, map_rsID),
      BETA   = fifelse(need_flip, -BETA, BETA),
      EA     = alt_allele,
      NEA    = ref_allele,
      REF    = ref_allele,
      ALT    = alt_allele
    )]
    if ("EAF" %in% names(matched)) {
      matched[need_flip & !is.na(EAF), EAF := 1 - EAF]
    }

    keep_cols <- intersect(
      c("VAR_ID", "SNP", "rsID", "CHR", "POS", "EA", "NEA", "REF", "ALT",
        "BETA", "SE", "P_VALUE", "N", "EAF"),
      names(matched)
    )
    result_list[[chr]] <- matched[, ..keep_cols]
  }

  if (length(result_list) == 0) return(data.table())
  out <- rbindlist(result_list, fill = TRUE)
  # Remove duplicate VAR_IDs, keeping the most significant row (NA p-values last)
  setorder(out, P_VALUE, na.last = TRUE)
  unique(out, by = "VAR_ID")
}


# =============================================================================
# MAIN PER-FILE FUNCTION
# =============================================================================

standardize_gwas <- function(row, input_dir, output_dir, rsid_map_dir,
                              rsid_map_file, chain_file, col_overrides) {
  trait     <- row$trait
  study     <- row$study
  build     <- tolower(trimws(row$build))
  ext_N     <- row$N
  is_or     <- isTRUE(row$odds_ratio)
  out_gz    <- file.path(output_dir, sprintf("%s_%s_formatted.txt.gz", trait, study))

  if (file.exists(out_gz) && !OVERWRITE) {
    message(sprintf("[%s] Output already exists — skipping (set OVERWRITE <- TRUE to regenerate).", trait))
    return(invisible(NULL))
  }

  file_path <- file.path(input_dir, row$file)
  if (!file.exists(file_path)) stop(sprintf("Input file not found: %s", file_path))

  message(sprintf("\n[%s] Reading %s ...", trait, basename(file_path)))

  # ---- 1. Read ----
  dt <- fread(file_path, data.table = TRUE, stringsAsFactors = FALSE)
  message(sprintf("    %i rows x %i columns", nrow(dt), ncol(dt)))

  # ---- 2. Detect / apply column mappings ----
  message("    Detecting columns...")
  cols <- resolve_cols(names(dt), col_overrides, trait)

  # ---- 3. Rename to standard names ----
  # CHR + POS: from a combined "CHR:POS" column, or from separate columns
  if (!is.null(cols$chrpos)) {
    message(sprintf("    [override] chrpos -> '%s' (splitting on ':')", cols$chrpos))
    tmp <- tstrsplit(dt[[cols$chrpos]], ":", fixed = TRUE)
    dt[, CHR := tmp[[1L]]]
    dt[, POS := tmp[[2L]]]
  } else {
    setnames(dt, cols$chr, "CHR")
    setnames(dt, cols$pos, "POS")
  }
  setnames(dt, cols$ea,   "EA")
  setnames(dt, cols$beta, "BETA")
  setnames(dt, cols$se,   "SE")
  if (!is.null(cols$nea))    setnames(dt, cols$nea,    "NEA")
  if (!is.null(cols$eaf))    setnames(dt, cols$eaf,    "EAF")
  if (!is.null(cols$n))      setnames(dt, cols$n,      "N")
  # NEA from VAR_ID (CHR_POS_REF_ALT): used when no NEA column exists in the file
  if (!is.null(cols$varid) && !"NEA" %in% names(dt)) {
    message(sprintf("    [override] varid -> '%s' (deriving NEA from REF/ALT)", cols$varid))
    tmp <- tstrsplit(dt[[cols$varid]], "_", fixed = TRUE)
    dt[, NEA := fifelse(toupper(EA) == toupper(tmp[[4L]]), tmp[[3L]], tmp[[4L]])]
  }

  # P-value: prefer direct p-value; fall back to LOG10P conversion
  if (!is.null(cols$p)) {
    setnames(dt, cols$p, "P_VALUE")
    dt[, P_VALUE := as.numeric(P_VALUE)]
  } else if (!is.null(cols$log10p)) {
    setnames(dt, cols$log10p, "LOG10P_tmp")
    dt[, P_VALUE := 10^(-as.numeric(LOG10P_tmp))]
    dt[, LOG10P_tmp := NULL]
    message("    Converted LOG10P -> P_VALUE")
  } else {
    dt[, P_VALUE := NA_real_]
    message("    WARNING: no p-value column found; P_VALUE set to NA")
  }

  # ---- 4. Type coercions ----
  dt[, CHR   := as.character(gsub("^chr", "", CHR, ignore.case = TRUE))]
  dt[, POS   := as.integer(POS)]
  dt[, EA    := toupper(trimws(EA))]
  if ("NEA"  %in% names(dt)) dt[, NEA  := toupper(trimws(NEA))]
  dt[, BETA  := as.numeric(BETA)]
  dt[, SE    := as.numeric(SE)]

  # OR -> log(OR)
  if (is_or) {
    if (any(dt$BETA <= 0, na.rm = TRUE))
      warning("    Some OR values <= 0; log(OR) will produce NaN — check input.")
    dt[, BETA := log(BETA)]
    message("    Converted OR -> log(OR)")
  }

  # Sample size
  if (!is.na(ext_N)) {
    dt[, N := as.numeric(ext_N)]
  } else if ("N" %in% names(dt)) {
    dt[, N := as.numeric(N)]
  } else {
    dt[, N := NA_real_]
  }

  # rsID: carry through if present in input (will be used preferentially by align_to_rsid_map)
  if (!"rsID" %in% names(dt)) dt[, rsID := NA_character_]

  # ---- 5. Basic filters (each step logged) ----
  n_before <- nrow(dt)
  log_drop <- function(label, n_prev) {
    message(sprintf("    filter [%-28s] dropped %i", label, n_prev - nrow(dt)))
    nrow(dt)
  }
  n_prev <- nrow(dt)
  dt[CHR == "23", CHR := "X"]   # X coded as 23
  valid_chroms <- as.character(c(1:22, "X", "Y"))
  dt <- dt[CHR %in% valid_chroms];                       n_prev <- log_drop("chromosome not 1-22/X/Y", n_prev)
  dt <- dt[!is.na(POS) & POS > 0];                       n_prev <- log_drop("missing/invalid POS", n_prev)
  dt <- dt[!is.na(EA) & !is.na(BETA) & !is.na(SE) & SE > 0]
  n_prev <- log_drop("missing EA/BETA/SE or SE<=0", n_prev)
  # SNPs only (single-base, ACGT alleles)
  dt <- dt[nchar(EA) == 1L & grepl("^[ACGT]$", EA)];     n_prev <- log_drop("EA not single-base ACGT", n_prev)
  if ("NEA" %in% names(dt)) {
    dt <- dt[nchar(NEA) == 1L & grepl("^[ACGT]$", NEA)]; n_prev <- log_drop("NEA not single-base ACGT", n_prev)
    dt <- dt[EA != NEA];                                 n_prev <- log_drop("EA == NEA", n_prev)
  }

  message(sprintf("    %i / %i variants retained after basic filters", nrow(dt), n_before))
  if (nrow(dt) == 0) stop("No variants remain after filtering.")

  # Replace P_VALUE = 0 with a small non-zero value to avoid -Inf log
  dt[P_VALUE == 0, P_VALUE := .Machine$double.xmin]

  # ---- 6. Liftover (hg38 -> hg19) ----
  if (build == "hg38") {
    if (is.null(chain_file)) stop("build='hg38' but CHAIN_FILE is NULL.")
    dt <- liftover_to_hg19(dt, chain_file)
  } else if (build != "hg19") {
    warning(sprintf("    Unrecognised build '%s'; assuming hg19.", build))
  }

  # ---- 7a. WITH rsID map: canonical VAR_ID + ALT-aligned BETA ----
  use_map <- !is.null(rsid_map_dir) || !is.null(rsid_map_file)
  if (use_map) {
    if (!"NEA" %in% names(dt))
      stop("rsID map alignment requires an NEA column. Add 'nea' to COL_OVERRIDES or supply a file with a non-effect allele column.")
    n_before <- nrow(dt)
    n_sig_before <- sum(dt$P_VALUE < 5e-8, na.rm = TRUE)
    dt <- align_to_rsid_map(dt, rsid_map_dir = rsid_map_dir,
                                rsid_map_file = rsid_map_file)
    message(sprintf("    %i / %i variants matched to rsID map", nrow(dt), n_before))
    if (nrow(dt) > 0) {
      n_sig_after <- sum(dt$P_VALUE < 5e-8, na.rm = TRUE)
      message(sprintf("    Genome-wide significant (P<5e-8): %i before map, %i after (%i lost)",
                      n_sig_before, n_sig_after, n_sig_before - n_sig_after))
    }
    if (nrow(dt) == 0) stop("No variants matched the rsID map.")

  # ---- 7b. WITHOUT rsID map: chr:pos SNP ID, preserve input orientation ----
  } else {
    dt[, SNP := paste0(CHR, ":", POS)]
    setorder(dt, P_VALUE, na.last = TRUE)
    dt <- unique(dt, by = intersect(c("SNP", "EA", "NEA"), names(dt)))
    keep_cols <- intersect(
      c("SNP", "EA", "NEA", "BETA", "SE", "P_VALUE", "N", "EAF"),
      names(dt)
    )
    dt <- dt[, ..keep_cols]
    message("    No rsID map — output uses chr:pos SNP IDs (suitable for trait files).")
    if (!"NEA" %in% names(dt))
      message("    WARNING: no NEA column in output. Pipeline allele-matching will be limited.")
  }

  # ---- 8. Write output ----
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_txt <- sub("\\.gz$", "", out_gz)
  write_delim(as.data.frame(dt), out_txt, delim = "\t")
  system(paste("gzip", shQuote(out_txt)))
  message(sprintf("    -> %s  (%i variants)", basename(out_gz), nrow(dt)))
}


# =============================================================================
# RUN
# =============================================================================

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(MANIFEST))) {
  tryCatch(
    standardize_gwas(
      row           = MANIFEST[i, ],
      input_dir     = INPUT_DIR,
      output_dir    = OUTPUT_DIR,
      rsid_map_dir  = RSID_MAP_DIR,
      rsid_map_file = RSID_MAP_FILE,
      chain_file    = CHAIN_FILE,
      col_overrides = COL_OVERRIDES
    ),
    error = function(e) {
      message(sprintf("\n  ERROR [%s]: %s", MANIFEST$trait[i], e$message))
    }
  )
}

message("\nFormatting complete.")
