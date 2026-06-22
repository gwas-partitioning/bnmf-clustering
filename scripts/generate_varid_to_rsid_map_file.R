# ==============================================================================
# Generate rsID to VAR_ID Mapping File for bNMF Clustering Pipeline
# ==============================================================================
# Downloads 1000 Genomes Phase 3 VCFs from Ensembl (GRCh37) and builds a
# VAR_ID <-> rsID mapping file for the bNMF pipeline.
#
# Output format: hg19_posID (chrN:POS) \t rsID \t ref_allele \t alt_allele (no header)
# Reference genome: GRCh37/hg19
#
# Requirements: curl, data.table
# Disk space: up to ~20 GB during download (files deleted after processing)
# ==============================================================================

library(curl)
library(data.table)

# Configuration
CURRENT    <- getwd()
base_url   <- "http://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/"
output_file <- file.path(CURRENT, "list_VARID_rsID_updated.txt")
temp_dir   <- file.path(CURRENT, "temp_vcf")

if (!dir.exists(temp_dir)) dir.create(temp_dir)

chromosomes  <- 1:22
filenames    <- sprintf("homo_sapiens-chr%s.vcf.gz", chromosomes)
urls         <- paste0(base_url, filenames)
local_files  <- file.path(temp_dir, filenames)

# ==============================================================================
# Step 1: Parallel download of all chromosomes
# ==============================================================================
to_download <- which(!file.exists(local_files))

if (length(to_download) > 0) {
  message(sprintf("Downloading %d chromosome VCFs in parallel...", length(to_download)))
  results <- multi_download(
    urls[to_download],
    local_files[to_download],
    progress  = TRUE,
    resume    = TRUE,
    timeout   = 7200
  )

  failed <- which(results$status_code != 200 | nchar(results$error) > 0)
  if (length(failed) > 0) {
    message("The following chromosomes failed to download:")
    print(results[failed, c("url", "status_code", "error")])
    for (f in results$destfile[failed]) {
      if (file.exists(f)) file.remove(f)
    }
  }
  message(sprintf("Downloads complete: %d/%d succeeded.",
                  length(to_download) - length(failed), length(to_download)))
} else {
  message("All chromosome files already present, skipping download.")
}

# ==============================================================================
# Step 2: Process each chromosome (fread via shell pipe — skips genotype columns)
# ==============================================================================
process_chromosome <- function(chr, local_file) {
  message(sprintf("Processing chromosome %s...", chr))

  if (!file.exists(local_file) || file.size(local_file) == 0) {
    message(sprintf("  Skipping chr%s — file not available.", chr))
    return(NULL)
  }

  tryCatch({
    vars <- fread(
      cmd          = sprintf("gzip -dc '%s' | grep -v '^#'", local_file),
      select       = 1:5,
      col.names    = c("CHROM", "POS", "ID", "REF", "ALT"),
      header       = FALSE,
      showProgress = FALSE
    )

    if (nrow(vars) == 0) {
      message(sprintf("  No variants in chr%s.", chr))
      return(NULL)
    }

    vars_filtered <- vars[
      ID != "." & grepl("^rs[0-9]+$", ID) &
      nchar(REF) == 1L & nchar(ALT) == 1L &
      REF %chin% c("A","C","G","T") &
      ALT %chin% c("A","C","G","T") &
      !grepl(",", ALT)
    ]

    vars_filtered[, hg19_posID := paste0("chr", CHROM, ":", POS)]
    vars_filtered <- vars_filtered[, .(hg19_posID, rsID = ID, ref_allele = REF, alt_allele = ALT)]
    vars_filtered <- unique(vars_filtered, by = "hg19_posID")
    vars_filtered <- unique(vars_filtered, by = "rsID")

    message(sprintf("  Chr%s: %s variants after filtering.",
                    chr, format(nrow(vars_filtered), big.mark = ",")))

    file.remove(local_file)
    return(vars_filtered)

  }, error = function(e) {
    message(sprintf("  Error processing chr%s: %s", chr, e$message))
    if (file.exists(local_file)) file.remove(local_file)
    return(NULL)
  })
}

all_variants <- vector("list", length(chromosomes))
for (i in seq_along(chromosomes)) {
  all_variants[[i]] <- process_chromosome(chromosomes[i], local_files[i])
}

# ==============================================================================
# Step 3: Combine and write output
# ==============================================================================
all_variants <- Filter(Negate(is.null), all_variants)

if (length(all_variants) > 0) {
  message("\nCombining results from all chromosomes...")
  vars_final <- rbindlist(all_variants)
  vars_final <- unique(vars_final, by = "hg19_posID")
  vars_final <- unique(vars_final, by = "rsID")
  setorder(vars_final, hg19_posID)

  message(sprintf("Total variants: %s", format(nrow(vars_final), big.mark = ",")))
  message(sprintf("Writing to: %s", output_file))

  fwrite(vars_final, file = output_file, sep = "\t", col.names = FALSE, quote = FALSE)

  message("\n=== SUMMARY ===")
  chr_counts <- vars_final[, .(count = .N), by = .(CHR = sub("chr(.*):.*", "\\1", hg19_posID))]
  chr_counts[, CHR := as.integer(CHR)]
  setorder(chr_counts, CHR)
  for (i in seq_len(nrow(chr_counts))) {
    message(sprintf("  Chr%s: %s variants",
                    chr_counts$CHR[i], format(chr_counts$count[i], big.mark = ",")))
  }
  message(sprintf("\nFile format: VAR_ID<tab>rsID (no header)\nSaved: %s", output_file))

} else {
  message("\nERROR: No variants successfully processed.")
  message("Check your internet connection and the Ensembl FTP site.")
}

# Cleanup
if (dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
  message("Temporary files cleaned up.")
}

message("\nScript completed!")
