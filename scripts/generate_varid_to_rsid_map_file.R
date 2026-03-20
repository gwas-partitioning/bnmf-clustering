# ==============================================================================
# Generate rsID to VAR_ID Mapping File for bNMF Clustering Pipeline
# ==============================================================================
# This script downloads 1000 Genomes Phase 3 data from Ensembl (by chromosome)
# and creates a mapping file compatible with the bNMF clustering pipeline
# 
# Output format: VAR_ID (CHR_POS_REF_ALT) \t rsID
# Reference genome: GRCh37/hg19
# ==============================================================================

library(vcfR)
library(dplyr)

# Configuration
CURRENT <- getwd()
base_url <- "http://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/"
output_file <- file.path(CURRENT, "list_VARID_rsID_updated.txt")
temp_dir <- file.path(CURRENT, "temp_vcf")

# Create temporary directory
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

# Define chromosomes to process
chromosomes <- c(1:22)

# Function to download and process a single chromosome
process_chromosome <- function(chr) {
  message(sprintf("Processing chromosome %s...", chr))
  
  # Construct URL and local filename
  filename <- sprintf("homo_sapiens-chr%s.vcf.gz", chr)
  url <- paste0(base_url, filename)
  local_file <- file.path(temp_dir, filename)
  
  # Download if not exists
  if (!file.exists(local_file)) {
    message(sprintf("  Downloading %s...", filename))
    tryCatch({
      download.file(url, local_file, mode = "wb", quiet = TRUE)
    }, error = function(e) {
      message(sprintf("  Error downloading chr%s: %s", chr, e$message))
      return(NULL)
    })
  }
  
  # Check if file exists and has content
  if (!file.exists(local_file) || file.size(local_file) == 0) {
    message(sprintf("  Skipping chr%s - file not available", chr))
    return(NULL)
  }
  
  # Read and process VCF
  tryCatch({
    message(sprintf("  Reading VCF for chr%s...", chr))
    vcf <- read.vcfR(local_file, verbose = FALSE)
    
    if (nrow(vcf@fix) == 0) {
      message(sprintf("  No variants found in chr%s", chr))
      return(NULL)
    }
    
    # Extract variant information
    vars <- as.data.frame(vcf@fix[, c("CHROM", "POS", "ID", "REF", "ALT")])
    
    # Apply quality control filters
    vars_filtered <- vars %>%
      # Valid rsIDs only
      filter(ID != "." & grepl("^rs[0-9]+$", ID)) %>%
      # Biallelic SNPs only
      filter(
        nchar(REF) == 1 & nchar(ALT) == 1 &
          REF %in% c("A", "C", "G", "T") & ALT %in% c("A", "C", "G", "T") &
          !grepl(",", ALT)
      ) %>%
      # Create VAR_ID
      mutate(VAR_ID = paste(CHROM, POS, REF, ALT, sep = "_")) %>%
      select(VAR_ID, rsID = ID) %>%
      # Remove duplicates within chromosome
      distinct(VAR_ID, .keep_all = TRUE) %>%
      distinct(rsID, .keep_all = TRUE)
    
    message(sprintf("  Chr%s: %s variants after filtering", 
                    chr, format(nrow(vars_filtered), big.mark = ",")))
    
    # Clean up chromosome file
    file.remove(local_file)
    
    return(vars_filtered)
    
  }, error = function(e) {
    message(sprintf("  Error processing chr%s: %s", chr, e$message))
    if (file.exists(local_file)) file.remove(local_file)
    return(NULL)
  })
}

# Main processing
message("Starting chromosome-by-chromosome processing...")
message(sprintf("Will process chromosomes: %s", paste(chromosomes, collapse = ", ")))

# Process all chromosomes
all_variants <- list()
total_downloaded <- 0

for (chr in chromosomes) {
  result <- process_chromosome(chr)
  if (!is.null(result) && nrow(result) > 0) {
    all_variants[[as.character(chr)]] <- result
    total_downloaded <- total_downloaded + 1
  }
}

# Combine all chromosomes
if (length(all_variants) > 0) {
  message("\nCombining results from all chromosomes...")
  vars_final <- bind_rows(all_variants)
  
  # Final deduplication across chromosomes (in case of overlaps)
  vars_final <- vars_final %>%
    distinct(VAR_ID, .keep_all = TRUE) %>%
    distinct(rsID, .keep_all = TRUE) %>%
    arrange(VAR_ID)
  
  message(sprintf("Combined total: %s variants", format(nrow(vars_final), big.mark = ",")))
  
  # Write output file
  message(sprintf("Writing mapping file to: %s", output_file))
  write.table(vars_final, 
              file = output_file, 
              append = FALSE, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  # Summary statistics
  message("\n=== SUMMARY ===")
  message(sprintf("Chromosomes successfully processed: %d/%d", total_downloaded, length(chromosomes)))
  message(sprintf("Total variants in final mapping: %s", format(nrow(vars_final), big.mark = ",")))
  
  # Chromosome breakdown
  chr_summary <- vars_final %>%
    separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
    count(CHR, name = "count") %>%
    arrange(as.numeric(CHR)) 
  
  message("Variants by chromosome:")
  for(i in 1:nrow(chr_summary)) {
    message(sprintf("  Chr %s: %s variants", 
                    chr_summary$CHR[i], 
                    format(chr_summary$count[i], big.mark = ",")))
  }
  
  message(sprintf("\nMapping file saved as: %s", basename(output_file)))
  message("File format: VAR_ID<tab>rsID (no header)")
  
} else {
  message("\nERROR: No variants were successfully processed!")
  message("Please check your internet connection and the Ensembl FTP site.")
}

# Clean up temporary directory
if (dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
  message("Temporary files cleaned up.")
}

# Additional notes for users
message("\n=== NOTES ===")
message("• This mapping is based on 1000 Genomes Phase 3 data (GRCh37/hg19)")
message("• Only biallelic SNPs with valid rsIDs are included")
message("• Multi-allelic variants and indels are excluded")
message("• If some chromosomes failed, you can re-run to retry those only")

message("\nScript completed!")