# =============================================================================
# post_bNMF_2025.R
# =============================================================================
# Post-hoc analysis functions applied after bNMF clustering.
#
# Functions:
#   calculate_cutoff()  — Determine optimal activity weight cutoff (elbow method)
#   v2g()               — Variant-to-gene mapping using Open Targets V2G scores
#                         (requires local V2G database files; see v2g_index_dir /
#                          v2g_scored_dir parameters)
#   do_post_analysis()  — Full post-hoc pipeline: load results, annotate nearest
#                         gene, export Excel workbook, liftover hg19→hg38,
#                         generate PRS score files, and optionally run V2G
#
# Dependencies:
#   GenomicRanges, Homo.sapiens, openxlsx, rtracklayer (liftover), vroom
# =============================================================================

library(tidyverse)
library(magrittr)
library(rtracklayer)
library(vroom)
library(openxlsx)
library(GenomicRanges)
library(Homo.sapiens)

# post-hoc functions

calculate_cutoff <- function(w, h) {
  
  w_mat <- data.frame(t(w)) %>%
    set_colnames(.["variant", ]) %>%
    subset(!rownames(.) %in% "variant")
  
  h_mat <- h
  
  # w_mat <- fread("./L2EU.H.mat.8.txt",stringsAsFactors = FALSE,data.table = F)
  
  dist_point_line <- function(a, slope, intercept) {
    b = c(1, intercept+slope)
    c = c(0, intercept)
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1,v2)
    return(abs(det(m))/sqrt(sum(v1*v1)))
  }
  
  # w_mat$cluster = NULL
  # names(w_mat) = NULL
  # weight = unlist(w_mat)
  w_list <- unlist(as.list(w_mat))
  h_list <- unlist(as.list(h_mat))
  
  weight = c(w_list, h_list)
  
  
  weights = as.data.frame(weight)
  weights$weight = as.numeric(as.character(weights$weight))
  weights = as.data.frame(weights[order(weights$weight,decreasing = T),])
  names(weights)[1] = "weight"
  
  numWeights <- dim(weights)[1]
  weights$x = c(1:numWeights)
  weights$x_scale = weights$x/numWeights
  weights$w_scale = weights$weight/max(weights$weight)
  
  cutoff_test <- quantile(weights$weight, prob=1-1/100)
  n = ifelse(sum(weights$weight > cutoff_test)>3, 1, 5)
  m = 80
  top1 = as.data.frame(weights[weights$weight > quantile(weights$weight, prob=1-n/100),])
  low80 = as.data.frame(weights[weights$weight <= quantile(weights$weight, prob=m/100),])
  weights$group = ifelse(weights$x %in% top1$x,paste0("top",n),
                         ifelse(weights$x %in% low80$x,paste0("bottom",m),"other"))
  
  top1line <- lm(w_scale ~ x_scale, data=top1)
  bottom80line <- lm(w_scale ~ x_scale, data=low80)
  
  weights$dist_top1 = rep(0,numWeights)
  for(i in 1:numWeights){
    x = weights$x_scale[i]
    y = weights$w_scale[i]
    weights$dist_top1[i] = dist_point_line(c(x,y),top1line$coefficients[2],top1line$coefficients[1])
  }
  
  weights$dist_low80 = rep(0,numWeights)
  for(i in 1:numWeights){
    x = weights$x_scale[i]
    y = weights$w_scale[i]
    weights$dist_low80[i] = dist_point_line(c(x,y),bottom80line$coefficients[2],bottom80line$coefficients[1])
  }
  
  weights$diff = weights$dist_top1 - weights$dist_low80
  cut = weights[weights$diff > 0,]
  cutoff = cut$weight[1]
  cutoff_num = cut$x[1]
  
  highlight <- weights %>% filter(x == cutoff_num)
  
  
  ggplot(weights, aes(x=x_scale, y=w_scale, color=group)) +
    geom_point() + 
    geom_abline(intercept = top1line$coefficients[1], slope = top1line$coefficients[2], color='blue') + 
    geom_abline(intercept = bottom80line$coefficients[1], slope = bottom80line$coefficients[2], color='red') +
    geom_point(data=highlight, aes(x=x_scale,y=w_scale), color='red',size=5)
  
  ggsave(file.path(main_dir,"bNMF_cutoff_calc.png"))
  
  return(cutoff)
}


v2g <- function(gwas_data, prefix, output_dir,
                v2g_index_dir  = NULL,
                v2g_scored_dir = NULL,
                gene_symbol_file = NULL) {
  # Map GWAS variants to genes using Open Targets V2G scores.
  #
  # Args:
  #   gwas_data        : Data frame with columns: chromosome, position, rsid
  #   prefix           : Output filename prefix
  #   output_dir       : Directory to write output file
  #   v2g_index_dir    : Directory containing per-chromosome V2G index files
  #                      (chr{N}.v2g_index.tsv.gz). Download from Open Targets.
  #   v2g_scored_dir   : Directory containing per-chromosome V2G scored files
  #                      (chr{N}.v2g_scored.tsv.gz). Download from Open Targets.
  #   gene_symbol_file : Path to Ensembl gene symbol table (grch38_df.tsv)
  #                      with columns: ensgene, symbol.
  #
  # Returns: writes <prefix>.v2g.tsv.gz to output_dir; returns all_combine invisibly.

  if (is.null(v2g_index_dir) || is.null(v2g_scored_dir) || is.null(gene_symbol_file)) {
    stop("v2g() requires v2g_index_dir, v2g_scored_dir, and gene_symbol_file to be specified.")
  }

  # Sanitize prefix for file naming
  prefix <- gsub("[^a-zA-Z0-9_-]", "_", prefix)

  # Read gene symbols once (not per-chromosome)
  gene_symbol <- vroom(gene_symbol_file, show_col_types = FALSE) %>%
    transmute(gene_id = ensgene, gene = symbol)

  chrs_present <- sort(unique(gwas_data$chromosome))
  all_combine <- data.frame()
  for (chr in chrs_present) {
    cat('Working on chr', chr, '\n')

    gwas_chr_data <- gwas_data %>% filter(chromosome == chr)

    # Read index; filter to only rsIDs in our variant set
    our_rsids <- gwas_chr_data$rsid
    index_file <- file.path(v2g_index_dir, sprintf("chr%d.v2g_index.tsv.gz", chr))
    index_data <- vroom(index_file, show_col_types = FALSE) %>%
      dplyr::rename(position_v2g = position, rsid = rs_id) %>%
      filter(rsid %in% our_rsids)

    # Use matched hg38 positions to filter the 97M-row scored file via shell
    our_positions <- unique(index_data$position_v2g)
    scored_file <- file.path(v2g_scored_dir, sprintf("chr%d.v2g_scored.tsv.gz", chr))
    if (length(our_positions) > 0) {
      pos_pattern <- paste(our_positions, collapse = "|")
      scored_data <- vroom(
        pipe(sprintf("gzip -dc '%s' | awk -F'\\t' 'NR==1 || ($2~/^(%s)$/)'", scored_file, pos_pattern)),
        show_col_types = FALSE
      ) %>%
        dplyr::transmute(chromosome = chr_id, position_hg38 = position, gene_id, overall_score, source_score_list)
    } else {
      scored_data <- data.frame(chromosome = integer(), position_hg38 = integer(),
                                gene_id = character(), overall_score = numeric(),
                                source_score_list = character())
    }

    merged_data <- left_join(gwas_chr_data,
                             index_data %>% dplyr::rename(chromosome = chr_id, position_hg38 = position_v2g),
                             by = c('chromosome', 'rsid'))
    merged_data2 <- left_join(merged_data, scored_data, by = c("chromosome", "position_hg38"))

    merged_v2g_data3 <- merged_data2 %>%
      group_by(rsid) %>%
      arrange(desc(overall_score)) %>%
      filter(!duplicated(rsid)) %>%
      ungroup()

    combine_df <- left_join(merged_v2g_data3, gene_symbol, by = 'gene_id')
    all_combine <- bind_rows(all_combine, combine_df)
  }
  
  # Ensure output directory exists
  # dir.create("output", showWarnings = FALSE, recursive = TRUE)
  
  # Write the combined data to a gzipped TSV file
  new_file = paste0(prefix, '.v2g.tsv.gz')
  write_tsv(all_combine, file = file.path(output_dir, new_file))
  print("Done with V2G!")
  return(invisible(all_combine))
}

do_post_analysis <- function(main_dir,
                            my_chain,
                            df_gwas,
                            k              = NULL,
                            cluster_names  = NULL,
                            loci_file      = "query",
                            make_excel     = TRUE,
                            do_liftover    = TRUE,
                            do_v2g         = FALSE,
                            existing_dir   = NULL,
                            v2g_index_dir  = NULL,
                            v2g_scored_dir = NULL,
                            gene_symbol_file = NULL) {
  # Full post-hoc pipeline: load bNMF outputs, annotate nearest genes,
  # create Excel workbook, perform hg19→hg38 liftover, write PRS score files,
  # and optionally run V2G mapping.
  #
  # Args:
  #   main_dir         : Output directory from summarize_bNMF()
  #   my_chain         : rtracklayer chain object for hg19→hg38 liftover
  #   df_gwas          : Aligned GWAS summary stats data frame (from Step 7)
  #   k                : Number of clusters (NULL = auto-select by majority vote)
  #   cluster_names    : Character vector of cluster labels (NULL = "X1","X2",...)
  #   loci_file        : "query" to annotate nearest gene; otherwise path to loci file
  #   make_excel       : Write per-cluster sorted weight Excel workbook
  #   do_liftover      : Perform hg19→hg38 coordinate liftover
  #   do_v2g           : Run Open Targets V2G mapping (requires v2g_* paths below)
  #   existing_dir     : Optional alternative directory for rsID_map.txt lookup
  #   v2g_index_dir    : Directory with per-chr V2G index files (for do_v2g=TRUE)
  #   v2g_scored_dir   : Directory with per-chr V2G scored files (for do_v2g=TRUE)
  #   gene_symbol_file : Path to Ensembl gene symbol table (for do_v2g=TRUE)
  
  # 1.) load data ---
  
  # find majority K from run summary
  df_run_summary <- fread(file.path(main_dir,"run_summary.txt"),
                          stringsAsFactors = F, data.table = F)
  
  k_counts <- df_run_summary %>%
    dplyr::count(K) %>%
    mutate(perc = n/nrow(df_run_summary))
  
  # ONLY calculate best K if user did not provide one
  if (is.null(k)) {
    k <- k_counts$K[which.max(k_counts$n)]
    message(sprintf("Automatically selected K=%d (majority vote)", k))
  } else {
    message(sprintf("Using manually requested K=%d", k))
  }
  
  if (is.null(cluster_names)) {
    cluster_names <- paste0("X",1:k)
  }
  
  w <- fread(sprintf("%s/L2EU.W.mat.%i.txt",main_dir, k),
             stringsAsFactors = FALSE, data.table = F) %>%
    dplyr::select(variant, all_of(cluster_names))
  n_variants <- nrow(w)
  
  h <- fread(sprintf("%s/L2EU.H.mat.%i.txt",main_dir,k),
             stringsAsFactors = FALSE, data.table = F) %>%
    rename_at(
      vars(starts_with("X2hr")), function(x) {gsub("X","",x) })
  
  # do cutoff
  cutoff <- calculate_cutoff(w, h)
  print(sprintf("Optimal weight cutoff is %.3f!", cutoff))

  # 1b.) Run V2G early so gene assignments are available for Excel output
  v2g_result <- NULL
  if (do_v2g) {
    v2g_input <- df_gwas %>%
      separate(ChrPos, into = c("chromosome", "position"), sep = ":", remove = TRUE) %>%
      mutate_at(vars(chromosome, position), as.integer) %>%
      dplyr::select(chromosome, position, rsid = rsID)
    v2g_result <- v2g(
      gwas_data        = v2g_input,
      prefix           = paste(basename(main_dir), "v2g", sep = "_"),
      output_dir       = main_dir,
      v2g_index_dir    = v2g_index_dir,
      v2g_scored_dir   = v2g_scored_dir,
      gene_symbol_file = gene_symbol_file
    )
  }

  # 2.) get nearest gene ---
  print("reading rsID map...")
  rsID_dir <- ifelse(is.null(existing_dir), main_dir, existing_dir)
  df_rsIDs <- fread(file.path(rsID_dir,"rsID_map.txt"),
                    data.table=F, stringsAsFactors=F)
  
  variants_have_prefix <- grepl("chr", w$variant[1])
  
  if (!'variant' %in% names(df_rsIDs)) {
    if (variants_have_prefix) {
      df_rsIDs <- df_rsIDs %>%
        separate(VAR_ID, into=c('CHR','POS','REF','ALT'),sep='_',remove = F) %>%
        mutate(variant = paste0('chr',CHR,":",POS))
    } else {
      df_rsIDs <- df_rsIDs %>%
        separate(VAR_ID, into=c('CHR','POS','REF','ALT'),sep='_',remove = F) %>%
        mutate(variant = paste0(CHR,":",POS))
    }
    
  }
  
  # IMPROVED VERSION - Handle variant filtering throughout the gene annotation pipeline
  
  print("Making excel file...")
  if (make_excel==T) {
    print("Making excel file...")
    
    if (loci_file=="query") {
      # Gene annotation logic with robust parsing
      geneRanges <- function(db, column="ENTREZID"){
        g <- genes(db, columns=column)
        col <- mcols(g)[[column]]
        genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
        mcols(genes)[[column]] <- as.character(unlist(col))
        genes
      }
      
      # Suppress the gene dropping messages
      suppressMessages({
        gns = geneRanges(Homo.sapiens, column="SYMBOL")
        df_entrez_data <- geneRanges(Homo.sapiens, column="ENTREZID")
      })
      
      df_gns <- data.frame(gns)
      gtf.gene <- df_gns %>%
        mutate(chr = gsub("chr","",seqnames)) %>%
        dplyr::select(chr, start, end, SYMBOL)
      
      df_entrez <- data.frame(df_entrez_data) %>%
        mutate(chr = gsub("chr","",seqnames)) %>%
        mutate(ChrPos = paste(chr,start,sep=":")) %>%
        dplyr::select(ChrPos, ENTREZID)
      
      range2GRanges <- function(df) {
        require(GenomicRanges)
        require(IRanges)
        gr <- GenomicRanges::GRanges(
          seqnames = df[,1],
          ranges=IRanges(start = df[,2], end = df[,3])
        )
        return(gr)
      }
      
      # ROBUST VARIANT PARSING
      snps <- w$variant
      message(sprintf("Annotating %d variants with nearest gene...", length(snps)))
      
      # SAFE PARSING - handles any format
      ranges_df <- data.frame(
        chr = character(length(snps)),
        start = numeric(length(snps)),
        end = numeric(length(snps)),
        original_variant = snps,  # Keep track of original variant names
        stringsAsFactors = FALSE
      )
      
      for (i in seq_along(snps)) {
        snp <- snps[i]
        
        tryCatch({
          if (grepl("^chr", snp)) {
            # Format: chr1:12345
            clean_snp <- gsub("^chr", "", snp)
            parts <- strsplit(clean_snp, ":")[[1]]
            ranges_df$chr[i] <- parts[1]
            ranges_df$start[i] <- as.numeric(parts[2])
            ranges_df$end[i] <- as.numeric(parts[2])
            
          } else if (grepl(":", snp)) {
            # Format: 1:12345
            parts <- strsplit(snp, ":")[[1]]
            ranges_df$chr[i] <- parts[1]
            ranges_df$start[i] <- as.numeric(parts[2])
            ranges_df$end[i] <- as.numeric(parts[2])
            
          } else if (grepl("_", snp)) {
            # Format: 1_12345_A_T
            parts <- strsplit(snp, "_")[[1]]
            ranges_df$chr[i] <- parts[1]
            ranges_df$start[i] <- as.numeric(parts[2])
            ranges_df$end[i] <- as.numeric(parts[2])
            
          } else {
            # Unknown format - skip
            ranges_df$chr[i] <- NA
            ranges_df$start[i] <- NA
            ranges_df$end[i] <- NA
          }
        }, error = function(e) {
          # Skip any problematic variants
          ranges_df$chr[i] <- NA
          ranges_df$start[i] <- NA
          ranges_df$end[i] <- NA
        })
      }
      
      # Remove failed parses
      valid_rows <- !is.na(ranges_df$chr) & !is.na(ranges_df$start) & !is.na(ranges_df$end)
      ranges_df_clean <- ranges_df[valid_rows, ]
      snps_clean <- ranges_df_clean$original_variant  # Use original variant names
      
      message(sprintf("  Parsed %d/%d variants successfully", nrow(ranges_df_clean), length(snps)))
      
      if (nrow(ranges_df_clean) == 0) {
        warning("No variants could be parsed for gene annotation - using variant IDs as gene names.")
        w_wLoci <- w %>%
          mutate(gene = variant)
      } else {
        
        # Create GRanges object with error checking
        snps.granges <- GenomicRanges::GRanges(
          seqnames = ranges_df_clean$chr,
          ranges = IRanges::IRanges(start = ranges_df_clean$start, end = ranges_df_clean$end)
        )
        names(snps.granges) <- snps_clean
        
        message(sprintf("  Created GRanges with %d variants", length(snps.granges)))
        
        # Convert genes to GRanges
        gtf.granges <- GenomicRanges::GRanges(
          seqnames = gtf.gene$chr,
          ranges = IRanges::IRanges(start = as.numeric(gtf.gene$start), end = as.numeric(gtf.gene$end))
        )
        names(gtf.granges) <- gtf.gene$SYMBOL
        
        # Find nearest genes
        hits <- GenomicRanges::nearest(snps.granges, gtf.granges)
        
        # Prepare gene annotation data
        df_gns_final <- df_gns %>%
          mutate(chr = gsub("chr","",seqnames)) %>%
          mutate(ChrPos = paste(chr,start,sep=":")) %>%
          dplyr::select(chr, start, end, ChrPos, SYMBOL) %>%
          merge(df_entrez, by="ChrPos")
        
        # Create hits data frame
        df_hits <- data.frame(
          gene = names(gtf.granges)[hits],
          variant = names(snps.granges),
          stringsAsFactors = FALSE
        ) %>%
          merge(df_gns_final[,c("SYMBOL","ENTREZID")], by.x="gene", by.y="SYMBOL", all.x=TRUE) %>%
          filter(!duplicated(variant))
        
        message(sprintf("  Gene annotation completed for %d variants", nrow(df_hits)))
        
        # IMPORTANT: Handle variants that couldn't be parsed
        # Include all original variants, with gene annotation where available
        w_wLoci <- w %>%
          left_join(df_hits, by="variant") %>%
          mutate(gene = dplyr::coalesce(gene, variant))  # Use variant name if no gene found
        
        message(sprintf("  Final annotation: %d variants", nrow(w_wLoci)))
      }
      
    } else {
      w_wLoci <- w %>%
        mutate(gene=variant)
    }
    
    if (exists("df_rsIDs") && nrow(df_rsIDs) > 0) {
      w_excel_ready <- w_wLoci %>%
        inner_join(df_rsIDs, by = "variant")
      n_missing <- nrow(w_wLoci) - nrow(w_excel_ready)
      if (n_missing > 0) {
        warning(sprintf("%d variants could not be matched to rsID map and will be dropped from Excel output.", n_missing))
      }
    } else {
      message("No rsID mapping available - using variant names as identifiers.")
      w_excel_ready <- w_wLoci %>%
        mutate(VAR_ID = variant, rsID = variant)
    }

    # Rename nearest gene column; join OpenTargets V2G gene if available
    w_excel_ready <- w_excel_ready %>%
      dplyr::rename(Nearest_Gene = gene)
    if (!is.null(v2g_result)) {
      v2g_genes <- v2g_result %>%
        dplyr::select(rsID = rsid, OpenTargets_Gene = gene) %>%
        dplyr::filter(!is.na(OpenTargets_Gene))
      w_excel_ready <- w_excel_ready %>%
        left_join(v2g_genes, by = "rsID")
    }

    # Make Excel with improved error handling
    tryCatch({
      h.exc <- t(h) %>%
        data.frame() %>%
        mutate(trait = row.names(.)) %>%
        relocate(trait)
      
      w.exc <- w_excel_ready %>%
        dplyr::select(-any_of(c('ENTREZID', 'variant'))) %>%
        relocate(any_of(c("VAR_ID", "rsID", "Nearest_Gene", "OpenTargets_Gene")))

      # Determine last ID column safely
      id_cols <- intersect(c("VAR_ID", "rsID", "Nearest_Gene", "OpenTargets_Gene"), colnames(w.exc))
      last_id_col <- length(id_cols)
      
      posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      
      wb <- openxlsx::createWorkbook()
      for (i in 1:length(cluster_names)){
        w_temp <- w.exc %>%
          dplyr::arrange(desc(!!as.symbol(cluster_names[i])))
        h_temp <- h.exc %>%
          arrange(desc(!!as.symbol(cluster_names[i])))
        
        cur_sheet = sprintf("cluster%i",i)
        addWorksheet(wb, cur_sheet)
        start_col1 <- 1
        start_col2 <- ncol(w_temp)+3
        
        writeData(wb=wb, sheet=cur_sheet, x=w_temp, startCol = start_col1, startRow = 1,
                  rowNames=F,colNames = T)
        writeData(wb=wb, sheet=cur_sheet, x=h_temp, startCol = start_col2, startRow = 1,
                  rowNames=F,colNames = T)
        
        # Safe conditional formatting
        if (last_id_col > 0 && last_id_col + i <= ncol(w_temp)) {
          conditionalFormatting(wb, cur_sheet,
                                cols = last_id_col + i, rows = 2:(nrow(w_temp)+1),
                                rule = sprintf("%s2>=%.7f",
                                               openxlsx::int2col(last_id_col + i), cutoff),
                                style = posStyle)
        }
        
        if (start_col2 + i <= ncol(h_temp) + start_col2) {
          conditionalFormatting(wb, cur_sheet,
                                cols = start_col2 + i, rows = 2:(nrow(h_temp)+1),
                                rule = sprintf("%s2>=%.7f",
                                               openxlsx::int2col(start_col2 + i), cutoff),
                                style = posStyle)
        }
      }
      
      fp <- file.path(main_dir, sprintf("sorted_cluster_weights_K%i_rev.xlsx",k))
      openxlsx::saveWorkbook(wb, fp, overwrite = T)
      
      cat(sprintf("Excel file saved: %s\n", fp))
      
    }, error = function(e) {
      cat("Error creating Excel file:", e$message, "\n")
      cat("Continuing without Excel output...\n")
    })
  }  
  
  # 3.) hg19-to-hg38 liftover --
  
  if (do_liftover==T) {
    print("Doing liftover...")
    if (is.character(my_chain)) my_chain <- rtracklayer::import.chain(my_chain)
    
    df_weights <- w_wLoci %>%
      inner_join(df_rsIDs, by=c('variant')) %>%
      # dplyr::select(VAR_ID, rsID, gene) %>%
      separate(VAR_ID, into=c("CHR", "POS_hg19", "REF", "ALT"), sep = "_", remove = F) %>%
      mutate(ChrPos = paste(CHR, POS_hg19, sep = ":")) %>%
      mutate(group = seq.int(nrow(.))) %>%
      mutate_at(.vars = c('CHR','POS_hg19'),as.integer)
    
    gr <- GRanges(seqnames = paste0("chr",df_weights$CHR),
                  strand = "*",
                  ranges = IRanges(start = df_weights$POS_hg19,
                                   width = 1))
    
    df_liftover <- as.data.frame(rtracklayer::liftOver(x = gr, chain = my_chain)) %>%
      dplyr::select(group, seqnames, Pos_hg38=start) %>%
      mutate(ChrPos_hg38=paste(seqnames, Pos_hg38, sep=":"))
    
    missing <- df_weights %>%
      filter(!group %in% df_liftover$group)
    if (nrow(missing)>0) {
      print("Variants not found in liftover...")
      print(missing$variant)
    }
    
    # merge weights with hg38 SNPs using rsID
    #   new file should have VAR_ID_hg38, rsID, Effect_Allele, BETA and the cluster columns (X1, X2, etc.)
    input_snps38 <- df_weights %>%
      inner_join(df_liftover[, c("group","Pos_hg38")], by="group") %>%
      mutate(VAR_ID_hg38 = paste(CHR, Pos_hg38, REF, ALT, sep = "_")) %>%
      inner_join(df_gwas[,c("ChrPos","Risk_Allele","BETA")], by="ChrPos")
    
    
    liftover_map <- input_snps38 %>%
      dplyr::select(VAR_ID_hg19=VAR_ID, Risk_Allele, rsID, VAR_ID_hg38)
    print("Writing liftover map to directory...")
    write_delim(liftover_map, file.path(main_dir,"liftover_map.txt"),
                delim = "\t",col_names = T)
    
    
    # make weight files
    input_snps38.1 <- input_snps38 %>%
      dplyr::select(VAR_ID_hg38, rsID, Risk_Allele, BETA, starts_with("X"))
    
    
    #---- CREATE WEIGHTS FILE FOR EACH CLUSTER-----
    
    
    cluster_weights <- input_snps38.1 %>%
      dplyr::select(VAR_ID_hg38, starts_with("X")) %>%
      column_to_rownames(var="VAR_ID_hg38")
    
    
    cluster_info = data.frame(weight=numeric(),VAR_ID_hg38=character(),cluster=character())
    
    for(i in 1:ncol(cluster_weights)){
      w_tmp = cluster_weights[i]
      names(w_tmp) = c("Weight")
      w_tmp$VAR_ID_hg38 = rownames(w_tmp)
      w_tmp = w_tmp[w_tmp$Weight>=cutoff,]
      w_tmp$cluster=cluster_names[i]
      cluster_info = rbind(cluster_info,w_tmp)
    }
    
    
    # for PRS, need following columns:
    # Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight (Effect_Allele=Risk_Allele here)
    #   Ex: 11, 92673828, C, T, rs1387153, T, 4.658804585
    input_df <- input_snps38 %>%
      merge(cluster_info, by="VAR_ID_hg38") %>%
      separate(VAR_ID_hg38, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
      arrange(as.integer(Chr), Pos) %>%
      dplyr::select(Chr, Pos, Ref, Alt, RSID=rsID, Risk_Allele, Weight, cluster)
    # check no issues with mapping
    sum(is.na(input_df))
    
    # create file with snps and weights for each cluster
    print("Writing score files...")
    dir.create(file.path(main_dir,"score_info"))
    for (c in cluster_names) {
      df <- input_df %>%
        subset(cluster==c) %>%
        dplyr::select(-c(cluster))
      print(nrow(df))
      write.csv(x=df,
                file=file.path(main_dir,"score_info",sprintf("%s.csv", c)),
                row.names=F, quote=F)
    }
    
    # repeat for all SNPs (for total PRS)
    all_snps <- input_snps38 %>%
      separate(VAR_ID_hg38, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
      arrange(as.integer(Chr), Pos) %>%
      mutate(Weight=abs(as.numeric(BETA))) %>%
      dplyr::select(Chr, Pos, Ref, Alt, RSID=rsID, Risk_Allele, Weight)
    write.csv(all_snps,
              file=file.path(main_dir,"score_info","Total_GRS.csv"),
              row.names=F, quote=F)
    
    print("Writing BCFtools input file...")
    for_VCF <- input_snps38 %>%
      # filter(VAR_ID_hg38 %in% cluster_info$VAR_ID_hg38) %>%
      separate(VAR_ID_hg38, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
      arrange(as.integer(Chr), as.integer(Pos)) %>%
      mutate(Chr = paste0("chr",Chr), start = Pos, end = Pos) %>%
      dplyr::select(Chr, start, end)
    write_delim(for_VCF, file.path(main_dir,"bcftools_input.txt"),
                delim = "\t",col_names = F)
  }
  
  
  return(list(df_run_summary=df_run_summary, k=k, cluster_names=cluster_names, w=w_wLoci, h=h, cutoff=cutoff))
}
