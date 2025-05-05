library(tidyverse)
library(furrr)
library(progressr)
library(rtracklayer)
library(vroom)
##########################################################################
# Copyright (c) 2017, Broad Institute
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#     Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#     Neither the name of the Broad Institute nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#########################################################################

######################################################################
# Bayesian NMF algorithms for clustering
######################################################################
# For implementation details see the ppaer 
# Udler MS, Kim J, von Grotthuss M,
# Bonàs-Guarch S, Cole JB, Chiou J, et al. (2018)
# Type 2 diabetes genetic loci informed by multi-trait
# associations point to disease mechanisms and
# subtypes: A soft clustering analysis. PLoS Med 15
# (9): e1002654.
###########################
# For details on the original algorithms 
# see Tan, V.Y. & Févotte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence.
# IEEE Trans. Pattern Anal. Mach. Intell. 35, 1592–1605 (2013).
######################################################################

BayesNMF.L2EU <- function(
  V0, n.iter=10000, a0=10, tol=1e-7, K=15, K0=15, phi=1.0 #20, 10
) {
  
  # Bayesian NMF with half-normal priors for W and H
  # V0: input z-score matrix (variants x traits)
  # n.iter: Number of iterations for parameter optimization
  # a0: Hyper-parameter for inverse gamma prior on ARD relevance weights
  # tol: Tolerance for convergence of fitting procedure
  # K: Number of clusters to be initialized (algorithm may drive some to zero)
  # K0: Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
  # phi: Scaling parameter
  
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0)
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  W <- matrix(runif(N * K) * Vmax, ncol=K)
  H <- matrix(runif(M * K) * Vmax, ncol=M)
  
  I <- array(1, dim=c(N, M))
  V.ap <- W %*% H + eps
  
  phi <- sd(V)^2 * phi
  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda.bound <- b0 / C
  lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
  lambda.cut <- lambda.bound * 1.5
  
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  count <- 1
  while (del >= tol & iter < n.iter) {
    H <- H * (t(W) %*% V) / 
      (t(W) %*% V.ap + phi * H * matrix(rep(1 / lambda, M), ncol=M) + eps)
    V.ap <- W %*% H + eps
    W <- W * (V %*% t(H)) / 
      (V.ap %*% t(H) + phi * W * t(matrix(rep(1 / lambda, N), ncol=N)) + eps)
    V.ap <- W %*% H + eps
    lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    like <- sum((V - V.ap)^2) / 2
    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi * sum((0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / 
                                         lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V - V.ap)^2)
    if (iter %% 100 == 0) {
      cat(iter, n.evid[[iter]], n.like[[iter]], n.error[[iter]], del, 
          sum(colSums(W) != 0), sum(lambda >= lambda.cut), '\n')
    }
    iter <- iter + 1
  }
  return(list(
    W,  # Variant weight matrix (N x K)
    H,  # Trait weight matrix (K x M)
    n.like,  # List of reconstruction errors (sum of squared errors / 2) per iteration
    n.evid,  # List of negative log-likelihoods per iteration
    n.lambda,  # List of lambda vectors (shared weights for each of K clusters, some ~0) per iteration
    n.error  # List of reconstruction errors (sum of squared errors) per iteration
  ))
}


run_bNMF <- function(z_mat, n_reps=10, random_seed=1, K=20, K0=10, tolerance=1e-7) {
  
  # Given an input matrix as created by prep_z_matrix(), run the bNMF procedure
  # a series of times to generate results and evaluate cluster stability
  
  print(paste0("Running bNMF clustering procedure (", n_reps, " iterations)..."))
  print(sprintf("Using tolerance of %.2e!",tolerance))
  
  set.seed(random_seed)
  
  bnmf_reps <- lapply(1:n_reps, function(r) {
    print(paste("ITERATION",r))
    res <- BayesNMF.L2EU(V0 = z_mat, K=K, K0=K0, tol=tolerance)
    names(res) <- c("W", "H", "n.like", "n.evid", "n.lambda", "n.error")
    res
  })
  bnmf_reps
}


# Activate a global handler for progress bars

run_bNMF_parallel <- function(z_mat, n_reps = 10, random_seed = 1, K = 20, K0 = 10, tolerance = 1e-7) {
  
  print(paste0("Running bNMF clustering procedure in parallel! (", n_reps, " iterations)..."))
  print(sprintf("Using tolerance of %.2e!", tolerance))
  
  # Run each repetition in parallel
  bnmf_reps <- future_map(1:n_reps, function(rep) {
    
    # Logging progress
    log_message <- paste("ITERATION", rep)
    cat(log_message, "\n")  # Print to console to verify
    
    # Run the Bayesian NMF function and store the result
    res <- BayesNMF.L2EU(V0 = z_mat, K = K, K0 = K0, tol = tolerance)
    names(res) <- c("W", "H", "n.like", "n.evid", "n.lambda", "n.error")
    
    return(res)  # Explicitly return `res` from the function
  },
  .options = furrr_options(seed = random_seed)  # Set the seed to a specific number
  )
  
  return(bnmf_reps)
}

summarize_bNMF <- function(bnmf_reps, dir_save=NULL) {
  
  # Given output from bNMF (list of length N_iterations),
  # generate summary tables and plots

  make_run_summary <- function(reps) {
    
    # Given a list of bNMF iteration outputs, summarize the K choices and associated likelihoods across runs
    
    run_summary <- map_dfr(1:length(reps), function(i) {
      res <- reps[[i]]
      final_lambdas <- res$n.lambda[[length(res$n.lambda)]]
      tibble(
        run=i,
        K=sum(final_lambdas > min(final_lambdas)),  # Assume that lambdas equal to the minimum lambda are ~ 0
        evid=res$n.evid[[length(res$n.evid)]]  # Evidence = -log_likelihood
      )
    }) %>%
      arrange(evid)
    
    unique.K <- table(run_summary$K)
    n.K <- length(unique.K)  # Number of distinct K
    MAP.K.run <- sapply(names(unique.K), function(k) {  # bNMF run index with the maximum posterior for given K
      tmp <- run_summary[run_summary$K == k, ]
      tmp$run[which.min(tmp$evid)]
    })
    
    list(run_tbl=run_summary, unique.K=unique.K, MAP.K.run=MAP.K.run)
  }
  if (!is.null(dir_save)) {
    dir.create(file.path(dir_save))
    dir_save=paste0(dir_save,"/")
  } else {dir_save="./"}
  
  print("Summarizing bNMF results...")
  
  print("Writing table of chosen K across iterations...")
  run_summary <- make_run_summary(bnmf_reps)
  write_tsv(run_summary$run_tbl, paste0(dir_save,"run_summary.txt"))

  n.K <- length(run_summary$unique.K)  # Number of distinct K
  
  get_W <- function(clustering) {
    W_raw <- clustering$W
    W_raw[, colSums(W_raw > 1e-10) > 0]
  }
  
  get_H <- function(clustering) {
    H_raw <- clustering$H
    H_raw[rowSums(H_raw > 1e-10) > 0, ]
  }
  
  print("Plotting variant and trait contributions...")
  silent <- sapply(names(run_summary$unique.K), function(k) {  # Create heatmaps for MAP iteration for each K
    res <- bnmf_reps[[run_summary$MAP.K.run[as.character(k)]]]
    W <- res$W[, colSums(res$W) != 0]  # feature-cluster association matrix
    H <- res$H[rowSums(res$H) != 0, ]  # cluster-gene association matrix
    W[W < 1.e-10] <- 0
    H[H < 1.e-10] <- 0
    
    W0 <- data.frame(W)
    W0[, "variant"] <- rownames(W)
    H0 <- data.frame(H)
    H0[, "cluster"] <- rownames(H)
    
    write_tsv(W0, file=paste0(dir_save,"L2EU.W.mat.", k, ".txt"))
    write_tsv(H0, file=paste0(dir_save,"L2EU.H.mat.", k, ".txt"))
    
    mat.reconstructed <- W %*% H   # reconstructed matrix == approximation for the input matrix 
    
    # Setup for plotting
    scale0 <- 0.8
    scale <- 1
    g.ordering <- paste("G", seq(1:ncol(W)), sep="")
    color.axis <- "black"
    .theme_ss <- theme_bw(base_size=12) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8 * scale, 
                                       family="mono", face='bold', color=color.axis),
            axis.text.y = element_text(hjust = 0.5,size=12 * scale, family="mono",face='bold',color=color.axis),
            axis.text = element_text(size = 12 * scale, family = "mono",color=color.axis),
            axis.title=element_text(face="bold", size=12 * scale,color="black"),
            plot.title=element_text(face="bold", size=12 * scale))
    
    # Plot W matrix (feature activities)
    W_hc <- hclust(dist(W, method="euclidean"), method="ward.D")
    W_variant.ordering <- W_hc$labels[W_hc$order]
    W_plt_df <- W %>%
      as.data.frame() %>%
      rownames_to_column(var="variant") %>%
      gather(key="cluster", value="activity", -variant) %>%
      mutate(variant=factor(variant, levels=W_variant.ordering),
             cluster=factor(cluster, 
                            levels=paste0("V", 1:ncol(W))))
    W_plt <- ggplot(W_plt_df, aes(x=variant, y=cluster, fill=activity)) + 
      geom_tile() +
      scale_fill_gradient2(low="white", high ="black", name=paste("Activity", sep="")) +
      #p = p + scale_fill_gradientn(values=c(0,0.1,0.2,0.5,0.7,1.0),colours=c("yellow","green","black","red","magenta"),limit=c(0,1.0))
      .theme_ss +
      ggtitle(paste0("Variant Association to Clusters (k=", k, ")")) +
      ylab("Cluster") + xlab("Variant") +
      theme(axis.title.x = element_text(face="bold",colour="black", size=12 * scale0)) +
      theme(axis.title.y = element_text(face="bold",colour="black", size=12 * scale0)) +
      theme(legend.position="right") +
      theme(legend.key.size = unit(0.5, "cm"))
    ggsave(paste0(dir_save,"W_plot_K", k, ".pdf"), plot=W_plt)
    
    H_hc <- hclust(dist(t(H), method="euclidean"), method="ward.D")
    H_trait.ordering <- H_hc$labels[H_hc$order]
    H_plt_df <- t(H) %>%
      as.data.frame() %>%
      rownames_to_column(var="trait") %>%
      gather(key="cluster", value="activity", -trait) %>%
      mutate(cluster=factor(cluster, levels=paste0("V", 1:nrow(H))),
             trait=factor(trait, levels=H_trait.ordering))
    H_plt <- ggplot(H_plt_df, aes(x=trait, y=cluster, fill=activity)) + 
      geom_tile() +
      scale_fill_gradient2(low="white", high ="black", name=paste("Activity", sep="")) +
      #p = p + scale_fill_gradientn(values=c(0,0.1,0.2,0.5,0.7,1.0),colours=c("yellow","green","black","red","magenta"),limit=c(0,1.0))
      .theme_ss +
      ggtitle(paste0("Variant Association to Clusters (k=", k, ")")) +
      ylab("Cluster") + xlab("Trait") +
      theme(axis.title.x=element_text(face="bold", colour="black", size=12 * scale0)) +
      theme(axis.title.y=element_text(face="bold", colour="black", size=12 * scale0)) +
      theme(legend.position="right") +
      theme(legend.key.size = unit(0.5, "cm"))
    ggsave(paste0(dir_save,"H_plot_K", k, ".pdf"), plot=H_plt)
  })
}


# post-hoc functions

calculate_cutoff <- function(w) {
  
  w_mat <- data.frame(t(w)) %>%
    set_colnames(.["variant", ]) %>%
    subset(!rownames(.) %in% "variant")
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
  w_list <- as.list(w_mat)
  weight = unlist(w_list)
  
  
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


v2g <- function(gwas_data, prefix, output_dir) {
  
  # gwas_data <- vroom(gwas_path)
  # prefix <- str_replace(basename(gwas_path), '.tsv.gz', '')
  
  # Sanitize prefix for file naming
  prefix <- gsub("[^a-zA-Z0-9_-]", "_", prefix)
  
  all_combine <- data.frame()
  for (chr in 1:22) {
    cat('Working on chr', chr, '\n')
    
    # Use a temporary variable to avoid modifying gwas_data in place
    gwas_chr_data <- gwas_data %>% filter(chromosome == chr)
    
    # Load the index data for the current chromosome
    index_file <- paste0("/humgen/diabetes2/users/satoshi/database/opentarget/V2G/v2g_index_perchr/chr", chr, ".v2g_index.tsv.gz")
    index_data <- vroom(index_file) %>%
      dplyr::rename(position_v2g = position, rsid = rs_id)
    
    # Load the scored data for the current chromosome
    scored_file <- paste0("/humgen/diabetes2/users/satoshi/database/opentarget/V2G/v2g_scored_perchr/chr", chr, ".v2g_scored.tsv.gz")
    scored_data <- vroom(scored_file) %>%
      dplyr::transmute(chromosome = chr_id, position_hg38 = position, gene_id, overall_score, source_score_list)
    
    # Merge the GWAS data with the index data by "rsid"
    merged_data <- left_join(gwas_chr_data, index_data %>%
                               dplyr::rename(chromosome = chr_id, position_hg38 = position_v2g), by = c('chromosome', 'rsid'))
    
    # Merge the merged data with the scored data by "chromosome" and "position_hg38"
    merged_data2 <- left_join(merged_data, scored_data, by = c("chromosome", "position_hg38"))
    
    # Select the top V2G data based on the highest overall score
    merged_v2g_data3 <- merged_data2 %>%
      group_by(rsid) %>%
      arrange(desc(overall_score)) %>%
      filter(!duplicated(rsid)) %>%
      ungroup()
    
    ########################
    # Finally, we'll add gene symbol based on ENSG
    ########################
    gene_symbol <- vroom('/humgen/diabetes2/users/satoshi/database/ensembl_gene/grch38_df.tsv')
    gene_symbol %<>% transmute(gene_id = ensgene, gene = symbol)
    
    # Join with gene symbol data
    combine_df <- left_join(merged_v2g_data3, gene_symbol, by = 'gene_id')
    
    # Append to all_combine using bind_rows for better performance
    all_combine <- bind_rows(all_combine, combine_df)
  }
  
  # Ensure output directory exists
  # dir.create("output", showWarnings = FALSE, recursive = TRUE)
  
  # Write the combined data to a gzipped TSV file
  new_file = paste0(prefix, '.v2g.tsv.gz')
  write_tsv(all_combine, file = file.path(output_dir, new_file))
  print("Done with V2G!")
}

do_post_analysis <- function(main_dir,
                             my_chain,
                             df_gwas,
                             k=NULL,
                             cluster_names=NULL,
                             loci_file="query",
                             make_excel = T,
                             do_liftover = T,
                             do_v2g = F) {
  
  # 1.) load data ---
  
  
  # find majority K from run summary
  df_run_summary <- fread(file.path(main_dir,"run_summary.txt"),
                          stringsAsFactors = F, data.table = F)
  
  k_counts <- df_run_summary %>%
    dplyr::count(K) %>%
    mutate(perc = n/nrow(df_run_summary))
  
  k <- k_counts$K[which.max(k_counts$n)]
  
  if (is.null(cluster_names)) {
    # cluster_names <- names(w)[names(w) %like% "X"]
    cluster_names <- paste0("X",1:k)
    
  } else {
    colnames(w) <- c(cluster_names, "variant")
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

  cutoff <- calculate_cutoff(w)
  print(sprintf("Optimal weight cutoff is %.3f!", cutoff))
  
  
  
  # 2.) get nearest gene ---
  
  df_rsIDs <- fread(file.path(main_dir,"rsID_map.txt"),
                    data.table=F, stringsAsFactors=F) %>%
    mutate(variant = gsub("_",":",str_before_nth(VAR_ID, "_", 2))) %>%
    filter(!duplicated(variant)) %>%
    dplyr::select(VAR_ID, rsID, variant)
  
  if (make_excel==T) {
    if (loci_file=="query") {
      # print("Querying for locus names...")
      library(GenomicRanges)
      library(Homo.sapiens)
      
      geneRanges <- function(db, column="ENTREZID"){
        g <- genes(db, columns=column)
        col <- mcols(g)[[column]]
        genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
        mcols(genes)[[column]] <- as.character(unlist(col))
        genes
      }
      
      gns = geneRanges(Homo.sapiens, column="SYMBOL")
      # need chr | start | end, with gene names for row names
      df_gns <- data.frame(gns)
      gtf.gene <- df_gns %>%
        mutate(chr = gsub("chr","",seqnames)) %>%
        dplyr::select(chr, start, end, SYMBOL)
      
      # need chr | start | end, with gene names for row names
      df_entrez <- geneRanges(Homo.sapiens, column="ENTREZID") %>%
        data.frame() %>%
        mutate(chr = gsub("chr","",seqnames)) %>%
        mutate(ChrPos = paste(chr,start,sep=":")) %>%
        dplyr::select(ChrPos, ENTREZID)
      
      #' Convert from string to range
      #'
      #' @param pos A vector of strings ex. chr1 2938302 2938329
      #' @param delim Delimiter for string splitting
      #' @param region Boolean of whether region or just one position
      #'
      #' @returns Dataframe of ranges
      #'
      string2range <- function(pos, delim=' ', region=TRUE) {
        posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
        posp[,1] <- posp[,1]
        posp[,2] <- as.numeric(as.character(posp[,2]))
        if(region) {
          posp[,3] <- as.numeric(as.character(posp[,3]))
        } else {
          posp[,3] <- posp[,2]
        }
        return(posp)
      }
      
      #' Convert from ranges to GRanges
      #'
      #' @param df Dataframe with columns as sequence name, start, and end
      #'
      #' @returns GRanges version
      #'
      range2GRanges <- function(df) {
        require(GenomicRanges)
        require(IRanges)
        gr <- GenomicRanges::GRanges(
          seqnames = df[,1],
          ranges=IRanges(start = df[,2], end = df[,3])
        )
        return(gr)
      }
      
      # convert SNPs to GRanges
      snps <- w$variant
      snps.ranges <- string2range(snps, delim=":", region=FALSE)
      snps.granges <- range2GRanges(snps.ranges)
      names(snps.granges) <- snps
      
      # convert genes to GRanges
      gtf.granges <- range2GRanges(gtf.gene)
      names(gtf.granges) <-  gtf.gene$SYMBOL  #gene.names
      
      hits <- GenomicRanges::nearest(snps.granges, gtf.granges)
      # make vector of SNPs to gene
      
      df_gns <- df_gns %>%
        mutate(chr = gsub("chr","",seqnames)) %>%
        mutate(ChrPos = paste(chr,start,sep=":")) %>%
        dplyr::select(chr, start, end, ChrPos, SYMBOL) %>%
        merge(df_entrez, by="ChrPos")
      
      df_hits <- data.frame(gene=names(gtf.granges)[hits]) %>%
        mutate(variant = names(snps.granges)) %>%
        merge(df_gns[,c("SYMBOL","ENTREZID")], by.x="gene", by.y="SYMBOL") %>%
        filter(!duplicated(variant))
      
      w_wLoci <- w %>%
        left_join(df_hits, by="variant") %>%
        mutate(gene = dplyr::coalesce(gene, variant)) # coalese fills in empty genes w/ variant
      
    } else if (is.character(loci_file) & !loci_file %in% c("query","ChrPos")) {
      
      # print("Mapping manually curated locus names to SNPs...")
      locus_map <- read.csv(loci_file, header = T, stringsAsFactors = F) %>%
        mutate(variant = paste(CHR,POS,sep=":"))
      
      w_wLoci <- w %>%
        merge(locus_map[c("variant","Locus")], by="variant") %>%
        dplyr::rename(gene=Locus) %>%
        mutate(gene=make.unique(gene))
    } else {
      # print("Using Chr:Pos for locus names...")
      w_wLoci <- w %>%
        mutate(gene=variant)
    }
    
    # make excel
    library(openxlsx)
    library(tidyverse)
    
    h.exc <- t(h) %>%
      data.frame() %>%
      mutate(trait = row.names(.)) %>%
      relocate(trait)
  
    w.exc <- w_wLoci %>%
      left_join(df_rsIDs, by = "variant") %>%
      dplyr::select_if(!names(.) %in% c('ENTREZID', 'variant')) %>%
      relocate(VAR_ID, rsID, gene) 
    last_id_col <- ncol(w.exc)-length(cluster_names)
    
    
    posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
    
    # for (i in 1:length(colnames(w.exc)){
    wb <- openxlsx::createWorkbook()
    for (i in 1:length(cluster_names)){
      w_temp <- w.exc %>%
        dplyr::arrange(-!!as.symbol(cluster_names[i]))
      h_temp <- h.exc %>%
        arrange(-!!as.symbol(cluster_names[i]))
      
      cur_sheet = sprintf("cluster%i",i)
      addWorksheet(wb, cur_sheet)
      start_col1 <- 1
      start_col2 <- ncol(w_temp)+3
      # cutoff_col <- 2*(ncol(w_temp)+2)+1
      
      writeData(wb=wb, sheet=cur_sheet, x=w_temp, startCol = start_col1, startRow = 1,
                rowNames=F,colNames = T)
      writeData(wb=wb, sheet=cur_sheet, x=h_temp, startCol = start_col2, startRow = 1,
                rowNames=F,colNames = T)
      # writeData(wb=wb, sheet="Foo", x=data.frame(cutoff=cutoff), startCol = cutoff_col, startRow = 1)
      conditionalFormatting(wb, cur_sheet,
                            cols = last_id_col, rows = 2:(nrow(w_temp)+1),
                            rule =sprintf("%s2>=%.7f",
                                          int2col(last_id_col+i),cutoff),
                            style = posStyle)
      conditionalFormatting(wb, cur_sheet,
                            cols = start_col2, rows = 2:(nrow(h_temp)+1),
                            rule = sprintf("%s2>=%.7f",
                                           int2col(start_col2+i),cutoff),
                            style = posStyle)
    }
    fp <- file.path(main_dir, sprintf("sorted_cluster_weights_K%i_rev.xlsx",k))
    openxlsx::saveWorkbook(wb, fp, overwrite = T)
  }
  
  # 3.) hg19-to-hg38 liftover --
  
  if (do_liftover==T) {
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
    
    print("Variants not found in liftover...")
    df_weights %>%
      filter(!group %in% df_liftover$group)
    
    
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
  
  
  # V2G

  if (do_v2g) {
    v2g_input <- df_gwas %>%
      separate(ChrPos, into=c("chromosome","position"),sep=":",remove = T) %>%
      mutate_at(vars(chromosome, position), as.integer) %>%
      dplyr::select(chromosome, position, rsid=rsID)
    
    v2g(gwas_data = v2g_input,
           prefix=paste(version,"v2g",sep = "_"),
           output_dir = main_dir
    )
  } 

  
  return(list(df_run_summary=df_run_summary, k=k, cluster_names=cluster_names, w=w_wLoci, h=h, cutoff=cutoff))
  
}
