library(tidyverse)


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
