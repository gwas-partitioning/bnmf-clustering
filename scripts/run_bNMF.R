library(tidyverse)


############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:
####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.
####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.
####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.
#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

######################################################################################################
####### Bayesian NMF algorithms for clustering
######################################################################################################
####### For implementation details see the ppaer 
####### Udler MS, Kim J, von Grotthuss M,
####### Bonàs-Guarch S, Cole JB, Chiou J, et al. (2018)
####### Type 2 diabetes genetic loci informed by multi-trait
####### associations point to disease mechanisms and
####### subtypes: A soft clustering analysis. PLoS Med 15
####### (9): e1002654.
#################################
####### For details on the original algorithms 
####### see Tan, V.Y. & Févotte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence.
####### IEEE Trans. Pattern Anal. Mach. Intell. 35, 1592–1605 (2013).
######################################################################################################

BayesNMF.L2EU <- function(
  V0, n.iter=10000, a0=10, tol=1e-7, K=10, K0=10, phi=1.0
) {
  
  # Bayesian NMF with half-normal priors for W and H
  # V0: input z-score matrix (variants x traits) *****
  # n.iter: Number of iterations for parameter optimization
  # a0: *****
  # tol: Tolerance for convergence of fitting procedure
  # K: Number of clusters to be initialized (algorithm may drive some to zero)
  # K0: *****
  # phi: *****
  
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0-min(V0)
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
    n.like,  # List of reconstruction errors (sum of squared errors / 2) per iteration *****
    n.evid,  # List of negative log-likelihoods per iteration
    n.lambda,  # List of lambda vectors (shared weights for each of K clusters, some ~0) per iteration
    n.error  # List of reconstruction errors (sum of squared errors) per iteration
  ))
}
