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

###########################
###########################
##### Bayesian NMF with half-normal priors for W and H
BayesNMF.L2EU <- function(V0,n.iter,a0,tol,K,K0,phi) {
        eps <- 1.e-50
        del <- 1.0
        active_nodes <- colSums(V0) != 0
        V0 <- V0[,active_nodes]
        V <- V0-min(V0)
        Vmin <- min(V)
        Vmax <- max(V)
        N <- dim(V)[1]
        M <- dim(V)[2]

        W <- matrix(runif(N * K)*Vmax,ncol=K)
        H <- matrix(runif(M * K)*Vmax,ncol=M)
        I <- array(1,dim=c(N,M))
        V.ap <- W%*%H+eps

        phi <- sd(V)^2*phi
        C <- (N+M)/2+a0+1
        b0 <- 3.14*(a0-1)*mean(V)/(2*K0)
        lambda.bound <- b0/C
        lambda <- (0.5*colSums(W^2)+0.5*rowSums(H^2)+b0)/C
        lambda.cut <- lambda.bound*1.5

        n.like <- list()
        n.evid <- list()
        n.error <- list()
        n.lambda <- list()
        n.lambda[[1]] <- lambda
        iter <- 2
        count <- 1
        while (del >= tol & iter < n.iter) {
                H <- H*(t(W)%*%V)/(t(W)%*%V.ap+phi*H*matrix(rep(1/lambda,M),ncol=M)+eps)
                V.ap <- W %*% H + eps
                W <- W*(V%*%t(H))/(V.ap%*%t(H)+phi*W*t(matrix(rep(1/lambda,N),ncol=N))+eps)
                V.ap <- W %*% H + eps
                lambda <- (0.5*colSums(W^2)+0.5*rowSums(H^2)+b0)/C
                del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
                like <- sum((V-V.ap)^2)/2
                n.like[[iter]] <- like
                n.evid[[iter]] <- like + phi*sum((0.5*colSums(W^2)+0.5*rowSums(H^2)+b0)/lambda+C*log(lambda))
                n.lambda[[iter]] <- lambda
                n.error[[iter]] <- sum((V-V.ap)^2)
                if (iter %% 100 == 0) {
                        cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
                }
                iter <- iter+1
        }
        return(list(W,H,n.like,n.evid,n.lambda,n.error))
}

plot.heatmap.ggplot.new <- function(mat) {
	scale0 <- 0.8
        scale <- 1
        g.ordering <- c("G4","G3","G2","G1")
        color.axis <- "black"
        .theme_ss <- theme_bw(base_size=12) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono",face='bold',color=color.axis),
                axis.text.y = element_text(hjust = 0.5,size=8*scale, family="mono",face='bold',color=color.axis),
                axis.text = element_text(size = 12*scale, family = "mono",color=color.axis),
                axis.title=element_text(face="bold",size=12*scale,color="black"),
                plot.title=element_text(face="bold",size=12*scale))
        mat[mat < 1.e-10] <- 0
        hc <- hclust(dist(mat,method="euclidean"),method="ward.D")
        feature.ordering <- hc$labels[hc$order]
        df <- melt(mat)
        colnames(df) <- c("feature","signature","activity")
        #df$feature <- factor(df$feature,levels=feature.ordering)
        #df$signature <- factor(df$signature,levels=c("W4","W3","W2","W1"))
        p = ggplot(df,aes(y=feature,x=signature,fill=activity))+geom_tile() #geom_tile(colour="yellow")
        p = p + scale_fill_gradient2(low="white",high ="black",name=paste("Activity",sep=""))
        #p = p + scale_fill_gradientn(values=c(0,0.1,0.2,0.5,0.7,1.0),colours=c("yellow","green","black","red","magenta"),limit=c(0,1.0))
        p = p + .theme_ss
        p = p + ggtitle("Feature Assoication to Clusters")
        #p = p + ylab("Contributions") + xlab("Feature")
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=12*scale0))
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=12*scale0))
        p = p + theme(legend.position="right")
        p = p + theme(legend.key.size = unit(0.5, "cm"))
        #pdf(file=paste(OUTPUT,"feature.association_to_clusters.pdf",sep=""),width=10,height=4)
                plot(p)
        #dev.off()
}

plot.heatmap.2 <- function(x,rowTF,colTF) {
        s1 <- 0.75
        s2 <- 1.0
        s3 <- 1.5
        mydist <- function(c) {dist(c,method="euclidean")}
        myclust <- function(c) {hclust(c,method="ward.D")}
        heatmap.2(as.matrix(x), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both",margins=c(8,8),
                Rowv=rowTF, Colv=colTF, symbreaks=F, key=TRUE, symkey=F,
                density.info="none", trace="none",labCol=colnames(x),labRow=rownames(x),col=greenred(40),cex.lab=s1,cexRow=0.7,cexCol=0.7,keysize=s1)
}

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(reshape2)

#CURRENT <- paste(getwd(),"/",sep="")
#OUTPUT <- paste(CURRENT,"OUTPUT/DIST_test/",sep="")
#system(paste("mkdir",OUTPUT,sep=" "))

##### mat.all: the trait by genotype matrix; the positive and negative association of each trait to genes are separately handled with distinct features.
##### mat.all in this scipt is not exactly same as the one used in the paper and we included it as an example input.
#load("t2d_data.example.RData")
#mat_all <- mat.all
# load data
library(tidyverse)
setwd("C:/Users/hk745/Dropbox (Partners HealthCare)")

CURRENT <- paste(getwd(),"/",sep="")
OUTPUT <- paste(CURRENT,"OUTPUT/DIST_test/T2D/trait_test/",sep="")

#inputtraits <- read.delim("C:/Users/hk745/Dropbox (Partners HealthCare)/traitlist_pval_0.05_snpcnt_filterGIANT.txt",stringsAsFactors = FALSE, header = F)
#inputtraits$V1 <- gsub("-", ".", inputtraits$V1)


#t2d_snps <- read.delim("C:/Users/hk745/Dropbox (Partners HealthCare)/T2D_98snps_23traits.txt",stringsAsFactors = FALSE)
#filter <- t2d_snps
#t2d_snps <- read.delim("C:/Users/hk745/Dropbox (Partners HealthCare)/T2D_98snps_noloop_minusCHOL.txt",stringsAsFactors = FALSE)
#t2d_snps <- t2d_snps[,names(t2d_snps) %in% names(filter)]
#gbd_snps <- t2d_snps


#filtered <- read.delim("C:/Users/hk745/Dropbox (Partners HealthCare)/T2D_filter.txt",stringsAsFactors = FALSE, header = F)
#filtered <- filtered[filtered$V2<0.05,]

gbd_snps <- read.delim("C:/Users/hk745/Dropbox (Partners HealthCare)/T2Donly_eu_45traits_412snps.txt",stringsAsFactors = FALSE)

#gbd_snps <- gbd_snps[gbd_snps$VAR_ID_hg19 %in% filtered$V1,]

#gbd_snps$GLGC_dv2.CHOL.ZN <- NULL
gbd_t <- gbd_snps
rownames(gbd_t) <- gbd_snps$locus
gbd_t <- gbd_t[,-c(1:4)]
gbd_t <- t(gbd_t)
gbd_t[gbd_t == "."] <- NA



#gbd_t <- gbd_t[rownames(gbd_t) %in% inputtraits$V1,]
#write.table(rownames(gbd_t),file=paste(CURRENT,paste("traitlist_0.05_snpcnt.txt",sep="."),sep=""), append = F, quote = F, sep = "\t",
#            eol = "\n", na = "NA", dec = ".", row.names = F,
 #           col.names = F, qmethod = c("escape", "double"))

gbd_t <- as.data.frame(gbd_t,stringsAsFactors = FALSE)
na <- as.data.frame(as.data.frame(map(gbd_t, ~sum(is.na(.)))))
nadrop <- na[,na>dim(gbd_t)[1]*0.2, drop = FALSE]  #### 
test <- gbd_t
test <- test[,!(colnames(test) %in% colnames(nadrop))]
gbd_t <- test

#gbd <- as.data.frame(t(gbd_t))
#nat <- as.data.frame(as.data.frame(map(gbd, ~sum(is.na(.)))))
#nadropt <- nat[,nat>dim(gbd)[1]*0.2]

gbd_t[is.na(gbd_t)] <- 0
mat.neg <- gbd_t
mat.pos <- as.data.frame(apply(mat.neg, 1, function(x) ifelse(as.numeric(x)>=0,as.numeric(x),0)))
#mat.pos <- as.data.frame(apply(mat.neg, 1, function(x) ifelse(as.numeric(x)<0,as.numeric(x),0)))
mat.pos <- t(mat.pos)
colnames(mat.pos) <- colnames(mat.neg)
mat.neg <- as.data.frame(apply(mat.neg, 1, function(x) ifelse(as.numeric(x)>=0,0,as.numeric(x))))
#mat.neg <- as.data.frame(apply(mat.neg, 1, function(x) ifelse(as.numeric(x)<0,0,as.numeric(x))))
mat.neg <- t(mat.neg)
colnames(mat.neg) <- colnames(mat.pos)
mat.neg <- mat.neg*(-1)
#mat.pos <- mat.pos*(-1)
rownames(mat.pos) <- paste(rownames(mat.pos), "pos", sep = "_")
rownames(mat.neg) <- paste(rownames(mat.neg), "neg", sep = "_")
mat.all <- rbind(mat.pos,mat.neg)



##### simple hierarchical clustering ############## not working
library(pheatmap)
hc1 <- hclust(dist(mat.all,method="euclidean"),method="ward.D")
hc1.ordering <- hc1$labels[hc1$order]
hc2 <- hclust(dist(t(mat.all),method="euclidean"),method="ward.D")
hc2.ordering <- hc2$labels[hc2$order]
order1 <- match(hc1.ordering,rownames(mat.all),nomatch=0)
order2 <- match(hc2.ordering,colnames(mat.all),nomatch=0)
pdf(file=paste(OUTPUT,"hierarchicalc.mat.pdf",sep=""),width=40,height=10)
fontsize_row = 8
fontsize_col = 0.5
pheatmap(mat.all[order1,order2], fontsize_col = fontsize_col, fontsize_row = fontsize_row)+ scale_fill_gradient2(low="white",high ="black",name=paste("Activity",sep=""))
#plot.heatmap.2(mat.all[order1,order2],F,F)
dev.off()



##### Running Bayesian NMF with half-normal priors for W and H 
if (TRUE) {
#n.iter <- 100 ### number of runs.. you can increase this number
n.iter <- 30 ### number of runs.. you can increase this number

for (i in 1:n.iter) {
        res <- BayesNMF.L2EU(as.matrix(mat.all),200000,10,1.e-07,10,10,1.0)
        save(res,file=paste(OUTPUT,paste("res.L2EU.Bayes",i,"RData",sep="."),sep=""))
}
}

tmpK <- rep(0,n.iter)
tmpE <- rep(0,n.iter)
tmpRUN <- rep(0,n.iter)
for (i in 1:n.iter) {
        load(file=paste(OUTPUT,paste("res.L2EU.Bayes",i,"RData",sep="."),sep=""))
        lambda <- res[[5]]
        lambda <- unlist(lambda[length(lambda)])
        lambda <- lambda-min(lambda)
        tmpK[i] <- sum(lambda > 0)
	evid <- res[[4]]
	tmpE[i] <- evid[length(evid)]
	tmpRUN[i] <- i
}
df.run <- data.frame(tmpK,unlist(tmpE),tmpRUN)
colnames(df.run) <- c("K","evid","run")
df.run <- df.run[order(df.run$evid,decreasing=T),] 
write.table(df.run,file=paste(OUTPUT,paste("summary.run.txt",sep="."),sep=""), append = F, quote = F, sep = "\t",
        eol = "\n", na = "NaN", dec = ".", row.names = F,
        col.names = T, qmethod = c("escape", "double"))

#### df.run - the summary data-frame for bNMF runs; K = number of clusters, evid = -log(posterior), run = the index of bNMF run
#### How to choose K: (i) We usually perfer the most probable K. For example, here 57% K=5 and 43% K=4, so we will consider K=5.
#### (ii) After selcting K then look at "evid" for all runs with the selected K (here K=5) and choose the run with the lowest "evid" 
#### corresponding to the maximum posterior solution
#### (iii) Sometimes you may need a manual inspection for other solutions based on your prior knowledge or biological consideration. 
#### Specially, when your most probable solution corresponds to the lowest K, it is recommended to examine the solution with (K+1) and check which solution 
#### is more biologically plausible. 

#### Below we will generate outputs of the maximum posterior solutions at different K
unique.K <- table(df.run$K)
n.K <- length(unique.K) ### number of distict K
MAP.K <- rep(0,n.K) ### bNMF run index with the maximum posterior for given K
for (i in 1:n.K) {
        tmp <- df.run[df.run$K==as.numeric(names(unique.K)[i]),]
        MAP.K[i] <- tmp$run[which.min(tmp$evid)]
}

for (m in 1:n.K) {

index.m <- as.numeric(names(unique.K)[m])

load(file=paste(OUTPUT,paste("res.L2EU.Bayes",MAP.K[m],"RData",sep="."),sep=""))
W <- res[[1]]
H <- res[[2]]
W <- W[,colSums(W)!=0]
H <- H[rowSums(H)!=0,]
colnames(W) <- paste("W",seq(1:ncol(W)),sep="")
rownames(H) <- colnames(W)
W[W < 1.e-10] <- 0 ### feature-cluster association matrix
H[H < 1.e-10] <- 0 ### cluster-gene association matrix

if (FALSE) {
	W.mid <- W
	H.mid <- H
	for (i in 1:ncol(W)) {
		H.mid[i,] <- H.mid[i,]*colSums(W)[i]
		W.mid[i,] <- W.mid[i,]*rowSums(H)[i]
	}
	W.norm <- apply(W.mid,2,function(x) x/sum(x))
	H.norm <- apply(H.mid,2,function(x) x/sum(x))
}

W0 <- data.frame(W)
W0[,"feature"] <- rownames(W)
H0 <- data.frame(H)
H0[,"cluster"] <- rownames(H)

if (TRUE) {
write.table(W0,file=paste(OUTPUT,paste("L2EU.W.mat",index.m,"txt",sep="."),sep=""), append = F, quote = F, sep = "\t",
        eol = "\n", na = "NaN", dec = ".", row.names = F,
        col.names = T, qmethod = c("escape", "double"))
write.table(H0,file=paste(OUTPUT,paste("L2EU.H.mat",index.m,"txt",sep="."),sep=""), append = F, quote = F, sep = "\t",
        eol = "\n", na = "NaN", dec = ".", row.names = F,
        col.names = T, qmethod = c("escape", "double"))
}

mat.reconstructed <- W%*%H ### reconstructed matrix == approximation for the input matrix 
#pdf(file=paste(OUTPUT,paste("L2EU.hc.mat.WH.0",index.m,"pdf",sep="."),sep=""),width=8,height=8)
#	plot.heatmap.2(mat.reconstructed[order1,order2],F,F)
#dev.off()

K <- ncol(W)
for (i in 1:K) {
	mat1 <- W[,i]%*%t(as.matrix(H[i,]))
	rownames(mat1) <- rownames(mat.all)
	#pdf(file=paste(OUTPUT,paste("hc.mat.WH",i,index.m,"pdf",sep="."),sep=""),width=8,height=8)
	#	plot.heatmap.2(mat1[order1,order2],F,F)
	#dev.off()
}

	scale0 <- 0.8
        scale <- 1
	g.ordering <- paste("G",seq(1:ncol(W)),sep="")
        color.axis <- "black"
        .theme_ss <- theme_bw(base_size=12) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono",face='bold',color=color.axis),
                axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono",face='bold',color=color.axis),
                axis.text = element_text(size = 12*scale, family = "mono",color=color.axis),
                axis.title=element_text(face="bold",size=12*scale,color="black"),
                plot.title=element_text(face="bold",size=12*scale))
        mat <- W
        mat[mat < 1.e-10] <- 0
        hc <- hclust(dist(mat,method="euclidean"),method="ward.D")
        feature.ordering <- hc$labels[hc$order]
        df <- melt(mat)
        colnames(df) <- c("feature","signature","activity")
        df$feature <- factor(df$feature,levels=feature.ordering)
        df$signature <- factor(df$signature,levels=paste("W",seq(1:ncol(W)),sep=""))
        p = ggplot(df,aes(x=feature,y=signature,fill=activity))+geom_tile() #geom_tile(colour="yellow")
	p = p + scale_fill_gradient2(low="white",high ="black",name=paste("Activity",sep=""))
        #p = p + scale_fill_gradientn(values=c(0,0.1,0.2,0.5,0.7,1.0),colours=c("yellow","green","black","red","magenta"),limit=c(0,1.0))
        p = p + .theme_ss
        p = p + ggtitle("Feature Assoication to Clusters")
        p = p + ylab("Contributions") + xlab("Feature")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=12*scale0))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=12*scale0))
        p = p + theme(legend.position="right")
        p = p + theme(legend.key.size = unit(0.5, "cm"))
	pdf(file=paste(OUTPUT,paste("L2EU.feature.association_to_clusters",index.m,"pdf",sep="."),sep=""),width=10,height=4)
		plot(p)
	dev.off()

	# size = 8*scale (original)
	scale0 <- 0.8
        scale <- 1
	g.ordering <- paste("G",seq(1:ncol(W)),sep="")
        color.axis <- "black"
        .theme_ss <- theme_bw(base_size=12) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono",face='bold',color=color.axis),
                axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono",face='bold',color=color.axis),
                axis.text = element_text(size = 12*scale, family = "mono",color=color.axis),
                axis.title=element_text(face="bold",size=12*scale,color="black"),
                plot.title=element_text(face="bold",size=12*scale))
        mat <- H
        hc <- hclust(dist(t(mat),method="euclidean"),method="ward.D")
        gene.ordering <- hc$labels[hc$order]
        df <- melt(mat)
        colnames(df) <- c("signature","gene","activity")
        df$signature <- factor(df$signature,levels=paste("W",seq(1:ncol(W)),sep=""))
        df$gene <- factor(df$gene,levels=gene.ordering)
        p = ggplot(df,aes(x=gene,y=signature,fill=activity))+geom_tile() #geom_tile(colour="yellow")
	p = p + scale_fill_gradient2(low="white",high ="black",name=paste("Activity",sep=""))
        #p = p + scale_fill_gradientn(values=c(0,q1,q2,q3,q4,1),colours=c("yellow","green","black","red","magenta"),limit=c(0,1))
        p = p + .theme_ss
        p = p + ggtitle("Gene Assoication to Clusters")
        p = p + ylab("Contributions") + xlab("Genes")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=12*scale0))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=12*scale0))
        p = p + theme(legend.position="right")
        p = p + theme(legend.key.size = unit(0.5, "cm"))
	pdf(file=paste(OUTPUT,paste("L2EU.gene.association_to_clusters",index.m,"pdf",sep="."),sep=""),width=10,height=4)
		plot(p)
	dev.off()
}

