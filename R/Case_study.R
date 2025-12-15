######################################################################################################################################################################################
########################################################################### Clean application in ASD dataset #########################################################################
######################################################################################################################################################################################

## Load required packages
library(doParallel)
library(lmtest)
library(pracma)
library(robustbase)
library(vars)
library(dynlm)
library(SuperCell)
library(MetacellAnalysisToolkit)
library(Seurat)
library(patterncausality)
library(dplyr)
library(igraph)
library(energy)
library(SCORPIUS)
library(e1071)    
library(caret)    
library(pROC)

## Load input data
load("ASD_lncR_mR.RData")

## Index of cause and effect for causal inference
cause <- 5134:5946
effect <- 1:5133

########################################################## Metacells construction #####################################################################
## Identifying metacells using SuperCell R package (https://github.com/GfellerLab/SuperCell)
# ASD
ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), ASD_SC10$membership)
colnames(ASD_SC10_GE) <- paste("MC", seq(ncol(ASD_SC10_GE)), sep = "")

ASD_SC20 <- supercell_rescale(ASD_SC10, gamma = 20)
ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), ASD_SC20$membership)
colnames(ASD_SC20_GE) <- paste("MC", seq(ncol(ASD_SC20_GE)), sep = "")

ASD_SC30 <- supercell_rescale(ASD_SC10, gamma = 30)
ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), ASD_SC30$membership)
colnames(ASD_SC30_GE) <- paste("MC", seq(ncol(ASD_SC30_GE)), sep = "")

ASD_SC40 <- supercell_rescale(ASD_SC10, gamma = 40)
ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), ASD_SC40$membership)
colnames(ASD_SC40_GE) <- paste("MC", seq(ncol(ASD_SC40_GE)), sep = "")

ASD_SC50 <- supercell_rescale(ASD_SC10, gamma = 50)
ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), ASD_SC50$membership)
colnames(ASD_SC50_GE) <- paste("MC", seq(ncol(ASD_SC50_GE)), sep = "")

# Normal
Normal_SC10 <- SCimplify(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)), Normal_SC10$membership)
colnames(Normal_SC10_GE) <- paste("MC", seq(ncol(Normal_SC10_GE)), sep = "")

Normal_SC20 <- supercell_rescale(Normal_SC10, gamma = 20)
Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)), Normal_SC20$membership)
colnames(Normal_SC20_GE) <- paste("MC", seq(ncol(Normal_SC20_GE)), sep = "")

Normal_SC30 <- supercell_rescale(Normal_SC10, gamma = 30)
Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)), Normal_SC30$membership)
colnames(Normal_SC30_GE) <- paste("MC", seq(ncol(Normal_SC30_GE)), sep = "")

Normal_SC40 <- supercell_rescale(Normal_SC10, gamma = 40)
Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)), Normal_SC40$membership)
colnames(Normal_SC40_GE) <- paste("MC", seq(ncol(Normal_SC40_GE)), sep = "")

Normal_SC50 <- supercell_rescale(Normal_SC10, gamma = 50)
Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)), Normal_SC50$membership)
colnames(Normal_SC50_GE) <- paste("MC", seq(ncol(Normal_SC50_GE)), sep = "")

# ACC_ASD
ACC_ASD_SC10 <- SCimplify(as.matrix(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ACC_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), ACC_ASD_SC10$membership)
colnames(ACC_ASD_SC10_GE) <- paste("MC", seq(ncol(ACC_ASD_SC10_GE)), sep = "")

ACC_ASD_SC20 <- supercell_rescale(ACC_ASD_SC10, gamma = 20)
ACC_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), ACC_ASD_SC20$membership)
colnames(ACC_ASD_SC20_GE) <- paste("MC", seq(ncol(ACC_ASD_SC20_GE)), sep = "")

ACC_ASD_SC30 <- supercell_rescale(ACC_ASD_SC10, gamma = 30)
ACC_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), ACC_ASD_SC30$membership)
colnames(ACC_ASD_SC30_GE) <- paste("MC", seq(ncol(ACC_ASD_SC30_GE)), sep = "")

ACC_ASD_SC40 <- supercell_rescale(ACC_ASD_SC10, gamma = 40)
ACC_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), ACC_ASD_SC40$membership)
colnames(ACC_ASD_SC40_GE) <- paste("MC", seq(ncol(ACC_ASD_SC40_GE)), sep = "")

ACC_ASD_SC50 <- supercell_rescale(ACC_ASD_SC10, gamma = 50)
ACC_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), ACC_ASD_SC50$membership)
colnames(ACC_ASD_SC50_GE) <- paste("MC", seq(ncol(ACC_ASD_SC50_GE)), sep = "")

# ACC_Normal
ACC_Normal_SC10 <- SCimplify(as.matrix(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ACC_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), ACC_Normal_SC10$membership)
colnames(ACC_Normal_SC10_GE) <- paste("MC", seq(ncol(ACC_Normal_SC10_GE)), sep = "")

ACC_Normal_SC20 <- supercell_rescale(ACC_Normal_SC10, gamma = 20)
ACC_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), ACC_Normal_SC20$membership)
colnames(ACC_Normal_SC20_GE) <- paste("MC", seq(ncol(ACC_Normal_SC20_GE)), sep = "")

ACC_Normal_SC30 <- supercell_rescale(ACC_Normal_SC10, gamma = 30)
ACC_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), ACC_Normal_SC30$membership)
colnames(ACC_Normal_SC30_GE) <- paste("MC", seq(ncol(ACC_Normal_SC30_GE)), sep = "")

ACC_Normal_SC40 <- supercell_rescale(ACC_Normal_SC10, gamma = 40)
ACC_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), ACC_Normal_SC40$membership)
colnames(ACC_Normal_SC40_GE) <- paste("MC", seq(ncol(ACC_Normal_SC40_GE)), sep = "")

ACC_Normal_SC50 <- supercell_rescale(ACC_Normal_SC10, gamma = 50)
ACC_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), ACC_Normal_SC50$membership)
colnames(ACC_Normal_SC50_GE) <- paste("MC", seq(ncol(ACC_Normal_SC50_GE)), sep = "")

# PFC_ASD
PFC_ASD_SC10 <- SCimplify(as.matrix(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
PFC_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), PFC_ASD_SC10$membership)
colnames(PFC_ASD_SC10_GE) <- paste("MC", seq(ncol(PFC_ASD_SC10_GE)), sep = "")

PFC_ASD_SC20 <- supercell_rescale(PFC_ASD_SC10, gamma = 20)
PFC_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), PFC_ASD_SC20$membership)
colnames(PFC_ASD_SC20_GE) <- paste("MC", seq(ncol(PFC_ASD_SC20_GE)), sep = "")

PFC_ASD_SC30 <- supercell_rescale(PFC_ASD_SC10, gamma = 30)
PFC_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), PFC_ASD_SC30$membership)
colnames(PFC_ASD_SC30_GE) <- paste("MC", seq(ncol(PFC_ASD_SC30_GE)), sep = "")

PFC_ASD_SC40 <- supercell_rescale(PFC_ASD_SC10, gamma = 40)
PFC_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), PFC_ASD_SC40$membership)
colnames(PFC_ASD_SC40_GE) <- paste("MC", seq(ncol(PFC_ASD_SC40_GE)), sep = "")

PFC_ASD_SC50 <- supercell_rescale(PFC_ASD_SC10, gamma = 50)
PFC_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), PFC_ASD_SC50$membership)
colnames(PFC_ASD_SC50_GE) <- paste("MC", seq(ncol(PFC_ASD_SC50_GE)), sep = "")

# PFC_Normal
PFC_Normal_SC10 <- SCimplify(as.matrix(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
PFC_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), PFC_Normal_SC10$membership)
colnames(PFC_Normal_SC10_GE) <- paste("MC", seq(ncol(PFC_Normal_SC10_GE)), sep = "")

PFC_Normal_SC20 <- supercell_rescale(PFC_Normal_SC10, gamma = 20)
PFC_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), PFC_Normal_SC20$membership)
colnames(PFC_Normal_SC20_GE) <- paste("MC", seq(ncol(PFC_Normal_SC20_GE)), sep = "")

PFC_Normal_SC30 <- supercell_rescale(PFC_Normal_SC10, gamma = 30)
PFC_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), PFC_Normal_SC30$membership)
colnames(PFC_Normal_SC30_GE) <- paste("MC", seq(ncol(PFC_Normal_SC30_GE)), sep = "")

PFC_Normal_SC40 <- supercell_rescale(PFC_Normal_SC10, gamma = 40)
PFC_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), PFC_Normal_SC40$membership)
colnames(PFC_Normal_SC40_GE) <- paste("MC", seq(ncol(PFC_Normal_SC40_GE)), sep = "")

PFC_Normal_SC50 <- supercell_rescale(PFC_Normal_SC10, gamma = 50)
PFC_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), PFC_Normal_SC50$membership)
colnames(PFC_Normal_SC50_GE) <- paste("MC", seq(ncol(PFC_Normal_SC50_GE)), sep = "")

# Lower18_ASD
Lower18_ASD_SC10 <- SCimplify(as.matrix(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Lower18_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), Lower18_ASD_SC10$membership)
colnames(Lower18_ASD_SC10_GE) <- paste("MC", seq(ncol(Lower18_ASD_SC10_GE)), sep = "")

Lower18_ASD_SC20 <- supercell_rescale(Lower18_ASD_SC10, gamma = 20)
Lower18_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), Lower18_ASD_SC20$membership)
colnames(Lower18_ASD_SC20_GE) <- paste("MC", seq(ncol(Lower18_ASD_SC20_GE)), sep = "")

Lower18_ASD_SC30 <- supercell_rescale(Lower18_ASD_SC10, gamma = 30)
Lower18_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), Lower18_ASD_SC30$membership)
colnames(Lower18_ASD_SC30_GE) <- paste("MC", seq(ncol(Lower18_ASD_SC30_GE)), sep = "")

Lower18_ASD_SC40 <- supercell_rescale(Lower18_ASD_SC10, gamma = 40)
Lower18_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), Lower18_ASD_SC40$membership)
colnames(Lower18_ASD_SC40_GE) <- paste("MC", seq(ncol(Lower18_ASD_SC40_GE)), sep = "")

Lower18_ASD_SC50 <- supercell_rescale(Lower18_ASD_SC10, gamma = 50)
Lower18_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), Lower18_ASD_SC50$membership)
colnames(Lower18_ASD_SC50_GE) <- paste("MC", seq(ncol(Lower18_ASD_SC50_GE)), sep = "")

# Lower18_Normal
Lower18_Normal_SC10 <- SCimplify(as.matrix(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Lower18_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), Lower18_Normal_SC10$membership)
colnames(Lower18_Normal_SC10_GE) <- paste("MC", seq(ncol(Lower18_Normal_SC10_GE)), sep = "")

Lower18_Normal_SC20 <- supercell_rescale(Lower18_Normal_SC10, gamma = 20)
Lower18_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), Lower18_Normal_SC20$membership)
colnames(Lower18_Normal_SC20_GE) <- paste("MC", seq(ncol(Lower18_Normal_SC20_GE)), sep = "")

Lower18_Normal_SC30 <- supercell_rescale(Lower18_Normal_SC10, gamma = 30)
Lower18_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), Lower18_Normal_SC30$membership)
colnames(Lower18_Normal_SC30_GE) <- paste("MC", seq(ncol(Lower18_Normal_SC30_GE)), sep = "")

Lower18_Normal_SC40 <- supercell_rescale(Lower18_Normal_SC10, gamma = 40)
Lower18_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), Lower18_Normal_SC40$membership)
colnames(Lower18_Normal_SC40_GE) <- paste("MC", seq(ncol(Lower18_Normal_SC40_GE)), sep = "")

Lower18_Normal_SC50 <- supercell_rescale(Lower18_Normal_SC10, gamma = 50)
Lower18_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), Lower18_Normal_SC50$membership)
colnames(Lower18_Normal_SC50_GE) <- paste("MC", seq(ncol(Lower18_Normal_SC50_GE)), sep = "")

# Larger18_ASD
Larger18_ASD_SC10 <- SCimplify(as.matrix(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Larger18_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), Larger18_ASD_SC10$membership)
colnames(Larger18_ASD_SC10_GE) <- paste("MC", seq(ncol(Larger18_ASD_SC10_GE)), sep = "")

Larger18_ASD_SC20 <- supercell_rescale(Larger18_ASD_SC10, gamma = 20)
Larger18_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), Larger18_ASD_SC20$membership)
colnames(Larger18_ASD_SC20_GE) <- paste("MC", seq(ncol(Larger18_ASD_SC20_GE)), sep = "")

Larger18_ASD_SC30 <- supercell_rescale(Larger18_ASD_SC10, gamma = 30)
Larger18_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), Larger18_ASD_SC30$membership)
colnames(Larger18_ASD_SC30_GE) <- paste("MC", seq(ncol(Larger18_ASD_SC30_GE)), sep = "")

Larger18_ASD_SC40 <- supercell_rescale(Larger18_ASD_SC10, gamma = 40)
Larger18_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), Larger18_ASD_SC40$membership)
colnames(Larger18_ASD_SC40_GE) <- paste("MC", seq(ncol(Larger18_ASD_SC40_GE)), sep = "")

Larger18_ASD_SC50 <- supercell_rescale(Larger18_ASD_SC10, gamma = 50)
Larger18_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), Larger18_ASD_SC50$membership)
colnames(Larger18_ASD_SC50_GE) <- paste("MC", seq(ncol(Larger18_ASD_SC50_GE)), sep = "")

# Larger18_Normal
Larger18_Normal_SC10 <- SCimplify(as.matrix(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Larger18_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), Larger18_Normal_SC10$membership)
colnames(Larger18_Normal_SC10_GE) <- paste("MC", seq(ncol(Larger18_Normal_SC10_GE)), sep = "")

Larger18_Normal_SC20 <- supercell_rescale(Larger18_Normal_SC10, gamma = 20)
Larger18_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), Larger18_Normal_SC20$membership)
colnames(Larger18_Normal_SC20_GE) <- paste("MC", seq(ncol(Larger18_Normal_SC20_GE)), sep = "")

Larger18_Normal_SC30 <- supercell_rescale(Larger18_Normal_SC10, gamma = 30)
Larger18_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), Larger18_Normal_SC30$membership)
colnames(Larger18_Normal_SC30_GE) <- paste("MC", seq(ncol(Larger18_Normal_SC30_GE)), sep = "")

Larger18_Normal_SC40 <- supercell_rescale(Larger18_Normal_SC10, gamma = 40)
Larger18_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), Larger18_Normal_SC40$membership)
colnames(Larger18_Normal_SC40_GE) <- paste("MC", seq(ncol(Larger18_Normal_SC40_GE)), sep = "")

Larger18_Normal_SC50 <- supercell_rescale(Larger18_Normal_SC10, gamma = 50)
Larger18_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), Larger18_Normal_SC50$membership)
colnames(Larger18_Normal_SC50_GE) <- paste("MC", seq(ncol(Larger18_Normal_SC50_GE)), sep = "")

# Male_ASD
Male_ASD_SC10 <- SCimplify(as.matrix(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Male_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), Male_ASD_SC10$membership)
colnames(Male_ASD_SC10_GE) <- paste("MC", seq(ncol(Male_ASD_SC10_GE)), sep = "")

Male_ASD_SC20 <- supercell_rescale(Male_ASD_SC10, gamma = 20)
Male_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), Male_ASD_SC20$membership)
colnames(Male_ASD_SC20_GE) <- paste("MC", seq(ncol(Male_ASD_SC20_GE)), sep = "")

Male_ASD_SC30 <- supercell_rescale(Male_ASD_SC10, gamma = 30)
Male_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), Male_ASD_SC30$membership)
colnames(Male_ASD_SC30_GE) <- paste("MC", seq(ncol(Male_ASD_SC30_GE)), sep = "")

Male_ASD_SC40 <- supercell_rescale(Male_ASD_SC10, gamma = 40)
Male_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), Male_ASD_SC40$membership)
colnames(Male_ASD_SC40_GE) <- paste("MC", seq(ncol(Male_ASD_SC40_GE)), sep = "")

Male_ASD_SC50 <- supercell_rescale(Male_ASD_SC10, gamma = 50)
Male_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), Male_ASD_SC50$membership)
colnames(Male_ASD_SC50_GE) <- paste("MC", seq(ncol(Male_ASD_SC50_GE)), sep = "")

# Male_Normal
Male_Normal_SC10 <- SCimplify(as.matrix(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Male_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), Male_Normal_SC10$membership)
colnames(Male_Normal_SC10_GE) <- paste("MC", seq(ncol(Male_Normal_SC10_GE)), sep = "")

Male_Normal_SC20 <- supercell_rescale(Male_Normal_SC10, gamma = 20)
Male_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), Male_Normal_SC20$membership)
colnames(Male_Normal_SC20_GE) <- paste("MC", seq(ncol(Male_Normal_SC20_GE)), sep = "")

Male_Normal_SC30 <- supercell_rescale(Male_Normal_SC10, gamma = 30)
Male_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), Male_Normal_SC30$membership)
colnames(Male_Normal_SC30_GE) <- paste("MC", seq(ncol(Male_Normal_SC30_GE)), sep = "")

Male_Normal_SC40 <- supercell_rescale(Male_Normal_SC10, gamma = 40)
Male_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), Male_Normal_SC40$membership)
colnames(Male_Normal_SC40_GE) <- paste("MC", seq(ncol(Male_Normal_SC40_GE)), sep = "")

Male_Normal_SC50 <- supercell_rescale(Male_Normal_SC10, gamma = 50)
Male_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), Male_Normal_SC50$membership)
colnames(Male_Normal_SC50_GE) <- paste("MC", seq(ncol(Male_Normal_SC50_GE)), sep = "")

# Female_ASD
Female_ASD_SC10 <- SCimplify(as.matrix(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Female_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), Female_ASD_SC10$membership)
colnames(Female_ASD_SC10_GE) <- paste("MC", seq(ncol(Female_ASD_SC10_GE)), sep = "")

Female_ASD_SC20 <- supercell_rescale(Female_ASD_SC10, gamma = 20)
Female_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), Female_ASD_SC20$membership)
colnames(Female_ASD_SC20_GE) <- paste("MC", seq(ncol(Female_ASD_SC20_GE)), sep = "")

Female_ASD_SC30 <- supercell_rescale(Female_ASD_SC10, gamma = 30)
Female_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), Female_ASD_SC30$membership)
colnames(Female_ASD_SC30_GE) <- paste("MC", seq(ncol(Female_ASD_SC30_GE)), sep = "")

Female_ASD_SC40 <- supercell_rescale(Female_ASD_SC10, gamma = 40)
Female_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), Female_ASD_SC40$membership)
colnames(Female_ASD_SC40_GE) <- paste("MC", seq(ncol(Female_ASD_SC40_GE)), sep = "")

Female_ASD_SC50 <- supercell_rescale(Female_ASD_SC10, gamma = 50)
Female_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), Female_ASD_SC50$membership)
colnames(Female_ASD_SC50_GE) <- paste("MC", seq(ncol(Female_ASD_SC50_GE)), sep = "")

# Female_Normal
Female_Normal_SC10 <- SCimplify(as.matrix(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Female_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), Female_Normal_SC10$membership)
colnames(Female_Normal_SC10_GE) <- paste("MC", seq(ncol(Female_Normal_SC10_GE)), sep = "")

Female_Normal_SC20 <- supercell_rescale(Female_Normal_SC10, gamma = 20)
Female_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), Female_Normal_SC20$membership)
colnames(Female_Normal_SC20_GE) <- paste("MC", seq(ncol(Female_Normal_SC20_GE)), sep = "")

Female_Normal_SC30 <- supercell_rescale(Female_Normal_SC10, gamma = 30)
Female_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), Female_Normal_SC30$membership)
colnames(Female_Normal_SC30_GE) <- paste("MC", seq(ncol(Female_Normal_SC30_GE)), sep = "")

Female_Normal_SC40 <- supercell_rescale(Female_Normal_SC10, gamma = 40)
Female_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), Female_Normal_SC40$membership)
colnames(Female_Normal_SC40_GE) <- paste("MC", seq(ncol(Female_Normal_SC40_GE)), sep = "")

Female_Normal_SC50 <- supercell_rescale(Female_Normal_SC10, gamma = 50)
Female_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), Female_Normal_SC50$membership)
colnames(Female_Normal_SC50_GE) <- paste("MC", seq(ncol(Female_Normal_SC50_GE)), sep = "")

# NeuNRGNII_ASD
NeuNRGNII_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
NeuNRGNII_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), NeuNRGNII_ASD_SC10$membership)
colnames(NeuNRGNII_ASD_SC10_GE) <- paste("MC", seq(ncol(NeuNRGNII_ASD_SC10_GE)), sep = "")

NeuNRGNII_ASD_SC20 <- supercell_rescale(NeuNRGNII_ASD_SC10, gamma = 20)
NeuNRGNII_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), NeuNRGNII_ASD_SC20$membership)
colnames(NeuNRGNII_ASD_SC20_GE) <- paste("MC", seq(ncol(NeuNRGNII_ASD_SC20_GE)), sep = "")

NeuNRGNII_ASD_SC30 <- supercell_rescale(NeuNRGNII_ASD_SC10, gamma = 30)
NeuNRGNII_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), NeuNRGNII_ASD_SC30$membership)
colnames(NeuNRGNII_ASD_SC30_GE) <- paste("MC", seq(ncol(NeuNRGNII_ASD_SC30_GE)), sep = "")

NeuNRGNII_ASD_SC40 <- supercell_rescale(NeuNRGNII_ASD_SC10, gamma = 40)
NeuNRGNII_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), NeuNRGNII_ASD_SC40$membership)
colnames(NeuNRGNII_ASD_SC40_GE) <- paste("MC", seq(ncol(NeuNRGNII_ASD_SC40_GE)), sep = "")

NeuNRGNII_ASD_SC50 <- supercell_rescale(NeuNRGNII_ASD_SC10, gamma = 50)
NeuNRGNII_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), NeuNRGNII_ASD_SC50$membership)
colnames(NeuNRGNII_ASD_SC50_GE) <- paste("MC", seq(ncol(NeuNRGNII_ASD_SC50_GE)), sep = "")

# NeuNRGNII_Normal
NeuNRGNII_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
NeuNRGNII_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), NeuNRGNII_Normal_SC10$membership)
colnames(NeuNRGNII_Normal_SC10_GE) <- paste("MC", seq(ncol(NeuNRGNII_Normal_SC10_GE)), sep = "")

NeuNRGNII_Normal_SC20 <- supercell_rescale(NeuNRGNII_Normal_SC10, gamma = 20)
NeuNRGNII_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), NeuNRGNII_Normal_SC20$membership)
colnames(NeuNRGNII_Normal_SC20_GE) <- paste("MC", seq(ncol(NeuNRGNII_Normal_SC20_GE)), sep = "")

NeuNRGNII_Normal_SC30 <- supercell_rescale(NeuNRGNII_Normal_SC10, gamma = 30)
NeuNRGNII_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), NeuNRGNII_Normal_SC30$membership)
colnames(NeuNRGNII_Normal_SC30_GE) <- paste("MC", seq(ncol(NeuNRGNII_Normal_SC30_GE)), sep = "")

NeuNRGNII_Normal_SC40 <- supercell_rescale(NeuNRGNII_Normal_SC10, gamma = 40)
NeuNRGNII_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), NeuNRGNII_Normal_SC40$membership)
colnames(NeuNRGNII_Normal_SC40_GE) <- paste("MC", seq(ncol(NeuNRGNII_Normal_SC40_GE)), sep = "")

NeuNRGNII_Normal_SC50 <- supercell_rescale(NeuNRGNII_Normal_SC10, gamma = 50)
NeuNRGNII_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), NeuNRGNII_Normal_SC50$membership)
colnames(NeuNRGNII_Normal_SC50_GE) <- paste("MC", seq(ncol(NeuNRGNII_Normal_SC50_GE)), sep = "")

# L56_ASD
L56_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L56_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), L56_ASD_SC10$membership)
colnames(L56_ASD_SC10_GE) <- paste("MC", seq(ncol(L56_ASD_SC10_GE)), sep = "")

L56_ASD_SC20 <- supercell_rescale(L56_ASD_SC10, gamma = 20)
L56_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), L56_ASD_SC20$membership)
colnames(L56_ASD_SC20_GE) <- paste("MC", seq(ncol(L56_ASD_SC20_GE)), sep = "")

L56_ASD_SC30 <- supercell_rescale(L56_ASD_SC10, gamma = 30)
L56_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), L56_ASD_SC30$membership)
colnames(L56_ASD_SC30_GE) <- paste("MC", seq(ncol(L56_ASD_SC30_GE)), sep = "")

L56_ASD_SC40 <- supercell_rescale(L56_ASD_SC10, gamma = 40)
L56_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), L56_ASD_SC40$membership)
colnames(L56_ASD_SC40_GE) <- paste("MC", seq(ncol(L56_ASD_SC40_GE)), sep = "")

L56_ASD_SC50 <- supercell_rescale(L56_ASD_SC10, gamma = 50)
L56_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), L56_ASD_SC50$membership)
colnames(L56_ASD_SC50_GE) <- paste("MC", seq(ncol(L56_ASD_SC50_GE)), sep = "")

# L56_Normal
L56_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L56_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), L56_Normal_SC10$membership)
colnames(L56_Normal_SC10_GE) <- paste("MC", seq(ncol(L56_Normal_SC10_GE)), sep = "")

L56_Normal_SC20 <- supercell_rescale(L56_Normal_SC10, gamma = 20)
L56_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), L56_Normal_SC20$membership)
colnames(L56_Normal_SC20_GE) <- paste("MC", seq(ncol(L56_Normal_SC20_GE)), sep = "")

L56_Normal_SC30 <- supercell_rescale(L56_Normal_SC10, gamma = 30)
L56_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), L56_Normal_SC30$membership)
colnames(L56_Normal_SC30_GE) <- paste("MC", seq(ncol(L56_Normal_SC30_GE)), sep = "")

L56_Normal_SC40 <- supercell_rescale(L56_Normal_SC10, gamma = 40)
L56_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), L56_Normal_SC40$membership)
colnames(L56_Normal_SC40_GE) <- paste("MC", seq(ncol(L56_Normal_SC40_GE)), sep = "")

L56_Normal_SC50 <- supercell_rescale(L56_Normal_SC10, gamma = 50)
L56_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), L56_Normal_SC50$membership)
colnames(L56_Normal_SC50_GE) <- paste("MC", seq(ncol(L56_Normal_SC50_GE)), sep = "")

# Oligodendrocytes_ASD
Oligodendrocytes_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Oligodendrocytes_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), Oligodendrocytes_ASD_SC10$membership)
colnames(Oligodendrocytes_ASD_SC10_GE) <- paste("MC", seq(ncol(Oligodendrocytes_ASD_SC10_GE)), sep = "")

Oligodendrocytes_ASD_SC20 <- supercell_rescale(Oligodendrocytes_ASD_SC10, gamma = 20)
Oligodendrocytes_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), Oligodendrocytes_ASD_SC20$membership)
colnames(Oligodendrocytes_ASD_SC20_GE) <- paste("MC", seq(ncol(Oligodendrocytes_ASD_SC20_GE)), sep = "")

Oligodendrocytes_ASD_SC30 <- supercell_rescale(Oligodendrocytes_ASD_SC10, gamma = 30)
Oligodendrocytes_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), Oligodendrocytes_ASD_SC30$membership)
colnames(Oligodendrocytes_ASD_SC30_GE) <- paste("MC", seq(ncol(Oligodendrocytes_ASD_SC30_GE)), sep = "")

Oligodendrocytes_ASD_SC40 <- supercell_rescale(Oligodendrocytes_ASD_SC10, gamma = 40)
Oligodendrocytes_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), Oligodendrocytes_ASD_SC40$membership)
colnames(Oligodendrocytes_ASD_SC40_GE) <- paste("MC", seq(ncol(Oligodendrocytes_ASD_SC40_GE)), sep = "")

Oligodendrocytes_ASD_SC50 <- supercell_rescale(Oligodendrocytes_ASD_SC10, gamma = 50)
Oligodendrocytes_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), Oligodendrocytes_ASD_SC50$membership)
colnames(Oligodendrocytes_ASD_SC50_GE) <- paste("MC", seq(ncol(Oligodendrocytes_ASD_SC50_GE)), sep = "")

# Oligodendrocytes_Normal
Oligodendrocytes_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Oligodendrocytes_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), Oligodendrocytes_Normal_SC10$membership)
colnames(Oligodendrocytes_Normal_SC10_GE) <- paste("MC", seq(ncol(Oligodendrocytes_Normal_SC10_GE)), sep = "")

Oligodendrocytes_Normal_SC20 <- supercell_rescale(Oligodendrocytes_Normal_SC10, gamma = 20)
Oligodendrocytes_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), Oligodendrocytes_Normal_SC20$membership)
colnames(Oligodendrocytes_Normal_SC20_GE) <- paste("MC", seq(ncol(Oligodendrocytes_Normal_SC20_GE)), sep = "")

Oligodendrocytes_Normal_SC30 <- supercell_rescale(Oligodendrocytes_Normal_SC10, gamma = 30)
Oligodendrocytes_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), Oligodendrocytes_Normal_SC30$membership)
colnames(Oligodendrocytes_Normal_SC30_GE) <- paste("MC", seq(ncol(Oligodendrocytes_Normal_SC30_GE)), sep = "")

Oligodendrocytes_Normal_SC40 <- supercell_rescale(Oligodendrocytes_Normal_SC10, gamma = 40)
Oligodendrocytes_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), Oligodendrocytes_Normal_SC40$membership)
colnames(Oligodendrocytes_Normal_SC40_GE) <- paste("MC", seq(ncol(Oligodendrocytes_Normal_SC40_GE)), sep = "")

Oligodendrocytes_Normal_SC50 <- supercell_rescale(Oligodendrocytes_Normal_SC10, gamma = 50)
Oligodendrocytes_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), Oligodendrocytes_Normal_SC50$membership)
colnames(Oligodendrocytes_Normal_SC50_GE) <- paste("MC", seq(ncol(Oligodendrocytes_Normal_SC50_GE)), sep = "")

# OPC_ASD
OPC_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
OPC_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), OPC_ASD_SC10$membership)
colnames(OPC_ASD_SC10_GE) <- paste("MC", seq(ncol(OPC_ASD_SC10_GE)), sep = "")

OPC_ASD_SC20 <- supercell_rescale(OPC_ASD_SC10, gamma = 20)
OPC_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), OPC_ASD_SC20$membership)
colnames(OPC_ASD_SC20_GE) <- paste("MC", seq(ncol(OPC_ASD_SC20_GE)), sep = "")

OPC_ASD_SC30 <- supercell_rescale(OPC_ASD_SC10, gamma = 30)
OPC_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), OPC_ASD_SC30$membership)
colnames(OPC_ASD_SC30_GE) <- paste("MC", seq(ncol(OPC_ASD_SC30_GE)), sep = "")

OPC_ASD_SC40 <- supercell_rescale(OPC_ASD_SC10, gamma = 40)
OPC_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), OPC_ASD_SC40$membership)
colnames(OPC_ASD_SC40_GE) <- paste("MC", seq(ncol(OPC_ASD_SC40_GE)), sep = "")

OPC_ASD_SC50 <- supercell_rescale(OPC_ASD_SC10, gamma = 50)
OPC_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), OPC_ASD_SC50$membership)
colnames(OPC_ASD_SC50_GE) <- paste("MC", seq(ncol(OPC_ASD_SC50_GE)), sep = "")

# OPC_Normal
OPC_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
OPC_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), OPC_Normal_SC10$membership)
colnames(OPC_Normal_SC10_GE) <- paste("MC", seq(ncol(OPC_Normal_SC10_GE)), sep = "")

OPC_Normal_SC20 <- supercell_rescale(OPC_Normal_SC10, gamma = 20)
OPC_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), OPC_Normal_SC20$membership)
colnames(OPC_Normal_SC20_GE) <- paste("MC", seq(ncol(OPC_Normal_SC20_GE)), sep = "")

OPC_Normal_SC30 <- supercell_rescale(OPC_Normal_SC10, gamma = 30)
OPC_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), OPC_Normal_SC30$membership)
colnames(OPC_Normal_SC30_GE) <- paste("MC", seq(ncol(OPC_Normal_SC30_GE)), sep = "")

OPC_Normal_SC40 <- supercell_rescale(OPC_Normal_SC10, gamma = 40)
OPC_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), OPC_Normal_SC40$membership)
colnames(OPC_Normal_SC40_GE) <- paste("MC", seq(ncol(OPC_Normal_SC40_GE)), sep = "")

OPC_Normal_SC50 <- supercell_rescale(OPC_Normal_SC10, gamma = 50)
OPC_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), OPC_Normal_SC50$membership)
colnames(OPC_Normal_SC50_GE) <- paste("MC", seq(ncol(OPC_Normal_SC50_GE)), sep = "")

# ASTFB_ASD
ASTFB_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ASTFB_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), ASTFB_ASD_SC10$membership)
colnames(ASTFB_ASD_SC10_GE) <- paste("MC", seq(ncol(ASTFB_ASD_SC10_GE)), sep = "")

ASTFB_ASD_SC20 <- supercell_rescale(ASTFB_ASD_SC10, gamma = 20)
ASTFB_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), ASTFB_ASD_SC20$membership)
colnames(ASTFB_ASD_SC20_GE) <- paste("MC", seq(ncol(ASTFB_ASD_SC20_GE)), sep = "")

ASTFB_ASD_SC30 <- supercell_rescale(ASTFB_ASD_SC10, gamma = 30)
ASTFB_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), ASTFB_ASD_SC30$membership)
colnames(ASTFB_ASD_SC30_GE) <- paste("MC", seq(ncol(ASTFB_ASD_SC30_GE)), sep = "")

ASTFB_ASD_SC40 <- supercell_rescale(ASTFB_ASD_SC10, gamma = 40)
ASTFB_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), ASTFB_ASD_SC40$membership)
colnames(ASTFB_ASD_SC40_GE) <- paste("MC", seq(ncol(ASTFB_ASD_SC40_GE)), sep = "")

ASTFB_ASD_SC50 <- supercell_rescale(ASTFB_ASD_SC10, gamma = 50)
ASTFB_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), ASTFB_ASD_SC50$membership)
colnames(ASTFB_ASD_SC50_GE) <- paste("MC", seq(ncol(ASTFB_ASD_SC50_GE)), sep = "")

# ASTFB_Normal
ASTFB_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ASTFB_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), ASTFB_Normal_SC10$membership)
colnames(ASTFB_Normal_SC10_GE) <- paste("MC", seq(ncol(ASTFB_Normal_SC10_GE)), sep = "")

ASTFB_Normal_SC20 <- supercell_rescale(ASTFB_Normal_SC10, gamma = 20)
ASTFB_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), ASTFB_Normal_SC20$membership)
colnames(ASTFB_Normal_SC20_GE) <- paste("MC", seq(ncol(ASTFB_Normal_SC20_GE)), sep = "")

ASTFB_Normal_SC30 <- supercell_rescale(ASTFB_Normal_SC10, gamma = 30)
ASTFB_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), ASTFB_Normal_SC30$membership)
colnames(ASTFB_Normal_SC30_GE) <- paste("MC", seq(ncol(ASTFB_Normal_SC30_GE)), sep = "")

ASTFB_Normal_SC40 <- supercell_rescale(ASTFB_Normal_SC10, gamma = 40)
ASTFB_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), ASTFB_Normal_SC40$membership)
colnames(ASTFB_Normal_SC40_GE) <- paste("MC", seq(ncol(ASTFB_Normal_SC40_GE)), sep = "")

ASTFB_Normal_SC50 <- supercell_rescale(ASTFB_Normal_SC10, gamma = 50)
ASTFB_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), ASTFB_Normal_SC50$membership)
colnames(ASTFB_Normal_SC50_GE) <- paste("MC", seq(ncol(ASTFB_Normal_SC50_GE)), sep = "")

# Endothelial_ASD
Endothelial_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Endothelial_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), Endothelial_ASD_SC10$membership)
colnames(Endothelial_ASD_SC10_GE) <- paste("MC", seq(ncol(Endothelial_ASD_SC10_GE)), sep = "")

Endothelial_ASD_SC20 <- supercell_rescale(Endothelial_ASD_SC10, gamma = 20)
Endothelial_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), Endothelial_ASD_SC20$membership)
colnames(Endothelial_ASD_SC20_GE) <- paste("MC", seq(ncol(Endothelial_ASD_SC20_GE)), sep = "")

Endothelial_ASD_SC30 <- supercell_rescale(Endothelial_ASD_SC10, gamma = 30)
Endothelial_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), Endothelial_ASD_SC30$membership)
colnames(Endothelial_ASD_SC30_GE) <- paste("MC", seq(ncol(Endothelial_ASD_SC30_GE)), sep = "")

Endothelial_ASD_SC40 <- supercell_rescale(Endothelial_ASD_SC10, gamma = 40)
Endothelial_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), Endothelial_ASD_SC40$membership)
colnames(Endothelial_ASD_SC40_GE) <- paste("MC", seq(ncol(Endothelial_ASD_SC40_GE)), sep = "")

Endothelial_ASD_SC50 <- supercell_rescale(Endothelial_ASD_SC10, gamma = 50)
Endothelial_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), Endothelial_ASD_SC50$membership)
colnames(Endothelial_ASD_SC50_GE) <- paste("MC", seq(ncol(Endothelial_ASD_SC50_GE)), sep = "")

# Endothelial_Normal
Endothelial_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Endothelial_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), Endothelial_Normal_SC10$membership)
colnames(Endothelial_Normal_SC10_GE) <- paste("MC", seq(ncol(Endothelial_Normal_SC10_GE)), sep = "")

Endothelial_Normal_SC20 <- supercell_rescale(Endothelial_Normal_SC10, gamma = 20)
Endothelial_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), Endothelial_Normal_SC20$membership)
colnames(Endothelial_Normal_SC20_GE) <- paste("MC", seq(ncol(Endothelial_Normal_SC20_GE)), sep = "")

Endothelial_Normal_SC30 <- supercell_rescale(Endothelial_Normal_SC10, gamma = 30)
Endothelial_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), Endothelial_Normal_SC30$membership)
colnames(Endothelial_Normal_SC30_GE) <- paste("MC", seq(ncol(Endothelial_Normal_SC30_GE)), sep = "")

Endothelial_Normal_SC40 <- supercell_rescale(Endothelial_Normal_SC10, gamma = 40)
Endothelial_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), Endothelial_Normal_SC40$membership)
colnames(Endothelial_Normal_SC40_GE) <- paste("MC", seq(ncol(Endothelial_Normal_SC40_GE)), sep = "")

Endothelial_Normal_SC50 <- supercell_rescale(Endothelial_Normal_SC10, gamma = 50)
Endothelial_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), Endothelial_Normal_SC50$membership)
colnames(Endothelial_Normal_SC50_GE) <- paste("MC", seq(ncol(Endothelial_Normal_SC50_GE)), sep = "")

# Microglia_ASD
Microglia_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Microglia_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), Microglia_ASD_SC10$membership)
colnames(Microglia_ASD_SC10_GE) <- paste("MC", seq(ncol(Microglia_ASD_SC10_GE)), sep = "")

Microglia_ASD_SC20 <- supercell_rescale(Microglia_ASD_SC10, gamma = 20)
Microglia_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), Microglia_ASD_SC20$membership)
colnames(Microglia_ASD_SC20_GE) <- paste("MC", seq(ncol(Microglia_ASD_SC20_GE)), sep = "")

Microglia_ASD_SC30 <- supercell_rescale(Microglia_ASD_SC10, gamma = 30)
Microglia_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), Microglia_ASD_SC30$membership)
colnames(Microglia_ASD_SC30_GE) <- paste("MC", seq(ncol(Microglia_ASD_SC30_GE)), sep = "")

Microglia_ASD_SC40 <- supercell_rescale(Microglia_ASD_SC10, gamma = 40)
Microglia_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), Microglia_ASD_SC40$membership)
colnames(Microglia_ASD_SC40_GE) <- paste("MC", seq(ncol(Microglia_ASD_SC40_GE)), sep = "")

Microglia_ASD_SC50 <- supercell_rescale(Microglia_ASD_SC10, gamma = 50)
Microglia_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), Microglia_ASD_SC50$membership)
colnames(Microglia_ASD_SC50_GE) <- paste("MC", seq(ncol(Microglia_ASD_SC50_GE)), sep = "")

# Microglia_Normal
Microglia_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Microglia_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), Microglia_Normal_SC10$membership)
colnames(Microglia_Normal_SC10_GE) <- paste("MC", seq(ncol(Microglia_Normal_SC10_GE)), sep = "")

Microglia_Normal_SC20 <- supercell_rescale(Microglia_Normal_SC10, gamma = 20)
Microglia_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), Microglia_Normal_SC20$membership)
colnames(Microglia_Normal_SC20_GE) <- paste("MC", seq(ncol(Microglia_Normal_SC20_GE)), sep = "")

Microglia_Normal_SC30 <- supercell_rescale(Microglia_Normal_SC10, gamma = 30)
Microglia_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), Microglia_Normal_SC30$membership)
colnames(Microglia_Normal_SC30_GE) <- paste("MC", seq(ncol(Microglia_Normal_SC30_GE)), sep = "")

Microglia_Normal_SC40 <- supercell_rescale(Microglia_Normal_SC10, gamma = 40)
Microglia_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), Microglia_Normal_SC40$membership)
colnames(Microglia_Normal_SC40_GE) <- paste("MC", seq(ncol(Microglia_Normal_SC40_GE)), sep = "")

Microglia_Normal_SC50 <- supercell_rescale(Microglia_Normal_SC10, gamma = 50)
Microglia_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), Microglia_Normal_SC50$membership)
colnames(Microglia_Normal_SC50_GE) <- paste("MC", seq(ncol(Microglia_Normal_SC50_GE)), sep = "")

# NeuNRGNI_ASD
NeuNRGNI_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
NeuNRGNI_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), NeuNRGNI_ASD_SC10$membership)
colnames(NeuNRGNI_ASD_SC10_GE) <- paste("MC", seq(ncol(NeuNRGNI_ASD_SC10_GE)), sep = "")

NeuNRGNI_ASD_SC20 <- supercell_rescale(NeuNRGNI_ASD_SC10, gamma = 20)
NeuNRGNI_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), NeuNRGNI_ASD_SC20$membership)
colnames(NeuNRGNI_ASD_SC20_GE) <- paste("MC", seq(ncol(NeuNRGNI_ASD_SC20_GE)), sep = "")

NeuNRGNI_ASD_SC30 <- supercell_rescale(NeuNRGNI_ASD_SC10, gamma = 30)
NeuNRGNI_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), NeuNRGNI_ASD_SC30$membership)
colnames(NeuNRGNI_ASD_SC30_GE) <- paste("MC", seq(ncol(NeuNRGNI_ASD_SC30_GE)), sep = "")

NeuNRGNI_ASD_SC40 <- supercell_rescale(NeuNRGNI_ASD_SC10, gamma = 40)
NeuNRGNI_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), NeuNRGNI_ASD_SC40$membership)
colnames(NeuNRGNI_ASD_SC40_GE) <- paste("MC", seq(ncol(NeuNRGNI_ASD_SC40_GE)), sep = "")

NeuNRGNI_ASD_SC50 <- supercell_rescale(NeuNRGNI_ASD_SC10, gamma = 50)
NeuNRGNI_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), NeuNRGNI_ASD_SC50$membership)
colnames(NeuNRGNI_ASD_SC50_GE) <- paste("MC", seq(ncol(NeuNRGNI_ASD_SC50_GE)), sep = "")

# NeuNRGNI_Normal
NeuNRGNI_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
NeuNRGNI_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), NeuNRGNI_Normal_SC10$membership)
colnames(NeuNRGNI_Normal_SC10_GE) <- paste("MC", seq(ncol(NeuNRGNI_Normal_SC10_GE)), sep = "")

NeuNRGNI_Normal_SC20 <- supercell_rescale(NeuNRGNI_Normal_SC10, gamma = 20)
NeuNRGNI_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), NeuNRGNI_Normal_SC20$membership)
colnames(NeuNRGNI_Normal_SC20_GE) <- paste("MC", seq(ncol(NeuNRGNI_Normal_SC20_GE)), sep = "")

NeuNRGNI_Normal_SC30 <- supercell_rescale(NeuNRGNI_Normal_SC10, gamma = 30)
NeuNRGNI_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), NeuNRGNI_Normal_SC30$membership)
colnames(NeuNRGNI_Normal_SC30_GE) <- paste("MC", seq(ncol(NeuNRGNI_Normal_SC30_GE)), sep = "")

NeuNRGNI_Normal_SC40 <- supercell_rescale(NeuNRGNI_Normal_SC10, gamma = 40)
NeuNRGNI_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), NeuNRGNI_Normal_SC40$membership)
colnames(NeuNRGNI_Normal_SC40_GE) <- paste("MC", seq(ncol(NeuNRGNI_Normal_SC40_GE)), sep = "")

NeuNRGNI_Normal_SC50 <- supercell_rescale(NeuNRGNI_Normal_SC10, gamma = 50)
NeuNRGNI_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), NeuNRGNI_Normal_SC50$membership)
colnames(NeuNRGNI_Normal_SC50_GE) <- paste("MC", seq(ncol(NeuNRGNI_Normal_SC50_GE)), sep = "")

# INVIP_ASD
INVIP_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INVIP_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), INVIP_ASD_SC10$membership)
colnames(INVIP_ASD_SC10_GE) <- paste("MC", seq(ncol(INVIP_ASD_SC10_GE)), sep = "")

INVIP_ASD_SC20 <- supercell_rescale(INVIP_ASD_SC10, gamma = 20)
INVIP_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), INVIP_ASD_SC20$membership)
colnames(INVIP_ASD_SC20_GE) <- paste("MC", seq(ncol(INVIP_ASD_SC20_GE)), sep = "")

INVIP_ASD_SC30 <- supercell_rescale(INVIP_ASD_SC10, gamma = 30)
INVIP_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), INVIP_ASD_SC30$membership)
colnames(INVIP_ASD_SC30_GE) <- paste("MC", seq(ncol(INVIP_ASD_SC30_GE)), sep = "")

INVIP_ASD_SC40 <- supercell_rescale(INVIP_ASD_SC10, gamma = 40)
INVIP_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), INVIP_ASD_SC40$membership)
colnames(INVIP_ASD_SC40_GE) <- paste("MC", seq(ncol(INVIP_ASD_SC40_GE)), sep = "")

INVIP_ASD_SC50 <- supercell_rescale(INVIP_ASD_SC10, gamma = 50)
INVIP_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), INVIP_ASD_SC50$membership)
colnames(INVIP_ASD_SC50_GE) <- paste("MC", seq(ncol(INVIP_ASD_SC50_GE)), sep = "")

# INVIP_Normal
INVIP_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INVIP_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), INVIP_Normal_SC10$membership)
colnames(INVIP_Normal_SC10_GE) <- paste("MC", seq(ncol(INVIP_Normal_SC10_GE)), sep = "")

INVIP_Normal_SC20 <- supercell_rescale(INVIP_Normal_SC10, gamma = 20)
INVIP_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), INVIP_Normal_SC20$membership)
colnames(INVIP_Normal_SC20_GE) <- paste("MC", seq(ncol(INVIP_Normal_SC20_GE)), sep = "")

INVIP_Normal_SC30 <- supercell_rescale(INVIP_Normal_SC10, gamma = 30)
INVIP_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), INVIP_Normal_SC30$membership)
colnames(INVIP_Normal_SC30_GE) <- paste("MC", seq(ncol(INVIP_Normal_SC30_GE)), sep = "")

INVIP_Normal_SC40 <- supercell_rescale(INVIP_Normal_SC10, gamma = 40)
INVIP_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), INVIP_Normal_SC40$membership)
colnames(INVIP_Normal_SC40_GE) <- paste("MC", seq(ncol(INVIP_Normal_SC40_GE)), sep = "")

INVIP_Normal_SC50 <- supercell_rescale(INVIP_Normal_SC10, gamma = 50)
INVIP_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), INVIP_Normal_SC50$membership)
colnames(INVIP_Normal_SC50_GE) <- paste("MC", seq(ncol(INVIP_Normal_SC50_GE)), sep = "")

# L56CC_ASD
L56CC_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L56CC_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), L56CC_ASD_SC10$membership)
colnames(L56CC_ASD_SC10_GE) <- paste("MC", seq(ncol(L56CC_ASD_SC10_GE)), sep = "")

L56CC_ASD_SC20 <- supercell_rescale(L56CC_ASD_SC10, gamma = 20)
L56CC_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), L56CC_ASD_SC20$membership)
colnames(L56CC_ASD_SC20_GE) <- paste("MC", seq(ncol(L56CC_ASD_SC20_GE)), sep = "")

L56CC_ASD_SC30 <- supercell_rescale(L56CC_ASD_SC10, gamma = 30)
L56CC_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), L56CC_ASD_SC30$membership)
colnames(L56CC_ASD_SC30_GE) <- paste("MC", seq(ncol(L56CC_ASD_SC30_GE)), sep = "")

L56CC_ASD_SC40 <- supercell_rescale(L56CC_ASD_SC10, gamma = 40)
L56CC_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), L56CC_ASD_SC40$membership)
colnames(L56CC_ASD_SC40_GE) <- paste("MC", seq(ncol(L56CC_ASD_SC40_GE)), sep = "")

L56CC_ASD_SC50 <- supercell_rescale(L56CC_ASD_SC10, gamma = 50)
L56CC_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), L56CC_ASD_SC50$membership)
colnames(L56CC_ASD_SC50_GE) <- paste("MC", seq(ncol(L56CC_ASD_SC50_GE)), sep = "")

# L56CC_Normal
L56CC_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L56CC_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), L56CC_Normal_SC10$membership)
colnames(L56CC_Normal_SC10_GE) <- paste("MC", seq(ncol(L56CC_Normal_SC10_GE)), sep = "")

L56CC_Normal_SC20 <- supercell_rescale(L56CC_Normal_SC10, gamma = 20)
L56CC_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), L56CC_Normal_SC20$membership)
colnames(L56CC_Normal_SC20_GE) <- paste("MC", seq(ncol(L56CC_Normal_SC20_GE)), sep = "")

L56CC_Normal_SC30 <- supercell_rescale(L56CC_Normal_SC10, gamma = 30)
L56CC_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), L56CC_Normal_SC30$membership)
colnames(L56CC_Normal_SC30_GE) <- paste("MC", seq(ncol(L56CC_Normal_SC30_GE)), sep = "")

L56CC_Normal_SC40 <- supercell_rescale(L56CC_Normal_SC10, gamma = 40)
L56CC_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), L56CC_Normal_SC40$membership)
colnames(L56CC_Normal_SC40_GE) <- paste("MC", seq(ncol(L56CC_Normal_SC40_GE)), sep = "")

L56CC_Normal_SC50 <- supercell_rescale(L56CC_Normal_SC10, gamma = 50)
L56CC_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), L56CC_Normal_SC50$membership)
colnames(L56CC_Normal_SC50_GE) <- paste("MC", seq(ncol(L56CC_Normal_SC50_GE)), sep = "")

# INSV2C_ASD
INSV2C_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INSV2C_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), INSV2C_ASD_SC10$membership)
colnames(INSV2C_ASD_SC10_GE) <- paste("MC", seq(ncol(INSV2C_ASD_SC10_GE)), sep = "")

INSV2C_ASD_SC20 <- supercell_rescale(INSV2C_ASD_SC10, gamma = 20)
INSV2C_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), INSV2C_ASD_SC20$membership)
colnames(INSV2C_ASD_SC20_GE) <- paste("MC", seq(ncol(INSV2C_ASD_SC20_GE)), sep = "")

INSV2C_ASD_SC30 <- supercell_rescale(INSV2C_ASD_SC10, gamma = 30)
INSV2C_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), INSV2C_ASD_SC30$membership)
colnames(INSV2C_ASD_SC30_GE) <- paste("MC", seq(ncol(INSV2C_ASD_SC30_GE)), sep = "")

INSV2C_ASD_SC40 <- supercell_rescale(INSV2C_ASD_SC10, gamma = 40)
INSV2C_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), INSV2C_ASD_SC40$membership)
colnames(INSV2C_ASD_SC40_GE) <- paste("MC", seq(ncol(INSV2C_ASD_SC40_GE)), sep = "")

INSV2C_ASD_SC50 <- supercell_rescale(INSV2C_ASD_SC10, gamma = 50)
INSV2C_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), INSV2C_ASD_SC50$membership)
colnames(INSV2C_ASD_SC50_GE) <- paste("MC", seq(ncol(INSV2C_ASD_SC50_GE)), sep = "")

# INSV2C_Normal
INSV2C_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INSV2C_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), INSV2C_Normal_SC10$membership)
colnames(INSV2C_Normal_SC10_GE) <- paste("MC", seq(ncol(INSV2C_Normal_SC10_GE)), sep = "")

INSV2C_Normal_SC20 <- supercell_rescale(INSV2C_Normal_SC10, gamma = 20)
INSV2C_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), INSV2C_Normal_SC20$membership)
colnames(INSV2C_Normal_SC20_GE) <- paste("MC", seq(ncol(INSV2C_Normal_SC20_GE)), sep = "")

INSV2C_Normal_SC30 <- supercell_rescale(INSV2C_Normal_SC10, gamma = 30)
INSV2C_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), INSV2C_Normal_SC30$membership)
colnames(INSV2C_Normal_SC30_GE) <- paste("MC", seq(ncol(INSV2C_Normal_SC30_GE)), sep = "")

INSV2C_Normal_SC40 <- supercell_rescale(INSV2C_Normal_SC10, gamma = 40)
INSV2C_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), INSV2C_Normal_SC40$membership)
colnames(INSV2C_Normal_SC40_GE) <- paste("MC", seq(ncol(INSV2C_Normal_SC40_GE)), sep = "")

INSV2C_Normal_SC50 <- supercell_rescale(INSV2C_Normal_SC10, gamma = 50)
INSV2C_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), INSV2C_Normal_SC50$membership)
colnames(INSV2C_Normal_SC50_GE) <- paste("MC", seq(ncol(INSV2C_Normal_SC50_GE)), sep = "")

# L23_ASD
L23_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L23_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), L23_ASD_SC10$membership)
colnames(L23_ASD_SC10_GE) <- paste("MC", seq(ncol(L23_ASD_SC10_GE)), sep = "")

L23_ASD_SC20 <- supercell_rescale(L23_ASD_SC10, gamma = 20)
L23_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), L23_ASD_SC20$membership)
colnames(L23_ASD_SC20_GE) <- paste("MC", seq(ncol(L23_ASD_SC20_GE)), sep = "")

L23_ASD_SC30 <- supercell_rescale(L23_ASD_SC10, gamma = 30)
L23_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), L23_ASD_SC30$membership)
colnames(L23_ASD_SC30_GE) <- paste("MC", seq(ncol(L23_ASD_SC30_GE)), sep = "")

L23_ASD_SC40 <- supercell_rescale(L23_ASD_SC10, gamma = 40)
L23_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), L23_ASD_SC40$membership)
colnames(L23_ASD_SC40_GE) <- paste("MC", seq(ncol(L23_ASD_SC40_GE)), sep = "")

L23_ASD_SC50 <- supercell_rescale(L23_ASD_SC10, gamma = 50)
L23_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), L23_ASD_SC50$membership)
colnames(L23_ASD_SC50_GE) <- paste("MC", seq(ncol(L23_ASD_SC50_GE)), sep = "")

# L23_Normal
L23_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L23_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), L23_Normal_SC10$membership)
colnames(L23_Normal_SC10_GE) <- paste("MC", seq(ncol(L23_Normal_SC10_GE)), sep = "")

L23_Normal_SC20 <- supercell_rescale(L23_Normal_SC10, gamma = 20)
L23_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), L23_Normal_SC20$membership)
colnames(L23_Normal_SC20_GE) <- paste("MC", seq(ncol(L23_Normal_SC20_GE)), sep = "")

L23_Normal_SC30 <- supercell_rescale(L23_Normal_SC10, gamma = 30)
L23_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), L23_Normal_SC30$membership)
colnames(L23_Normal_SC30_GE) <- paste("MC", seq(ncol(L23_Normal_SC30_GE)), sep = "")

L23_Normal_SC40 <- supercell_rescale(L23_Normal_SC10, gamma = 40)
L23_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), L23_Normal_SC40$membership)
colnames(L23_Normal_SC40_GE) <- paste("MC", seq(ncol(L23_Normal_SC40_GE)), sep = "")

L23_Normal_SC50 <- supercell_rescale(L23_Normal_SC10, gamma = 50)
L23_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), L23_Normal_SC50$membership)
colnames(L23_Normal_SC50_GE) <- paste("MC", seq(ncol(L23_Normal_SC50_GE)), sep = "")

# INPV_ASD
INPV_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INPV_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), INPV_ASD_SC10$membership)
colnames(INPV_ASD_SC10_GE) <- paste("MC", seq(ncol(INPV_ASD_SC10_GE)), sep = "")

INPV_ASD_SC20 <- supercell_rescale(INPV_ASD_SC10, gamma = 20)
INPV_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), INPV_ASD_SC20$membership)
colnames(INPV_ASD_SC20_GE) <- paste("MC", seq(ncol(INPV_ASD_SC20_GE)), sep = "")

INPV_ASD_SC30 <- supercell_rescale(INPV_ASD_SC10, gamma = 30)
INPV_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), INPV_ASD_SC30$membership)
colnames(INPV_ASD_SC30_GE) <- paste("MC", seq(ncol(INPV_ASD_SC30_GE)), sep = "")

INPV_ASD_SC40 <- supercell_rescale(INPV_ASD_SC10, gamma = 40)
INPV_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), INPV_ASD_SC40$membership)
colnames(INPV_ASD_SC40_GE) <- paste("MC", seq(ncol(INPV_ASD_SC40_GE)), sep = "")

INPV_ASD_SC50 <- supercell_rescale(INPV_ASD_SC10, gamma = 50)
INPV_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), INPV_ASD_SC50$membership)
colnames(INPV_ASD_SC50_GE) <- paste("MC", seq(ncol(INPV_ASD_SC50_GE)), sep = "")

# INPV_Normal
INPV_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INPV_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), INPV_Normal_SC10$membership)
colnames(INPV_Normal_SC10_GE) <- paste("MC", seq(ncol(INPV_Normal_SC10_GE)), sep = "")

INPV_Normal_SC20 <- supercell_rescale(INPV_Normal_SC10, gamma = 20)
INPV_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), INPV_Normal_SC20$membership)
colnames(INPV_Normal_SC20_GE) <- paste("MC", seq(ncol(INPV_Normal_SC20_GE)), sep = "")

INPV_Normal_SC30 <- supercell_rescale(INPV_Normal_SC10, gamma = 30)
INPV_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), INPV_Normal_SC30$membership)
colnames(INPV_Normal_SC30_GE) <- paste("MC", seq(ncol(INPV_Normal_SC30_GE)), sep = "")

INPV_Normal_SC40 <- supercell_rescale(INPV_Normal_SC10, gamma = 40)
INPV_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), INPV_Normal_SC40$membership)
colnames(INPV_Normal_SC40_GE) <- paste("MC", seq(ncol(INPV_Normal_SC40_GE)), sep = "")

INPV_Normal_SC50 <- supercell_rescale(INPV_Normal_SC10, gamma = 50)
INPV_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), INPV_Normal_SC50$membership)
colnames(INPV_Normal_SC50_GE) <- paste("MC", seq(ncol(INPV_Normal_SC50_GE)), sep = "")

# L4_ASD
L4_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L4_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), L4_ASD_SC10$membership)
colnames(L4_ASD_SC10_GE) <- paste("MC", seq(ncol(L4_ASD_SC10_GE)), sep = "")

L4_ASD_SC20 <- supercell_rescale(L4_ASD_SC10, gamma = 20)
L4_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), L4_ASD_SC20$membership)
colnames(L4_ASD_SC20_GE) <- paste("MC", seq(ncol(L4_ASD_SC20_GE)), sep = "")

L4_ASD_SC30 <- supercell_rescale(L4_ASD_SC10, gamma = 30)
L4_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), L4_ASD_SC30$membership)
colnames(L4_ASD_SC30_GE) <- paste("MC", seq(ncol(L4_ASD_SC30_GE)), sep = "")

L4_ASD_SC40 <- supercell_rescale(L4_ASD_SC10, gamma = 40)
L4_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), L4_ASD_SC40$membership)
colnames(L4_ASD_SC40_GE) <- paste("MC", seq(ncol(L4_ASD_SC40_GE)), sep = "")

L4_ASD_SC50 <- supercell_rescale(L4_ASD_SC10, gamma = 50)
L4_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), L4_ASD_SC50$membership)
colnames(L4_ASD_SC50_GE) <- paste("MC", seq(ncol(L4_ASD_SC50_GE)), sep = "")

# L4_Normal
L4_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
L4_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), L4_Normal_SC10$membership)
colnames(L4_Normal_SC10_GE) <- paste("MC", seq(ncol(L4_Normal_SC10_GE)), sep = "")

L4_Normal_SC20 <- supercell_rescale(L4_Normal_SC10, gamma = 20)
L4_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), L4_Normal_SC20$membership)
colnames(L4_Normal_SC20_GE) <- paste("MC", seq(ncol(L4_Normal_SC20_GE)), sep = "")

L4_Normal_SC30 <- supercell_rescale(L4_Normal_SC10, gamma = 30)
L4_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), L4_Normal_SC30$membership)
colnames(L4_Normal_SC30_GE) <- paste("MC", seq(ncol(L4_Normal_SC30_GE)), sep = "")

L4_Normal_SC40 <- supercell_rescale(L4_Normal_SC10, gamma = 40)
L4_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), L4_Normal_SC40$membership)
colnames(L4_Normal_SC40_GE) <- paste("MC", seq(ncol(L4_Normal_SC40_GE)), sep = "")

L4_Normal_SC50 <- supercell_rescale(L4_Normal_SC10, gamma = 50)
L4_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), L4_Normal_SC50$membership)
colnames(L4_Normal_SC50_GE) <- paste("MC", seq(ncol(L4_Normal_SC50_GE)), sep = "")

# INSST_ASD
INSST_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INSST_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), INSST_ASD_SC10$membership)
colnames(INSST_ASD_SC10_GE) <- paste("MC", seq(ncol(INSST_ASD_SC10_GE)), sep = "")

INSST_ASD_SC20 <- supercell_rescale(INSST_ASD_SC10, gamma = 20)
INSST_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), INSST_ASD_SC20$membership)
colnames(INSST_ASD_SC20_GE) <- paste("MC", seq(ncol(INSST_ASD_SC20_GE)), sep = "")

INSST_ASD_SC30 <- supercell_rescale(INSST_ASD_SC10, gamma = 30)
INSST_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), INSST_ASD_SC30$membership)
colnames(INSST_ASD_SC30_GE) <- paste("MC", seq(ncol(INSST_ASD_SC30_GE)), sep = "")

INSST_ASD_SC40 <- supercell_rescale(INSST_ASD_SC10, gamma = 40)
INSST_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), INSST_ASD_SC40$membership)
colnames(INSST_ASD_SC40_GE) <- paste("MC", seq(ncol(INSST_ASD_SC40_GE)), sep = "")

INSST_ASD_SC50 <- supercell_rescale(INSST_ASD_SC10, gamma = 50)
INSST_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), INSST_ASD_SC50$membership)
colnames(INSST_ASD_SC50_GE) <- paste("MC", seq(ncol(INSST_ASD_SC50_GE)), sep = "")

# INSST_Normal
INSST_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
INSST_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), INSST_Normal_SC10$membership)
colnames(INSST_Normal_SC10_GE) <- paste("MC", seq(ncol(INSST_Normal_SC10_GE)), sep = "")

INSST_Normal_SC20 <- supercell_rescale(INSST_Normal_SC10, gamma = 20)
INSST_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), INSST_Normal_SC20$membership)
colnames(INSST_Normal_SC20_GE) <- paste("MC", seq(ncol(INSST_Normal_SC20_GE)), sep = "")

INSST_Normal_SC30 <- supercell_rescale(INSST_Normal_SC10, gamma = 30)
INSST_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), INSST_Normal_SC30$membership)
colnames(INSST_Normal_SC30_GE) <- paste("MC", seq(ncol(INSST_Normal_SC30_GE)), sep = "")

INSST_Normal_SC40 <- supercell_rescale(INSST_Normal_SC10, gamma = 40)
INSST_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), INSST_Normal_SC40$membership)
colnames(INSST_Normal_SC40_GE) <- paste("MC", seq(ncol(INSST_Normal_SC40_GE)), sep = "")

INSST_Normal_SC50 <- supercell_rescale(INSST_Normal_SC10, gamma = 50)
INSST_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), INSST_Normal_SC50$membership)
colnames(INSST_Normal_SC50_GE) <- paste("MC", seq(ncol(INSST_Normal_SC50_GE)), sep = "")

# Neumat_ASD
Neumat_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Neumat_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), Neumat_ASD_SC10$membership)
colnames(Neumat_ASD_SC10_GE) <- paste("MC", seq(ncol(Neumat_ASD_SC10_GE)), sep = "")

Neumat_ASD_SC20 <- supercell_rescale(Neumat_ASD_SC10, gamma = 20)
Neumat_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), Neumat_ASD_SC20$membership)
colnames(Neumat_ASD_SC20_GE) <- paste("MC", seq(ncol(Neumat_ASD_SC20_GE)), sep = "")

Neumat_ASD_SC30 <- supercell_rescale(Neumat_ASD_SC10, gamma = 30)
Neumat_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), Neumat_ASD_SC30$membership)
colnames(Neumat_ASD_SC30_GE) <- paste("MC", seq(ncol(Neumat_ASD_SC30_GE)), sep = "")

Neumat_ASD_SC40 <- supercell_rescale(Neumat_ASD_SC10, gamma = 40)
Neumat_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), Neumat_ASD_SC40$membership)
colnames(Neumat_ASD_SC40_GE) <- paste("MC", seq(ncol(Neumat_ASD_SC40_GE)), sep = "")

Neumat_ASD_SC50 <- supercell_rescale(Neumat_ASD_SC10, gamma = 50)
Neumat_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), Neumat_ASD_SC50$membership)
colnames(Neumat_ASD_SC50_GE) <- paste("MC", seq(ncol(Neumat_ASD_SC50_GE)), sep = "")

# Neumat_Normal
Neumat_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Neumat_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), Neumat_Normal_SC10$membership)
colnames(Neumat_Normal_SC10_GE) <- paste("MC", seq(ncol(Neumat_Normal_SC10_GE)), sep = "")

Neumat_Normal_SC20 <- supercell_rescale(Neumat_Normal_SC10, gamma = 20)
Neumat_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), Neumat_Normal_SC20$membership)
colnames(Neumat_Normal_SC20_GE) <- paste("MC", seq(ncol(Neumat_Normal_SC20_GE)), sep = "")

Neumat_Normal_SC30 <- supercell_rescale(Neumat_Normal_SC10, gamma = 30)
Neumat_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), Neumat_Normal_SC30$membership)
colnames(Neumat_Normal_SC30_GE) <- paste("MC", seq(ncol(Neumat_Normal_SC30_GE)), sep = "")

Neumat_Normal_SC40 <- supercell_rescale(Neumat_Normal_SC10, gamma = 40)
Neumat_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), Neumat_Normal_SC40$membership)
colnames(Neumat_Normal_SC40_GE) <- paste("MC", seq(ncol(Neumat_Normal_SC40_GE)), sep = "")

Neumat_Normal_SC50 <- supercell_rescale(Neumat_Normal_SC10, gamma = 50)
Neumat_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), Neumat_Normal_SC50$membership)
colnames(Neumat_Normal_SC50_GE) <- paste("MC", seq(ncol(Neumat_Normal_SC50_GE)), sep = "")

# ASTPP_ASD
ASTPP_ASD_SC10 <- SCimplify(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ASTPP_ASD_SC10_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), ASTPP_ASD_SC10$membership)
colnames(ASTPP_ASD_SC10_GE) <- paste("MC", seq(ncol(ASTPP_ASD_SC10_GE)), sep = "")

ASTPP_ASD_SC20 <- supercell_rescale(ASTPP_ASD_SC10, gamma = 20)
ASTPP_ASD_SC20_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), ASTPP_ASD_SC20$membership)
colnames(ASTPP_ASD_SC20_GE) <- paste("MC", seq(ncol(ASTPP_ASD_SC20_GE)), sep = "")

ASTPP_ASD_SC30 <- supercell_rescale(ASTPP_ASD_SC10, gamma = 30)
ASTPP_ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), ASTPP_ASD_SC30$membership)
colnames(ASTPP_ASD_SC30_GE) <- paste("MC", seq(ncol(ASTPP_ASD_SC30_GE)), sep = "")

ASTPP_ASD_SC40 <- supercell_rescale(ASTPP_ASD_SC10, gamma = 40)
ASTPP_ASD_SC40_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), ASTPP_ASD_SC40$membership)
colnames(ASTPP_ASD_SC40_GE) <- paste("MC", seq(ncol(ASTPP_ASD_SC40_GE)), sep = "")

ASTPP_ASD_SC50 <- supercell_rescale(ASTPP_ASD_SC10, gamma = 50)
ASTPP_ASD_SC50_GE <- supercell_GE(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), ASTPP_ASD_SC50$membership)
colnames(ASTPP_ASD_SC50_GE) <- paste("MC", seq(ncol(ASTPP_ASD_SC50_GE)), sep = "")

# ASTPP_Normal
ASTPP_Normal_SC10 <- SCimplify(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 10, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ASTPP_Normal_SC10_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), ASTPP_Normal_SC10$membership)
colnames(ASTPP_Normal_SC10_GE) <- paste("MC", seq(ncol(ASTPP_Normal_SC10_GE)), sep = "")

ASTPP_Normal_SC20 <- supercell_rescale(ASTPP_Normal_SC10, gamma = 20)
ASTPP_Normal_SC20_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), ASTPP_Normal_SC20$membership)
colnames(ASTPP_Normal_SC20_GE) <- paste("MC", seq(ncol(ASTPP_Normal_SC20_GE)), sep = "")

ASTPP_Normal_SC30 <- supercell_rescale(ASTPP_Normal_SC10, gamma = 30)
ASTPP_Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), ASTPP_Normal_SC30$membership)
colnames(ASTPP_Normal_SC30_GE) <- paste("MC", seq(ncol(ASTPP_Normal_SC30_GE)), sep = "")

ASTPP_Normal_SC40 <- supercell_rescale(ASTPP_Normal_SC10, gamma = 40)
ASTPP_Normal_SC40_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), ASTPP_Normal_SC40$membership)
colnames(ASTPP_Normal_SC40_GE) <- paste("MC", seq(ncol(ASTPP_Normal_SC40_GE)), sep = "")

ASTPP_Normal_SC50 <- supercell_rescale(ASTPP_Normal_SC10, gamma = 50)
ASTPP_Normal_SC50_GE <- supercell_GE(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), ASTPP_Normal_SC50$membership)
colnames(ASTPP_Normal_SC50_GE) <- paste("MC", seq(ncol(ASTPP_Normal_SC50_GE)), sep = "")

## Performing quality controls on metacells using MetacellAnalysisToolkit R package (https://github.com/GfellerLab/MetacellAnalysisToolkit)
# Four metrics are used, including Purity, Compactness, Separation, INV (inner normalized variance). 
# Purity is the proportion of the most abundant cell type within the metacell. The larger the purity value the better.
# ASD
ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)))), 2]
ASD_SC10_purity <- mc_purity(membership = ASD_SC10$membership, annotation = ASD_annotation)
ASD_SC20_purity <- mc_purity(membership = ASD_SC20$membership, annotation = ASD_annotation)
ASD_SC30_purity <- mc_purity(membership = ASD_SC30$membership, annotation = ASD_annotation)
ASD_SC40_purity <- mc_purity(membership = ASD_SC40$membership, annotation = ASD_annotation)
ASD_SC50_purity <- mc_purity(membership = ASD_SC50$membership, annotation = ASD_annotation)

# Normal
Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)))), 2]
Normal_SC10_purity <- mc_purity(membership = Normal_SC10$membership, annotation = Normal_annotation)
Normal_SC20_purity <- mc_purity(membership = Normal_SC20$membership, annotation = Normal_annotation)
Normal_SC30_purity <- mc_purity(membership = Normal_SC30$membership, annotation = Normal_annotation)
Normal_SC40_purity <- mc_purity(membership = Normal_SC40$membership, annotation = Normal_annotation)
Normal_SC50_purity <- mc_purity(membership = Normal_SC50$membership, annotation = Normal_annotation)

# ACC_ASD
ACC_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)))), 2]
ACC_ASD_SC10_purity <- mc_purity(membership = ACC_ASD_SC10$membership, annotation = ACC_ASD_annotation)
ACC_ASD_SC20_purity <- mc_purity(membership = ACC_ASD_SC20$membership, annotation = ACC_ASD_annotation)
ACC_ASD_SC30_purity <- mc_purity(membership = ACC_ASD_SC30$membership, annotation = ACC_ASD_annotation)
ACC_ASD_SC40_purity <- mc_purity(membership = ACC_ASD_SC40$membership, annotation = ACC_ASD_annotation)
ACC_ASD_SC50_purity <- mc_purity(membership = ACC_ASD_SC50$membership, annotation = ACC_ASD_annotation)

# ACC_Normal
ACC_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)))), 2]
ACC_Normal_SC10_purity <- mc_purity(membership = ACC_Normal_SC10$membership, annotation = ACC_Normal_annotation)
ACC_Normal_SC20_purity <- mc_purity(membership = ACC_Normal_SC20$membership, annotation = ACC_Normal_annotation)
ACC_Normal_SC30_purity <- mc_purity(membership = ACC_Normal_SC30$membership, annotation = ACC_Normal_annotation)
ACC_Normal_SC40_purity <- mc_purity(membership = ACC_Normal_SC40$membership, annotation = ACC_Normal_annotation)
ACC_Normal_SC50_purity <- mc_purity(membership = ACC_Normal_SC50$membership, annotation = ACC_Normal_annotation)

# PFC_ASD
PFC_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)))), 2]
PFC_ASD_SC10_purity <- mc_purity(membership = PFC_ASD_SC10$membership, annotation = PFC_ASD_annotation)
PFC_ASD_SC20_purity <- mc_purity(membership = PFC_ASD_SC20$membership, annotation = PFC_ASD_annotation)
PFC_ASD_SC30_purity <- mc_purity(membership = PFC_ASD_SC30$membership, annotation = PFC_ASD_annotation)
PFC_ASD_SC40_purity <- mc_purity(membership = PFC_ASD_SC40$membership, annotation = PFC_ASD_annotation)
PFC_ASD_SC50_purity <- mc_purity(membership = PFC_ASD_SC50$membership, annotation = PFC_ASD_annotation)

# PFC_Normal
PFC_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)))), 2]
PFC_Normal_SC10_purity <- mc_purity(membership = PFC_Normal_SC10$membership, annotation = PFC_Normal_annotation)
PFC_Normal_SC20_purity <- mc_purity(membership = PFC_Normal_SC20$membership, annotation = PFC_Normal_annotation)
PFC_Normal_SC30_purity <- mc_purity(membership = PFC_Normal_SC30$membership, annotation = PFC_Normal_annotation)
PFC_Normal_SC40_purity <- mc_purity(membership = PFC_Normal_SC40$membership, annotation = PFC_Normal_annotation)
PFC_Normal_SC50_purity <- mc_purity(membership = PFC_Normal_SC50$membership, annotation = PFC_Normal_annotation)

# Lower18_ASD
Lower18_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)))), 2]
Lower18_ASD_SC10_purity <- mc_purity(membership = Lower18_ASD_SC10$membership, annotation = Lower18_ASD_annotation)
Lower18_ASD_SC20_purity <- mc_purity(membership = Lower18_ASD_SC20$membership, annotation = Lower18_ASD_annotation)
Lower18_ASD_SC30_purity <- mc_purity(membership = Lower18_ASD_SC30$membership, annotation = Lower18_ASD_annotation)
Lower18_ASD_SC40_purity <- mc_purity(membership = Lower18_ASD_SC40$membership, annotation = Lower18_ASD_annotation)
Lower18_ASD_SC50_purity <- mc_purity(membership = Lower18_ASD_SC50$membership, annotation = Lower18_ASD_annotation)

# Lower18_Normal
Lower18_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)))), 2]
Lower18_Normal_SC10_purity <- mc_purity(membership = Lower18_Normal_SC10$membership, annotation = Lower18_Normal_annotation)
Lower18_Normal_SC20_purity <- mc_purity(membership = Lower18_Normal_SC20$membership, annotation = Lower18_Normal_annotation)
Lower18_Normal_SC30_purity <- mc_purity(membership = Lower18_Normal_SC30$membership, annotation = Lower18_Normal_annotation)
Lower18_Normal_SC40_purity <- mc_purity(membership = Lower18_Normal_SC40$membership, annotation = Lower18_Normal_annotation)
Lower18_Normal_SC50_purity <- mc_purity(membership = Lower18_Normal_SC50$membership, annotation = Lower18_Normal_annotation)

# Larger18_ASD
Larger18_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)))), 2]
Larger18_ASD_SC10_purity <- mc_purity(membership = Larger18_ASD_SC10$membership, annotation = Larger18_ASD_annotation)
Larger18_ASD_SC20_purity <- mc_purity(membership = Larger18_ASD_SC20$membership, annotation = Larger18_ASD_annotation)
Larger18_ASD_SC30_purity <- mc_purity(membership = Larger18_ASD_SC30$membership, annotation = Larger18_ASD_annotation)
Larger18_ASD_SC40_purity <- mc_purity(membership = Larger18_ASD_SC40$membership, annotation = Larger18_ASD_annotation)
Larger18_ASD_SC50_purity <- mc_purity(membership = Larger18_ASD_SC50$membership, annotation = Larger18_ASD_annotation)

# Larger18_Normal
Larger18_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)))), 2]
Larger18_Normal_SC10_purity <- mc_purity(membership = Larger18_Normal_SC10$membership, annotation = Larger18_Normal_annotation)
Larger18_Normal_SC20_purity <- mc_purity(membership = Larger18_Normal_SC20$membership, annotation = Larger18_Normal_annotation)
Larger18_Normal_SC30_purity <- mc_purity(membership = Larger18_Normal_SC30$membership, annotation = Larger18_Normal_annotation)
Larger18_Normal_SC40_purity <- mc_purity(membership = Larger18_Normal_SC40$membership, annotation = Larger18_Normal_annotation)
Larger18_Normal_SC50_purity <- mc_purity(membership = Larger18_Normal_SC50$membership, annotation = Larger18_Normal_annotation)

# Male_ASD
Male_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)))), 2]
Male_ASD_SC10_purity <- mc_purity(membership = Male_ASD_SC10$membership, annotation = Male_ASD_annotation)
Male_ASD_SC20_purity <- mc_purity(membership = Male_ASD_SC20$membership, annotation = Male_ASD_annotation)
Male_ASD_SC30_purity <- mc_purity(membership = Male_ASD_SC30$membership, annotation = Male_ASD_annotation)
Male_ASD_SC40_purity <- mc_purity(membership = Male_ASD_SC40$membership, annotation = Male_ASD_annotation)
Male_ASD_SC50_purity <- mc_purity(membership = Male_ASD_SC50$membership, annotation = Male_ASD_annotation)

# Male_Normal
Male_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)))), 2]
Male_Normal_SC10_purity <- mc_purity(membership = Male_Normal_SC10$membership, annotation = Male_Normal_annotation)
Male_Normal_SC20_purity <- mc_purity(membership = Male_Normal_SC20$membership, annotation = Male_Normal_annotation)
Male_Normal_SC30_purity <- mc_purity(membership = Male_Normal_SC30$membership, annotation = Male_Normal_annotation)
Male_Normal_SC40_purity <- mc_purity(membership = Male_Normal_SC40$membership, annotation = Male_Normal_annotation)
Male_Normal_SC50_purity <- mc_purity(membership = Male_Normal_SC50$membership, annotation = Male_Normal_annotation)

# Female_ASD
Female_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data))) ), 2]
Female_ASD_SC10_purity <- mc_purity(membership = Female_ASD_SC10$membership, annotation = Female_ASD_annotation)
Female_ASD_SC20_purity <- mc_purity(membership = Female_ASD_SC20$membership, annotation = Female_ASD_annotation)
Female_ASD_SC30_purity <- mc_purity(membership = Female_ASD_SC30$membership, annotation = Female_ASD_annotation)
Female_ASD_SC40_purity <- mc_purity(membership = Female_ASD_SC40$membership, annotation = Female_ASD_annotation)
Female_ASD_SC50_purity <- mc_purity(membership = Female_ASD_SC50$membership, annotation = Female_ASD_annotation)

# Female_Normal
Female_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data))) ), 2]
Female_Normal_SC10_purity <- mc_purity(membership = Female_Normal_SC10$membership, annotation = Female_Normal_annotation)
Female_Normal_SC20_purity <- mc_purity(membership = Female_Normal_SC20$membership, annotation = Female_Normal_annotation)
Female_Normal_SC30_purity <- mc_purity(membership = Female_Normal_SC30$membership, annotation = Female_Normal_annotation)
Female_Normal_SC40_purity <- mc_purity(membership = Female_Normal_SC40$membership, annotation = Female_Normal_annotation)
Female_Normal_SC50_purity <- mc_purity(membership = Female_Normal_SC50$membership, annotation = Female_Normal_annotation)

# NeuNRGNII_ASD
NeuNRGNII_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])))), 2]
NeuNRGNII_ASD_SC10_purity <- mc_purity(membership = NeuNRGNII_ASD_SC10$membership, annotation = NeuNRGNII_ASD_annotation)
NeuNRGNII_ASD_SC20_purity <- mc_purity(membership = NeuNRGNII_ASD_SC20$membership, annotation = NeuNRGNII_ASD_annotation)
NeuNRGNII_ASD_SC30_purity <- mc_purity(membership = NeuNRGNII_ASD_SC30$membership, annotation = NeuNRGNII_ASD_annotation)
NeuNRGNII_ASD_SC40_purity <- mc_purity(membership = NeuNRGNII_ASD_SC40$membership, annotation = NeuNRGNII_ASD_annotation)
NeuNRGNII_ASD_SC50_purity <- mc_purity(membership = NeuNRGNII_ASD_SC50$membership, annotation = NeuNRGNII_ASD_annotation)

# NeuNRGNII_Normal
NeuNRGNII_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])))), 2]
NeuNRGNII_Normal_SC10_purity <- mc_purity(membership = NeuNRGNII_Normal_SC10$membership, annotation = NeuNRGNII_Normal_annotation)
NeuNRGNII_Normal_SC20_purity <- mc_purity(membership = NeuNRGNII_Normal_SC20$membership, annotation = NeuNRGNII_Normal_annotation)
NeuNRGNII_Normal_SC30_purity <- mc_purity(membership = NeuNRGNII_Normal_SC30$membership, annotation = NeuNRGNII_Normal_annotation)
NeuNRGNII_Normal_SC40_purity <- mc_purity(membership = NeuNRGNII_Normal_SC40$membership, annotation = NeuNRGNII_Normal_annotation)
NeuNRGNII_Normal_SC50_purity <- mc_purity(membership = NeuNRGNII_Normal_SC50$membership, annotation = NeuNRGNII_Normal_annotation)

# L56_ASD
L56_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])))), 2]
L56_ASD_SC10_purity <- mc_purity(membership = L56_ASD_SC10$membership, annotation = L56_ASD_annotation)
L56_ASD_SC20_purity <- mc_purity(membership = L56_ASD_SC20$membership, annotation = L56_ASD_annotation)
L56_ASD_SC30_purity <- mc_purity(membership = L56_ASD_SC30$membership, annotation = L56_ASD_annotation)
L56_ASD_SC40_purity <- mc_purity(membership = L56_ASD_SC40$membership, annotation = L56_ASD_annotation)
L56_ASD_SC50_purity <- mc_purity(membership = L56_ASD_SC50$membership, annotation = L56_ASD_annotation)

# L56_Normal
L56_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])))), 2]
L56_Normal_SC10_purity <- mc_purity(membership = L56_Normal_SC10$membership, annotation = L56_Normal_annotation)
L56_Normal_SC20_purity <- mc_purity(membership = L56_Normal_SC20$membership, annotation = L56_Normal_annotation)
L56_Normal_SC30_purity <- mc_purity(membership = L56_Normal_SC30$membership, annotation = L56_Normal_annotation)
L56_Normal_SC40_purity <- mc_purity(membership = L56_Normal_SC40$membership, annotation = L56_Normal_annotation)
L56_Normal_SC50_purity <- mc_purity(membership = L56_Normal_SC50$membership, annotation = L56_Normal_annotation)

# Oligodendrocytes_ASD
Oligodendrocytes_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])))), 2]
Oligodendrocytes_ASD_SC10_purity <- mc_purity(membership = Oligodendrocytes_ASD_SC10$membership, annotation = Oligodendrocytes_ASD_annotation)
Oligodendrocytes_ASD_SC20_purity <- mc_purity(membership = Oligodendrocytes_ASD_SC20$membership, annotation = Oligodendrocytes_ASD_annotation)
Oligodendrocytes_ASD_SC30_purity <- mc_purity(membership = Oligodendrocytes_ASD_SC30$membership, annotation = Oligodendrocytes_ASD_annotation)
Oligodendrocytes_ASD_SC40_purity <- mc_purity(membership = Oligodendrocytes_ASD_SC40$membership, annotation = Oligodendrocytes_ASD_annotation)
Oligodendrocytes_ASD_SC50_purity <- mc_purity(membership = Oligodendrocytes_ASD_SC50$membership, annotation = Oligodendrocytes_ASD_annotation)

# Oligodendrocytes_Normal
Oligodendrocytes_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])))), 2]
Oligodendrocytes_Normal_SC10_purity <- mc_purity(membership = Oligodendrocytes_Normal_SC10$membership, annotation = Oligodendrocytes_Normal_annotation)
Oligodendrocytes_Normal_SC20_purity <- mc_purity(membership = Oligodendrocytes_Normal_SC20$membership, annotation = Oligodendrocytes_Normal_annotation)
Oligodendrocytes_Normal_SC30_purity <- mc_purity(membership = Oligodendrocytes_Normal_SC30$membership, annotation = Oligodendrocytes_Normal_annotation)
Oligodendrocytes_Normal_SC40_purity <- mc_purity(membership = Oligodendrocytes_Normal_SC40$membership, annotation = Oligodendrocytes_Normal_annotation)
Oligodendrocytes_Normal_SC50_purity <- mc_purity(membership = Oligodendrocytes_Normal_SC50$membership, annotation = Oligodendrocytes_Normal_annotation)

# OPC_ASD
OPC_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])))), 2]
OPC_ASD_SC10_purity <- mc_purity(membership = OPC_ASD_SC10$membership, annotation = OPC_ASD_annotation)
OPC_ASD_SC20_purity <- mc_purity(membership = OPC_ASD_SC20$membership, annotation = OPC_ASD_annotation)
OPC_ASD_SC30_purity <- mc_purity(membership = OPC_ASD_SC30$membership, annotation = OPC_ASD_annotation)
OPC_ASD_SC40_purity <- mc_purity(membership = OPC_ASD_SC40$membership, annotation = OPC_ASD_annotation)
OPC_ASD_SC50_purity <- mc_purity(membership = OPC_ASD_SC50$membership, annotation = OPC_ASD_annotation)

# OPC_Normal
OPC_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])))), 2]
OPC_Normal_SC10_purity <- mc_purity(membership = OPC_Normal_SC10$membership, annotation = OPC_Normal_annotation)
OPC_Normal_SC20_purity <- mc_purity(membership = OPC_Normal_SC20$membership, annotation = OPC_Normal_annotation)
OPC_Normal_SC30_purity <- mc_purity(membership = OPC_Normal_SC30$membership, annotation = OPC_Normal_annotation)
OPC_Normal_SC40_purity <- mc_purity(membership = OPC_Normal_SC40$membership, annotation = OPC_Normal_annotation)
OPC_Normal_SC50_purity <- mc_purity(membership = OPC_Normal_SC50$membership, annotation = OPC_Normal_annotation)

# ASTFB_ASD
ASTFB_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])))), 2]
ASTFB_ASD_SC10_purity <- mc_purity(membership = ASTFB_ASD_SC10$membership, annotation = ASTFB_ASD_annotation)
ASTFB_ASD_SC20_purity <- mc_purity(membership = ASTFB_ASD_SC20$membership, annotation = ASTFB_ASD_annotation)
ASTFB_ASD_SC30_purity <- mc_purity(membership = ASTFB_ASD_SC30$membership, annotation = ASTFB_ASD_annotation)
ASTFB_ASD_SC40_purity <- mc_purity(membership = ASTFB_ASD_SC40$membership, annotation = ASTFB_ASD_annotation)
ASTFB_ASD_SC50_purity <- mc_purity(membership = ASTFB_ASD_SC50$membership, annotation = ASTFB_ASD_annotation)

# ASTFB_Normal
ASTFB_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])))), 2]
ASTFB_Normal_SC10_purity <- mc_purity(membership = ASTFB_Normal_SC10$membership, annotation = ASTFB_Normal_annotation)
ASTFB_Normal_SC20_purity <- mc_purity(membership = ASTFB_Normal_SC20$membership, annotation = ASTFB_Normal_annotation)
ASTFB_Normal_SC30_purity <- mc_purity(membership = ASTFB_Normal_SC30$membership, annotation = ASTFB_Normal_annotation)
ASTFB_Normal_SC40_purity <- mc_purity(membership = ASTFB_Normal_SC40$membership, annotation = ASTFB_Normal_annotation)
ASTFB_Normal_SC50_purity <- mc_purity(membership = ASTFB_Normal_SC50$membership, annotation = ASTFB_Normal_annotation)

# Endothelial_ASD
Endothelial_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])))), 2]
Endothelial_ASD_SC10_purity <- mc_purity(membership = Endothelial_ASD_SC10$membership, annotation = Endothelial_ASD_annotation)
Endothelial_ASD_SC20_purity <- mc_purity(membership = Endothelial_ASD_SC20$membership, annotation = Endothelial_ASD_annotation)
Endothelial_ASD_SC30_purity <- mc_purity(membership = Endothelial_ASD_SC30$membership, annotation = Endothelial_ASD_annotation)
Endothelial_ASD_SC40_purity <- mc_purity(membership = Endothelial_ASD_SC40$membership, annotation = Endothelial_ASD_annotation)
Endothelial_ASD_SC50_purity <- mc_purity(membership = Endothelial_ASD_SC50$membership, annotation = Endothelial_ASD_annotation)

# Endothelial_Normal
Endothelial_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])))), 2]
Endothelial_Normal_SC10_purity <- mc_purity(membership = Endothelial_Normal_SC10$membership, annotation = Endothelial_Normal_annotation)
Endothelial_Normal_SC20_purity <- mc_purity(membership = Endothelial_Normal_SC20$membership, annotation = Endothelial_Normal_annotation)
Endothelial_Normal_SC30_purity <- mc_purity(membership = Endothelial_Normal_SC30$membership, annotation = Endothelial_Normal_annotation)
Endothelial_Normal_SC40_purity <- mc_purity(membership = Endothelial_Normal_SC40$membership, annotation = Endothelial_Normal_annotation)
Endothelial_Normal_SC50_purity <- mc_purity(membership = Endothelial_Normal_SC50$membership, annotation = Endothelial_Normal_annotation)

# Microglia_ASD
Microglia_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])))), 2]
Microglia_ASD_SC10_purity <- mc_purity(membership = Microglia_ASD_SC10$membership, annotation = Microglia_ASD_annotation)
Microglia_ASD_SC20_purity <- mc_purity(membership = Microglia_ASD_SC20$membership, annotation = Microglia_ASD_annotation)
Microglia_ASD_SC30_purity <- mc_purity(membership = Microglia_ASD_SC30$membership, annotation = Microglia_ASD_annotation)
Microglia_ASD_SC40_purity <- mc_purity(membership = Microglia_ASD_SC40$membership, annotation = Microglia_ASD_annotation)
Microglia_ASD_SC50_purity <- mc_purity(membership = Microglia_ASD_SC50$membership, annotation = Microglia_ASD_annotation)

# Microglia_Normal
Microglia_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])))), 2]
Microglia_Normal_SC10_purity <- mc_purity(membership = Microglia_Normal_SC10$membership, annotation = Microglia_Normal_annotation)
Microglia_Normal_SC20_purity <- mc_purity(membership = Microglia_Normal_SC20$membership, annotation = Microglia_Normal_annotation)
Microglia_Normal_SC30_purity <- mc_purity(membership = Microglia_Normal_SC30$membership, annotation = Microglia_Normal_annotation)
Microglia_Normal_SC40_purity <- mc_purity(membership = Microglia_Normal_SC40$membership, annotation = Microglia_Normal_annotation)
Microglia_Normal_SC50_purity <- mc_purity(membership = Microglia_Normal_SC50$membership, annotation = Microglia_Normal_annotation)

# NeuNRGNI_ASD
NeuNRGNI_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])))), 2]
NeuNRGNI_ASD_SC10_purity <- mc_purity(membership = NeuNRGNI_ASD_SC10$membership, annotation = NeuNRGNI_ASD_annotation)
NeuNRGNI_ASD_SC20_purity <- mc_purity(membership = NeuNRGNI_ASD_SC20$membership, annotation = NeuNRGNI_ASD_annotation)
NeuNRGNI_ASD_SC30_purity <- mc_purity(membership = NeuNRGNI_ASD_SC30$membership, annotation = NeuNRGNI_ASD_annotation)
NeuNRGNI_ASD_SC40_purity <- mc_purity(membership = NeuNRGNI_ASD_SC40$membership, annotation = NeuNRGNI_ASD_annotation)
NeuNRGNI_ASD_SC50_purity <- mc_purity(membership = NeuNRGNI_ASD_SC50$membership, annotation = NeuNRGNI_ASD_annotation)

# NeuNRGNI_Normal
NeuNRGNI_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])))), 2]
NeuNRGNI_Normal_SC10_purity <- mc_purity(membership = NeuNRGNI_Normal_SC10$membership, annotation = NeuNRGNI_Normal_annotation)
NeuNRGNI_Normal_SC20_purity <- mc_purity(membership = NeuNRGNI_Normal_SC20$membership, annotation = NeuNRGNI_Normal_annotation)
NeuNRGNI_Normal_SC30_purity <- mc_purity(membership = NeuNRGNI_Normal_SC30$membership, annotation = NeuNRGNI_Normal_annotation)
NeuNRGNI_Normal_SC40_purity <- mc_purity(membership = NeuNRGNI_Normal_SC40$membership, annotation = NeuNRGNI_Normal_annotation)
NeuNRGNI_Normal_SC50_purity <- mc_purity(membership = NeuNRGNI_Normal_SC50$membership, annotation = NeuNRGNI_Normal_annotation)

# INVIP_ASD
INVIP_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])))), 2]
INVIP_ASD_SC10_purity <- mc_purity(membership = INVIP_ASD_SC10$membership, annotation = INVIP_ASD_annotation)
INVIP_ASD_SC20_purity <- mc_purity(membership = INVIP_ASD_SC20$membership, annotation = INVIP_ASD_annotation)
INVIP_ASD_SC30_purity <- mc_purity(membership = INVIP_ASD_SC30$membership, annotation = INVIP_ASD_annotation)
INVIP_ASD_SC40_purity <- mc_purity(membership = INVIP_ASD_SC40$membership, annotation = INVIP_ASD_annotation)
INVIP_ASD_SC50_purity <- mc_purity(membership = INVIP_ASD_SC50$membership, annotation = INVIP_ASD_annotation)

# INVIP_Normal
INVIP_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])))), 2]
INVIP_Normal_SC10_purity <- mc_purity(membership = INVIP_Normal_SC10$membership, annotation = INVIP_Normal_annotation)
INVIP_Normal_SC20_purity <- mc_purity(membership = INVIP_Normal_SC20$membership, annotation = INVIP_Normal_annotation)
INVIP_Normal_SC30_purity <- mc_purity(membership = INVIP_Normal_SC30$membership, annotation = INVIP_Normal_annotation)
INVIP_Normal_SC40_purity <- mc_purity(membership = INVIP_Normal_SC40$membership, annotation = INVIP_Normal_annotation)
INVIP_Normal_SC50_purity <- mc_purity(membership = INVIP_Normal_SC50$membership, annotation = INVIP_Normal_annotation)

# L56CC_ASD
L56CC_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])))), 2]
L56CC_ASD_SC10_purity <- mc_purity(membership = L56CC_ASD_SC10$membership, annotation = L56CC_ASD_annotation)
L56CC_ASD_SC20_purity <- mc_purity(membership = L56CC_ASD_SC20$membership, annotation = L56CC_ASD_annotation)
L56CC_ASD_SC30_purity <- mc_purity(membership = L56CC_ASD_SC30$membership, annotation = L56CC_ASD_annotation)
L56CC_ASD_SC40_purity <- mc_purity(membership = L56CC_ASD_SC40$membership, annotation = L56CC_ASD_annotation)
L56CC_ASD_SC50_purity <- mc_purity(membership = L56CC_ASD_SC50$membership, annotation = L56CC_ASD_annotation)

# L56CC_Normal
L56CC_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])))), 2]
L56CC_Normal_SC10_purity <- mc_purity(membership = L56CC_Normal_SC10$membership, annotation = L56CC_Normal_annotation)
L56CC_Normal_SC20_purity <- mc_purity(membership = L56CC_Normal_SC20$membership, annotation = L56CC_Normal_annotation)
L56CC_Normal_SC30_purity <- mc_purity(membership = L56CC_Normal_SC30$membership, annotation = L56CC_Normal_annotation)
L56CC_Normal_SC40_purity <- mc_purity(membership = L56CC_Normal_SC40$membership, annotation = L56CC_Normal_annotation)
L56CC_Normal_SC50_purity <- mc_purity(membership = L56CC_Normal_SC50$membership, annotation = L56CC_Normal_annotation)

# INSV2C_ASD
INSV2C_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])))), 2]
INSV2C_ASD_SC10_purity <- mc_purity(membership = INSV2C_ASD_SC10$membership, annotation = INSV2C_ASD_annotation)
INSV2C_ASD_SC20_purity <- mc_purity(membership = INSV2C_ASD_SC20$membership, annotation = INSV2C_ASD_annotation)
INSV2C_ASD_SC30_purity <- mc_purity(membership = INSV2C_ASD_SC30$membership, annotation = INSV2C_ASD_annotation)
INSV2C_ASD_SC40_purity <- mc_purity(membership = INSV2C_ASD_SC40$membership, annotation = INSV2C_ASD_annotation)
INSV2C_ASD_SC50_purity <- mc_purity(membership = INSV2C_ASD_SC50$membership, annotation = INSV2C_ASD_annotation)

# INSV2C_Normal
INSV2C_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])))), 2]
INSV2C_Normal_SC10_purity <- mc_purity(membership = INSV2C_Normal_SC10$membership, annotation = INSV2C_Normal_annotation)
INSV2C_Normal_SC20_purity <- mc_purity(membership = INSV2C_Normal_SC20$membership, annotation = INSV2C_Normal_annotation)
INSV2C_Normal_SC30_purity <- mc_purity(membership = INSV2C_Normal_SC30$membership, annotation = INSV2C_Normal_annotation)
INSV2C_Normal_SC40_purity <- mc_purity(membership = INSV2C_Normal_SC40$membership, annotation = INSV2C_Normal_annotation)
INSV2C_Normal_SC50_purity <- mc_purity(membership = INSV2C_Normal_SC50$membership, annotation = INSV2C_Normal_annotation)

# L23_ASD
L23_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])))), 2]
L23_ASD_SC10_purity <- mc_purity(membership = L23_ASD_SC10$membership, annotation = L23_ASD_annotation)
L23_ASD_SC20_purity <- mc_purity(membership = L23_ASD_SC20$membership, annotation = L23_ASD_annotation)
L23_ASD_SC30_purity <- mc_purity(membership = L23_ASD_SC30$membership, annotation = L23_ASD_annotation)
L23_ASD_SC40_purity <- mc_purity(membership = L23_ASD_SC40$membership, annotation = L23_ASD_annotation)
L23_ASD_SC50_purity <- mc_purity(membership = L23_ASD_SC50$membership, annotation = L23_ASD_annotation)

# L23_Normal
L23_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])))), 2]
L23_Normal_SC10_purity <- mc_purity(membership = L23_Normal_SC10$membership, annotation = L23_Normal_annotation)
L23_Normal_SC20_purity <- mc_purity(membership = L23_Normal_SC20$membership, annotation = L23_Normal_annotation)
L23_Normal_SC30_purity <- mc_purity(membership = L23_Normal_SC30$membership, annotation = L23_Normal_annotation)
L23_Normal_SC40_purity <- mc_purity(membership = L23_Normal_SC40$membership, annotation = L23_Normal_annotation)
L23_Normal_SC50_purity <- mc_purity(membership = L23_Normal_SC50$membership, annotation = L23_Normal_annotation)

# INPV_ASD
INPV_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])))), 2]
INPV_ASD_SC10_purity <- mc_purity(membership = INPV_ASD_SC10$membership, annotation = INPV_ASD_annotation)
INPV_ASD_SC20_purity <- mc_purity(membership = INPV_ASD_SC20$membership, annotation = INPV_ASD_annotation)
INPV_ASD_SC30_purity <- mc_purity(membership = INPV_ASD_SC30$membership, annotation = INPV_ASD_annotation)
INPV_ASD_SC40_purity <- mc_purity(membership = INPV_ASD_SC40$membership, annotation = INPV_ASD_annotation)
INPV_ASD_SC50_purity <- mc_purity(membership = INPV_ASD_SC50$membership, annotation = INPV_ASD_annotation)

# INPV_Normal
INPV_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])))), 2]
INPV_Normal_SC10_purity <- mc_purity(membership = INPV_Normal_SC10$membership, annotation = INPV_Normal_annotation)
INPV_Normal_SC20_purity <- mc_purity(membership = INPV_Normal_SC20$membership, annotation = INPV_Normal_annotation)
INPV_Normal_SC30_purity <- mc_purity(membership = INPV_Normal_SC30$membership, annotation = INPV_Normal_annotation)
INPV_Normal_SC40_purity <- mc_purity(membership = INPV_Normal_SC40$membership, annotation = INPV_Normal_annotation)
INPV_Normal_SC50_purity <- mc_purity(membership = INPV_Normal_SC50$membership, annotation = INPV_Normal_annotation)

# L4_ASD
L4_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])))), 2]
L4_ASD_SC10_purity <- mc_purity(membership = L4_ASD_SC10$membership, annotation = L4_ASD_annotation)
L4_ASD_SC20_purity <- mc_purity(membership = L4_ASD_SC20$membership, annotation = L4_ASD_annotation)
L4_ASD_SC30_purity <- mc_purity(membership = L4_ASD_SC30$membership, annotation = L4_ASD_annotation)
L4_ASD_SC40_purity <- mc_purity(membership = L4_ASD_SC40$membership, annotation = L4_ASD_annotation)
L4_ASD_SC50_purity <- mc_purity(membership = L4_ASD_SC50$membership, annotation = L4_ASD_annotation)

# L4_Normal
L4_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])))), 2]
L4_Normal_SC10_purity <- mc_purity(membership = L4_Normal_SC10$membership, annotation = L4_Normal_annotation)
L4_Normal_SC20_purity <- mc_purity(membership = L4_Normal_SC20$membership, annotation = L4_Normal_annotation)
L4_Normal_SC30_purity <- mc_purity(membership = L4_Normal_SC30$membership, annotation = L4_Normal_annotation)
L4_Normal_SC40_purity <- mc_purity(membership = L4_Normal_SC40$membership, annotation = L4_Normal_annotation)
L4_Normal_SC50_purity <- mc_purity(membership = L4_Normal_SC50$membership, annotation = L4_Normal_annotation)

# INSST_ASD
INSST_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])))), 2]
INSST_ASD_SC10_purity <- mc_purity(membership = INSST_ASD_SC10$membership, annotation = INSST_ASD_annotation)
INSST_ASD_SC20_purity <- mc_purity(membership = INSST_ASD_SC20$membership, annotation = INSST_ASD_annotation)
INSST_ASD_SC30_purity <- mc_purity(membership = INSST_ASD_SC30$membership, annotation = INSST_ASD_annotation)
INSST_ASD_SC40_purity <- mc_purity(membership = INSST_ASD_SC40$membership, annotation = INSST_ASD_annotation)
INSST_ASD_SC50_purity <- mc_purity(membership = INSST_ASD_SC50$membership, annotation = INSST_ASD_annotation)

# INSST_Normal
INSST_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])))), 2]
INSST_Normal_SC10_purity <- mc_purity(membership = INSST_Normal_SC10$membership, annotation = INSST_Normal_annotation)
INSST_Normal_SC20_purity <- mc_purity(membership = INSST_Normal_SC20$membership, annotation = INSST_Normal_annotation)
INSST_Normal_SC30_purity <- mc_purity(membership = INSST_Normal_SC30$membership, annotation = INSST_Normal_annotation)
INSST_Normal_SC40_purity <- mc_purity(membership = INSST_Normal_SC40$membership, annotation = INSST_Normal_annotation)
INSST_Normal_SC50_purity <- mc_purity(membership = INSST_Normal_SC50$membership, annotation = INSST_Normal_annotation)

# Neumat_ASD
Neumat_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])))), 2]
Neumat_ASD_SC10_purity <- mc_purity(membership = Neumat_ASD_SC10$membership, annotation = Neumat_ASD_annotation)
Neumat_ASD_SC20_purity <- mc_purity(membership = Neumat_ASD_SC20$membership, annotation = Neumat_ASD_annotation)
Neumat_ASD_SC30_purity <- mc_purity(membership = Neumat_ASD_SC30$membership, annotation = Neumat_ASD_annotation)
Neumat_ASD_SC40_purity <- mc_purity(membership = Neumat_ASD_SC40$membership, annotation = Neumat_ASD_annotation)
Neumat_ASD_SC50_purity <- mc_purity(membership = Neumat_ASD_SC50$membership, annotation = Neumat_ASD_annotation)

# Neumat_Normal
Neumat_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])))), 2]
Neumat_Normal_SC10_purity <- mc_purity(membership = Neumat_Normal_SC10$membership, annotation = Neumat_Normal_annotation)
Neumat_Normal_SC20_purity <- mc_purity(membership = Neumat_Normal_SC20$membership, annotation = Neumat_Normal_annotation)
Neumat_Normal_SC30_purity <- mc_purity(membership = Neumat_Normal_SC30$membership, annotation = Neumat_Normal_annotation)
Neumat_Normal_SC40_purity <- mc_purity(membership = Neumat_Normal_SC40$membership, annotation = Neumat_Normal_annotation)
Neumat_Normal_SC50_purity <- mc_purity(membership = Neumat_Normal_SC50$membership, annotation = Neumat_Normal_annotation)

# ASTPP_ASD
ASTPP_ASD_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])))), 2]
ASTPP_ASD_SC10_purity <- mc_purity(membership = ASTPP_ASD_SC10$membership, annotation = ASTPP_ASD_annotation)
ASTPP_ASD_SC20_purity <- mc_purity(membership = ASTPP_ASD_SC20$membership, annotation = ASTPP_ASD_annotation)
ASTPP_ASD_SC30_purity <- mc_purity(membership = ASTPP_ASD_SC30$membership, annotation = ASTPP_ASD_annotation)
ASTPP_ASD_SC40_purity <- mc_purity(membership = ASTPP_ASD_SC40$membership, annotation = ASTPP_ASD_annotation)
ASTPP_ASD_SC50_purity <- mc_purity(membership = ASTPP_ASD_SC50$membership, annotation = ASTPP_ASD_annotation)

# ASTPP_Normal
ASTPP_Normal_annotation <- meta.cell[which(meta.cell[, 1] %in% colnames(as.matrix(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])))), 2]
ASTPP_Normal_SC10_purity <- mc_purity(membership = ASTPP_Normal_SC10$membership, annotation = ASTPP_Normal_annotation)
ASTPP_Normal_SC20_purity <- mc_purity(membership = ASTPP_Normal_SC20$membership, annotation = ASTPP_Normal_annotation)
ASTPP_Normal_SC30_purity <- mc_purity(membership = ASTPP_Normal_SC30$membership, annotation = ASTPP_Normal_annotation)
ASTPP_Normal_SC40_purity <- mc_purity(membership = ASTPP_Normal_SC40$membership, annotation = ASTPP_Normal_annotation)
ASTPP_Normal_SC50_purity <- mc_purity(membership = ASTPP_Normal_SC50$membership, annotation = ASTPP_Normal_annotation)

# Compactness is the variance of the components within the metacell. The lower the compactness value the better.
# ASD
ASD_SC10_membership <- as.matrix(ASD_SC10$membership)
colnames(ASD_SC10_membership) <- c("membership")
ASD_SC20_membership <- as.matrix(ASD_SC20$membership)
colnames(ASD_SC20_membership) <- c("membership")
ASD_SC30_membership <- as.matrix(ASD_SC30$membership)
colnames(ASD_SC30_membership) <- c("membership")
ASD_SC40_membership <- as.matrix(ASD_SC40$membership)
colnames(ASD_SC40_membership) <- c("membership")
ASD_SC50_membership <- as.matrix(ASD_SC50$membership)
colnames(ASD_SC50_membership) <- c("membership")

# Normal
Normal_SC10_membership <- as.matrix(Normal_SC10$membership)
colnames(Normal_SC10_membership) <- c("membership")
Normal_SC20_membership <- as.matrix(Normal_SC20$membership)
colnames(Normal_SC20_membership) <- c("membership")
Normal_SC30_membership <- as.matrix(Normal_SC30$membership)
colnames(Normal_SC30_membership) <- c("membership")
Normal_SC40_membership <- as.matrix(Normal_SC40$membership)
colnames(Normal_SC40_membership) <- c("membership")
Normal_SC50_membership <- as.matrix(Normal_SC50$membership)
colnames(Normal_SC50_membership) <- c("membership")

# ACC_ASD
ACC_ASD_SC10_membership <- as.matrix(ACC_ASD_SC10$membership)
colnames(ACC_ASD_SC10_membership) <- c("membership")
ACC_ASD_SC20_membership <- as.matrix(ACC_ASD_SC20$membership)
colnames(ACC_ASD_SC20_membership) <- c("membership")
ACC_ASD_SC30_membership <- as.matrix(ACC_ASD_SC30$membership)
colnames(ACC_ASD_SC30_membership) <- c("membership")
ACC_ASD_SC40_membership <- as.matrix(ACC_ASD_SC40$membership)
colnames(ACC_ASD_SC40_membership) <- c("membership")
ACC_ASD_SC50_membership <- as.matrix(ACC_ASD_SC50$membership)
colnames(ACC_ASD_SC50_membership) <- c("membership")

# ACC_Normal
ACC_Normal_SC10_membership <- as.matrix(ACC_Normal_SC10$membership)
colnames(ACC_Normal_SC10_membership) <- c("membership")
ACC_Normal_SC20_membership <- as.matrix(ACC_Normal_SC20$membership)
colnames(ACC_Normal_SC20_membership) <- c("membership")
ACC_Normal_SC30_membership <- as.matrix(ACC_Normal_SC30$membership)
colnames(ACC_Normal_SC30_membership) <- c("membership")
ACC_Normal_SC40_membership <- as.matrix(ACC_Normal_SC40$membership)
colnames(ACC_Normal_SC40_membership) <- c("membership")
ACC_Normal_SC50_membership <- as.matrix(ACC_Normal_SC50$membership)
colnames(ACC_Normal_SC50_membership) <- c("membership")

# PFC_ASD
PFC_ASD_SC10_membership <- as.matrix(PFC_ASD_SC10$membership)
colnames(PFC_ASD_SC10_membership) <- c("membership")
PFC_ASD_SC20_membership <- as.matrix(PFC_ASD_SC20$membership)
colnames(PFC_ASD_SC20_membership) <- c("membership")
PFC_ASD_SC30_membership <- as.matrix(PFC_ASD_SC30$membership)
colnames(PFC_ASD_SC30_membership) <- c("membership")
PFC_ASD_SC40_membership <- as.matrix(PFC_ASD_SC40$membership)
colnames(PFC_ASD_SC40_membership) <- c("membership")
PFC_ASD_SC50_membership <- as.matrix(PFC_ASD_SC50$membership)
colnames(PFC_ASD_SC50_membership) <- c("membership")

# PFC_Normal
PFC_Normal_SC10_membership <- as.matrix(PFC_Normal_SC10$membership)
colnames(PFC_Normal_SC10_membership) <- c("membership")
PFC_Normal_SC20_membership <- as.matrix(PFC_Normal_SC20$membership)
colnames(PFC_Normal_SC20_membership) <- c("membership")
PFC_Normal_SC30_membership <- as.matrix(PFC_Normal_SC30$membership)
colnames(PFC_Normal_SC30_membership) <- c("membership")
PFC_Normal_SC40_membership <- as.matrix(PFC_Normal_SC40$membership)
colnames(PFC_Normal_SC40_membership) <- c("membership")
PFC_Normal_SC50_membership <- as.matrix(PFC_Normal_SC50$membership)
colnames(PFC_Normal_SC50_membership) <- c("membership")

# Lower18_ASD
Lower18_ASD_SC10_membership <- as.matrix(Lower18_ASD_SC10$membership)
colnames(Lower18_ASD_SC10_membership) <- c("membership")
Lower18_ASD_SC20_membership <- as.matrix(Lower18_ASD_SC20$membership)
colnames(Lower18_ASD_SC20_membership) <- c("membership")
Lower18_ASD_SC30_membership <- as.matrix(Lower18_ASD_SC30$membership)
colnames(Lower18_ASD_SC30_membership) <- c("membership")
Lower18_ASD_SC40_membership <- as.matrix(Lower18_ASD_SC40$membership)
colnames(Lower18_ASD_SC40_membership) <- c("membership")
Lower18_ASD_SC50_membership <- as.matrix(Lower18_ASD_SC50$membership)
colnames(Lower18_ASD_SC50_membership) <- c("membership")

# Lower18_Normal
Lower18_Normal_SC10_membership <- as.matrix(Lower18_Normal_SC10$membership)
colnames(Lower18_Normal_SC10_membership) <- c("membership")
Lower18_Normal_SC20_membership <- as.matrix(Lower18_Normal_SC20$membership)
colnames(Lower18_Normal_SC20_membership) <- c("membership")
Lower18_Normal_SC30_membership <- as.matrix(Lower18_Normal_SC30$membership)
colnames(Lower18_Normal_SC30_membership) <- c("membership")
Lower18_Normal_SC40_membership <- as.matrix(Lower18_Normal_SC40$membership)
colnames(Lower18_Normal_SC40_membership) <- c("membership")
Lower18_Normal_SC50_membership <- as.matrix(Lower18_Normal_SC50$membership)
colnames(Lower18_Normal_SC50_membership) <- c("membership")

# Larger18_ASD
Larger18_ASD_SC10_membership <- as.matrix(Larger18_ASD_SC10$membership)
colnames(Larger18_ASD_SC10_membership) <- c("membership")
Larger18_ASD_SC20_membership <- as.matrix(Larger18_ASD_SC20$membership)
colnames(Larger18_ASD_SC20_membership) <- c("membership")
Larger18_ASD_SC30_membership <- as.matrix(Larger18_ASD_SC30$membership)
colnames(Larger18_ASD_SC30_membership) <- c("membership")
Larger18_ASD_SC40_membership <- as.matrix(Larger18_ASD_SC40$membership)
colnames(Larger18_ASD_SC40_membership) <- c("membership")
Larger18_ASD_SC50_membership <- as.matrix(Larger18_ASD_SC50$membership)
colnames(Larger18_ASD_SC50_membership) <- c("membership")

# Larger18_Normal
Larger18_Normal_SC10_membership <- as.matrix(Larger18_Normal_SC10$membership)
colnames(Larger18_Normal_SC10_membership) <- c("membership")
Larger18_Normal_SC20_membership <- as.matrix(Larger18_Normal_SC20$membership)
colnames(Larger18_Normal_SC20_membership) <- c("membership")
Larger18_Normal_SC30_membership <- as.matrix(Larger18_Normal_SC30$membership)
colnames(Larger18_Normal_SC30_membership) <- c("membership")
Larger18_Normal_SC40_membership <- as.matrix(Larger18_Normal_SC40$membership)
colnames(Larger18_Normal_SC40_membership) <- c("membership")
Larger18_Normal_SC50_membership <- as.matrix(Larger18_Normal_SC50$membership)
colnames(Larger18_Normal_SC50_membership) <- c("membership")

# Male_ASD
Male_ASD_SC10_membership <- as.matrix(Male_ASD_SC10$membership)
colnames(Male_ASD_SC10_membership) <- c("membership")
Male_ASD_SC20_membership <- as.matrix(Male_ASD_SC20$membership)
colnames(Male_ASD_SC20_membership) <- c("membership")
Male_ASD_SC30_membership <- as.matrix(Male_ASD_SC30$membership)
colnames(Male_ASD_SC30_membership) <- c("membership")
Male_ASD_SC40_membership <- as.matrix(Male_ASD_SC40$membership)
colnames(Male_ASD_SC40_membership) <- c("membership")
Male_ASD_SC50_membership <- as.matrix(Male_ASD_SC50$membership)
colnames(Male_ASD_SC50_membership) <- c("membership")

# Male_Normal
Male_Normal_SC10_membership <- as.matrix(Male_Normal_SC10$membership)
colnames(Male_Normal_SC10_membership) <- c("membership")
Male_Normal_SC20_membership <- as.matrix(Male_Normal_SC20$membership)
colnames(Male_Normal_SC20_membership) <- c("membership")
Male_Normal_SC30_membership <- as.matrix(Male_Normal_SC30$membership)
colnames(Male_Normal_SC30_membership) <- c("membership")
Male_Normal_SC40_membership <- as.matrix(Male_Normal_SC40$membership)
colnames(Male_Normal_SC40_membership) <- c("membership")
Male_Normal_SC50_membership <- as.matrix(Male_Normal_SC50$membership)
colnames(Male_Normal_SC50_membership) <- c("membership")

# Female_ASD
Female_ASD_SC10_membership <- as.matrix(Female_ASD_SC10$membership)
colnames(Female_ASD_SC10_membership) <- c("membership")
Female_ASD_SC20_membership <- as.matrix(Female_ASD_SC20$membership)
colnames(Female_ASD_SC20_membership) <- c("membership")
Female_ASD_SC30_membership <- as.matrix(Female_ASD_SC30$membership)
colnames(Female_ASD_SC30_membership) <- c("membership")
Female_ASD_SC40_membership <- as.matrix(Female_ASD_SC40$membership)
colnames(Female_ASD_SC40_membership) <- c("membership")
Female_ASD_SC50_membership <- as.matrix(Female_ASD_SC50$membership)
colnames(Female_ASD_SC50_membership) <- c("membership")

# Female_Normal
Female_Normal_SC10_membership <- as.matrix(Female_Normal_SC10$membership)
colnames(Female_Normal_SC10_membership) <- c("membership")
Female_Normal_SC20_membership <- as.matrix(Female_Normal_SC20$membership)
colnames(Female_Normal_SC20_membership) <- c("membership")
Female_Normal_SC30_membership <- as.matrix(Female_Normal_SC30$membership)
colnames(Female_Normal_SC30_membership) <- c("membership")
Female_Normal_SC40_membership <- as.matrix(Female_Normal_SC40$membership)
colnames(Female_Normal_SC40_membership) <- c("membership")
Female_Normal_SC50_membership <- as.matrix(Female_Normal_SC50$membership)
colnames(Female_Normal_SC50_membership) <- c("membership")

# NeuNRGNII_ASD
NeuNRGNII_ASD_SC10_membership <- as.matrix(NeuNRGNII_ASD_SC10$membership)
colnames(NeuNRGNII_ASD_SC10_membership) <- c("membership")
NeuNRGNII_ASD_SC20_membership <- as.matrix(NeuNRGNII_ASD_SC20$membership)
colnames(NeuNRGNII_ASD_SC20_membership) <- c("membership")
NeuNRGNII_ASD_SC30_membership <- as.matrix(NeuNRGNII_ASD_SC30$membership)
colnames(NeuNRGNII_ASD_SC30_membership) <- c("membership")
NeuNRGNII_ASD_SC40_membership <- as.matrix(NeuNRGNII_ASD_SC40$membership)
colnames(NeuNRGNII_ASD_SC40_membership) <- c("membership")
NeuNRGNII_ASD_SC50_membership <- as.matrix(NeuNRGNII_ASD_SC50$membership)
colnames(NeuNRGNII_ASD_SC50_membership) <- c("membership")

# NeuNRGNII_Normal
NeuNRGNII_Normal_SC10_membership <- as.matrix(NeuNRGNII_Normal_SC10$membership)
colnames(NeuNRGNII_Normal_SC10_membership) <- c("membership")
NeuNRGNII_Normal_SC20_membership <- as.matrix(NeuNRGNII_Normal_SC20$membership)
colnames(NeuNRGNII_Normal_SC20_membership) <- c("membership")
NeuNRGNII_Normal_SC30_membership <- as.matrix(NeuNRGNII_Normal_SC30$membership)
colnames(NeuNRGNII_Normal_SC30_membership) <- c("membership")
NeuNRGNII_Normal_SC40_membership <- as.matrix(NeuNRGNII_Normal_SC40$membership)
colnames(NeuNRGNII_Normal_SC40_membership) <- c("membership")
NeuNRGNII_Normal_SC50_membership <- as.matrix(NeuNRGNII_Normal_SC50$membership)
colnames(NeuNRGNII_Normal_SC50_membership) <- c("membership")

# L56_ASD
L56_ASD_SC10_membership <- as.matrix(L56_ASD_SC10$membership)
colnames(L56_ASD_SC10_membership) <- c("membership")
L56_ASD_SC20_membership <- as.matrix(L56_ASD_SC20$membership)
colnames(L56_ASD_SC20_membership) <- c("membership")
L56_ASD_SC30_membership <- as.matrix(L56_ASD_SC30$membership)
colnames(L56_ASD_SC30_membership) <- c("membership")
L56_ASD_SC40_membership <- as.matrix(L56_ASD_SC40$membership)
colnames(L56_ASD_SC40_membership) <- c("membership")
L56_ASD_SC50_membership <- as.matrix(L56_ASD_SC50$membership)
colnames(L56_ASD_SC50_membership) <- c("membership")

# L56_Normal
L56_Normal_SC10_membership <- as.matrix(L56_Normal_SC10$membership)
colnames(L56_Normal_SC10_membership) <- c("membership")
L56_Normal_SC20_membership <- as.matrix(L56_Normal_SC20$membership)
colnames(L56_Normal_SC20_membership) <- c("membership")
L56_Normal_SC30_membership <- as.matrix(L56_Normal_SC30$membership)
colnames(L56_Normal_SC30_membership) <- c("membership")
L56_Normal_SC40_membership <- as.matrix(L56_Normal_SC40$membership)
colnames(L56_Normal_SC40_membership) <- c("membership")
L56_Normal_SC50_membership <- as.matrix(L56_Normal_SC50$membership)
colnames(L56_Normal_SC50_membership) <- c("membership")

# Oligodendrocytes_ASD
Oligodendrocytes_ASD_SC10_membership <- as.matrix(Oligodendrocytes_ASD_SC10$membership)
colnames(Oligodendrocytes_ASD_SC10_membership) <- c("membership")
Oligodendrocytes_ASD_SC20_membership <- as.matrix(Oligodendrocytes_ASD_SC20$membership)
colnames(Oligodendrocytes_ASD_SC20_membership) <- c("membership")
Oligodendrocytes_ASD_SC30_membership <- as.matrix(Oligodendrocytes_ASD_SC30$membership)
colnames(Oligodendrocytes_ASD_SC30_membership) <- c("membership")
Oligodendrocytes_ASD_SC40_membership <- as.matrix(Oligodendrocytes_ASD_SC40$membership)
colnames(Oligodendrocytes_ASD_SC40_membership) <- c("membership")
Oligodendrocytes_ASD_SC50_membership <- as.matrix(Oligodendrocytes_ASD_SC50$membership)
colnames(Oligodendrocytes_ASD_SC50_membership) <- c("membership")

# Oligodendrocytes_Normal
Oligodendrocytes_Normal_SC10_membership <- as.matrix(Oligodendrocytes_Normal_SC10$membership)
colnames(Oligodendrocytes_Normal_SC10_membership) <- c("membership")
Oligodendrocytes_Normal_SC20_membership <- as.matrix(Oligodendrocytes_Normal_SC20$membership)
colnames(Oligodendrocytes_Normal_SC20_membership) <- c("membership")
Oligodendrocytes_Normal_SC30_membership <- as.matrix(Oligodendrocytes_Normal_SC30$membership)
colnames(Oligodendrocytes_Normal_SC30_membership) <- c("membership")
Oligodendrocytes_Normal_SC40_membership <- as.matrix(Oligodendrocytes_Normal_SC40$membership)
colnames(Oligodendrocytes_Normal_SC40_membership) <- c("membership")
Oligodendrocytes_Normal_SC50_membership <- as.matrix(Oligodendrocytes_Normal_SC50$membership)
colnames(Oligodendrocytes_Normal_SC50_membership) <- c("membership")

# OPC_ASD
OPC_ASD_SC10_membership <- as.matrix(OPC_ASD_SC10$membership)
colnames(OPC_ASD_SC10_membership) <- c("membership")
OPC_ASD_SC20_membership <- as.matrix(OPC_ASD_SC20$membership)
colnames(OPC_ASD_SC20_membership) <- c("membership")
OPC_ASD_SC30_membership <- as.matrix(OPC_ASD_SC30$membership)
colnames(OPC_ASD_SC30_membership) <- c("membership")
OPC_ASD_SC40_membership <- as.matrix(OPC_ASD_SC40$membership)
colnames(OPC_ASD_SC40_membership) <- c("membership")
OPC_ASD_SC50_membership <- as.matrix(OPC_ASD_SC50$membership)
colnames(OPC_ASD_SC50_membership) <- c("membership")

# OPC_Normal
OPC_Normal_SC10_membership <- as.matrix(OPC_Normal_SC10$membership)
colnames(OPC_Normal_SC10_membership) <- c("membership")
OPC_Normal_SC20_membership <- as.matrix(OPC_Normal_SC20$membership)
colnames(OPC_Normal_SC20_membership) <- c("membership")
OPC_Normal_SC30_membership <- as.matrix(OPC_Normal_SC30$membership)
colnames(OPC_Normal_SC30_membership) <- c("membership")
OPC_Normal_SC40_membership <- as.matrix(OPC_Normal_SC40$membership)
colnames(OPC_Normal_SC40_membership) <- c("membership")
OPC_Normal_SC50_membership <- as.matrix(OPC_Normal_SC50$membership)
colnames(OPC_Normal_SC50_membership) <- c("membership")

# ASTFB_ASD
ASTFB_ASD_SC10_membership <- as.matrix(ASTFB_ASD_SC10$membership)
colnames(ASTFB_ASD_SC10_membership) <- c("membership")
ASTFB_ASD_SC20_membership <- as.matrix(ASTFB_ASD_SC20$membership)
colnames(ASTFB_ASD_SC20_membership) <- c("membership")
ASTFB_ASD_SC30_membership <- as.matrix(ASTFB_ASD_SC30$membership)
colnames(ASTFB_ASD_SC30_membership) <- c("membership")
ASTFB_ASD_SC40_membership <- as.matrix(ASTFB_ASD_SC40$membership)
colnames(ASTFB_ASD_SC40_membership) <- c("membership")
ASTFB_ASD_SC50_membership <- as.matrix(ASTFB_ASD_SC50$membership)
colnames(ASTFB_ASD_SC50_membership) <- c("membership")

# ASTFB_Normal
ASTFB_Normal_SC10_membership <- as.matrix(ASTFB_Normal_SC10$membership)
colnames(ASTFB_Normal_SC10_membership) <- c("membership")
ASTFB_Normal_SC20_membership <- as.matrix(ASTFB_Normal_SC20$membership)
colnames(ASTFB_Normal_SC20_membership) <- c("membership")
ASTFB_Normal_SC30_membership <- as.matrix(ASTFB_Normal_SC30$membership)
colnames(ASTFB_Normal_SC30_membership) <- c("membership")
ASTFB_Normal_SC40_membership <- as.matrix(ASTFB_Normal_SC40$membership)
colnames(ASTFB_Normal_SC40_membership) <- c("membership")
ASTFB_Normal_SC50_membership <- as.matrix(ASTFB_Normal_SC50$membership)
colnames(ASTFB_Normal_SC50_membership) <- c("membership")

# Endothelial_ASD
Endothelial_ASD_SC10_membership <- as.matrix(Endothelial_ASD_SC10$membership)
colnames(Endothelial_ASD_SC10_membership) <- c("membership")
Endothelial_ASD_SC20_membership <- as.matrix(Endothelial_ASD_SC20$membership)
colnames(Endothelial_ASD_SC20_membership) <- c("membership")
Endothelial_ASD_SC30_membership <- as.matrix(Endothelial_ASD_SC30$membership)
colnames(Endothelial_ASD_SC30_membership) <- c("membership")
Endothelial_ASD_SC40_membership <- as.matrix(Endothelial_ASD_SC40$membership)
colnames(Endothelial_ASD_SC40_membership) <- c("membership")
Endothelial_ASD_SC50_membership <- as.matrix(Endothelial_ASD_SC50$membership)
colnames(Endothelial_ASD_SC50_membership) <- c("membership")

# Endothelial_Normal
Endothelial_Normal_SC10_membership <- as.matrix(Endothelial_Normal_SC10$membership)
colnames(Endothelial_Normal_SC10_membership) <- c("membership")
Endothelial_Normal_SC20_membership <- as.matrix(Endothelial_Normal_SC20$membership)
colnames(Endothelial_Normal_SC20_membership) <- c("membership")
Endothelial_Normal_SC30_membership <- as.matrix(Endothelial_Normal_SC30$membership)
colnames(Endothelial_Normal_SC30_membership) <- c("membership")
Endothelial_Normal_SC40_membership <- as.matrix(Endothelial_Normal_SC40$membership)
colnames(Endothelial_Normal_SC40_membership) <- c("membership")
Endothelial_Normal_SC50_membership <- as.matrix(Endothelial_Normal_SC50$membership)
colnames(Endothelial_Normal_SC50_membership) <- c("membership")

# Microglia_ASD
Microglia_ASD_SC10_membership <- as.matrix(Microglia_ASD_SC10$membership)
colnames(Microglia_ASD_SC10_membership) <- c("membership")
Microglia_ASD_SC20_membership <- as.matrix(Microglia_ASD_SC20$membership)
colnames(Microglia_ASD_SC20_membership) <- c("membership")
Microglia_ASD_SC30_membership <- as.matrix(Microglia_ASD_SC30$membership)
colnames(Microglia_ASD_SC30_membership) <- c("membership")
Microglia_ASD_SC40_membership <- as.matrix(Microglia_ASD_SC40$membership)
colnames(Microglia_ASD_SC40_membership) <- c("membership")
Microglia_ASD_SC50_membership <- as.matrix(Microglia_ASD_SC50$membership)
colnames(Microglia_ASD_SC50_membership) <- c("membership")

# Microglia_Normal
Microglia_Normal_SC10_membership <- as.matrix(Microglia_Normal_SC10$membership)
colnames(Microglia_Normal_SC10_membership) <- c("membership")
Microglia_Normal_SC20_membership <- as.matrix(Microglia_Normal_SC20$membership)
colnames(Microglia_Normal_SC20_membership) <- c("membership")
Microglia_Normal_SC30_membership <- as.matrix(Microglia_Normal_SC30$membership)
colnames(Microglia_Normal_SC30_membership) <- c("membership")
Microglia_Normal_SC40_membership <- as.matrix(Microglia_Normal_SC40$membership)
colnames(Microglia_Normal_SC40_membership) <- c("membership")
Microglia_Normal_SC50_membership <- as.matrix(Microglia_Normal_SC50$membership)
colnames(Microglia_Normal_SC50_membership) <- c("membership")

# NeuNRGNI_ASD
NeuNRGNI_ASD_SC10_membership <- as.matrix(NeuNRGNI_ASD_SC10$membership)
colnames(NeuNRGNI_ASD_SC10_membership) <- c("membership")
NeuNRGNI_ASD_SC20_membership <- as.matrix(NeuNRGNI_ASD_SC20$membership)
colnames(NeuNRGNI_ASD_SC20_membership) <- c("membership")
NeuNRGNI_ASD_SC30_membership <- as.matrix(NeuNRGNI_ASD_SC30$membership)
colnames(NeuNRGNI_ASD_SC30_membership) <- c("membership")
NeuNRGNI_ASD_SC40_membership <- as.matrix(NeuNRGNI_ASD_SC40$membership)
colnames(NeuNRGNI_ASD_SC40_membership) <- c("membership")
NeuNRGNI_ASD_SC50_membership <- as.matrix(NeuNRGNI_ASD_SC50$membership)
colnames(NeuNRGNI_ASD_SC50_membership) <- c("membership")

# NeuNRGNI_Normal
NeuNRGNI_Normal_SC10_membership <- as.matrix(NeuNRGNI_Normal_SC10$membership)
colnames(NeuNRGNI_Normal_SC10_membership) <- c("membership")
NeuNRGNI_Normal_SC20_membership <- as.matrix(NeuNRGNI_Normal_SC20$membership)
colnames(NeuNRGNI_Normal_SC20_membership) <- c("membership")
NeuNRGNI_Normal_SC30_membership <- as.matrix(NeuNRGNI_Normal_SC30$membership)
colnames(NeuNRGNI_Normal_SC30_membership) <- c("membership")
NeuNRGNI_Normal_SC40_membership <- as.matrix(NeuNRGNI_Normal_SC40$membership)
colnames(NeuNRGNI_Normal_SC40_membership) <- c("membership")
NeuNRGNI_Normal_SC50_membership <- as.matrix(NeuNRGNI_Normal_SC50$membership)
colnames(NeuNRGNI_Normal_SC50_membership) <- c("membership")

# INVIP_ASD
INVIP_ASD_SC10_membership <- as.matrix(INVIP_ASD_SC10$membership)
colnames(INVIP_ASD_SC10_membership) <- c("membership")
INVIP_ASD_SC20_membership <- as.matrix(INVIP_ASD_SC20$membership)
colnames(INVIP_ASD_SC20_membership) <- c("membership")
INVIP_ASD_SC30_membership <- as.matrix(INVIP_ASD_SC30$membership)
colnames(INVIP_ASD_SC30_membership) <- c("membership")
INVIP_ASD_SC40_membership <- as.matrix(INVIP_ASD_SC40$membership)
colnames(INVIP_ASD_SC40_membership) <- c("membership")
INVIP_ASD_SC50_membership <- as.matrix(INVIP_ASD_SC50$membership)
colnames(INVIP_ASD_SC50_membership) <- c("membership")

# INVIP_Normal
INVIP_Normal_SC10_membership <- as.matrix(INVIP_Normal_SC10$membership)
colnames(INVIP_Normal_SC10_membership) <- c("membership")
INVIP_Normal_SC20_membership <- as.matrix(INVIP_Normal_SC20$membership)
colnames(INVIP_Normal_SC20_membership) <- c("membership")
INVIP_Normal_SC30_membership <- as.matrix(INVIP_Normal_SC30$membership)
colnames(INVIP_Normal_SC30_membership) <- c("membership")
INVIP_Normal_SC40_membership <- as.matrix(INVIP_Normal_SC40$membership)
colnames(INVIP_Normal_SC40_membership) <- c("membership")
INVIP_Normal_SC50_membership <- as.matrix(INVIP_Normal_SC50$membership)
colnames(INVIP_Normal_SC50_membership) <- c("membership")

# L56CC_ASD
L56CC_ASD_SC10_membership <- as.matrix(L56CC_ASD_SC10$membership)
colnames(L56CC_ASD_SC10_membership) <- c("membership")
L56CC_ASD_SC20_membership <- as.matrix(L56CC_ASD_SC20$membership)
colnames(L56CC_ASD_SC20_membership) <- c("membership")
L56CC_ASD_SC30_membership <- as.matrix(L56CC_ASD_SC30$membership)
colnames(L56CC_ASD_SC30_membership) <- c("membership")
L56CC_ASD_SC40_membership <- as.matrix(L56CC_ASD_SC40$membership)
colnames(L56CC_ASD_SC40_membership) <- c("membership")
L56CC_ASD_SC50_membership <- as.matrix(L56CC_ASD_SC50$membership)
colnames(L56CC_ASD_SC50_membership) <- c("membership")

# L56CC_Normal
L56CC_Normal_SC10_membership <- as.matrix(L56CC_Normal_SC10$membership)
colnames(L56CC_Normal_SC10_membership) <- c("membership")
L56CC_Normal_SC20_membership <- as.matrix(L56CC_Normal_SC20$membership)
colnames(L56CC_Normal_SC20_membership) <- c("membership")
L56CC_Normal_SC30_membership <- as.matrix(L56CC_Normal_SC30$membership)
colnames(L56CC_Normal_SC30_membership) <- c("membership")
L56CC_Normal_SC40_membership <- as.matrix(L56CC_Normal_SC40$membership)
colnames(L56CC_Normal_SC40_membership) <- c("membership")
L56CC_Normal_SC50_membership <- as.matrix(L56CC_Normal_SC50$membership)
colnames(L56CC_Normal_SC50_membership) <- c("membership")

# INSV2C_ASD
INSV2C_ASD_SC10_membership <- as.matrix(INSV2C_ASD_SC10$membership)
colnames(INSV2C_ASD_SC10_membership) <- c("membership")
INSV2C_ASD_SC20_membership <- as.matrix(INSV2C_ASD_SC20$membership)
colnames(INSV2C_ASD_SC20_membership) <- c("membership")
INSV2C_ASD_SC30_membership <- as.matrix(INSV2C_ASD_SC30$membership)
colnames(INSV2C_ASD_SC30_membership) <- c("membership")
INSV2C_ASD_SC40_membership <- as.matrix(INSV2C_ASD_SC40$membership)
colnames(INSV2C_ASD_SC40_membership) <- c("membership")
INSV2C_ASD_SC50_membership <- as.matrix(INSV2C_ASD_SC50$membership)
colnames(INSV2C_ASD_SC50_membership) <- c("membership")

# INSV2C_Normal
INSV2C_Normal_SC10_membership <- as.matrix(INSV2C_Normal_SC10$membership)
colnames(INSV2C_Normal_SC10_membership) <- c("membership")
INSV2C_Normal_SC20_membership <- as.matrix(INSV2C_Normal_SC20$membership)
colnames(INSV2C_Normal_SC20_membership) <- c("membership")
INSV2C_Normal_SC30_membership <- as.matrix(INSV2C_Normal_SC30$membership)
colnames(INSV2C_Normal_SC30_membership) <- c("membership")
INSV2C_Normal_SC40_membership <- as.matrix(INSV2C_Normal_SC40$membership)
colnames(INSV2C_Normal_SC40_membership) <- c("membership")
INSV2C_Normal_SC50_membership <- as.matrix(INSV2C_Normal_SC50$membership)
colnames(INSV2C_Normal_SC50_membership) <- c("membership")

# L23_ASD
L23_ASD_SC10_membership <- as.matrix(L23_ASD_SC10$membership)
colnames(L23_ASD_SC10_membership) <- c("membership")
L23_ASD_SC20_membership <- as.matrix(L23_ASD_SC20$membership)
colnames(L23_ASD_SC20_membership) <- c("membership")
L23_ASD_SC30_membership <- as.matrix(L23_ASD_SC30$membership)
colnames(L23_ASD_SC30_membership) <- c("membership")
L23_ASD_SC40_membership <- as.matrix(L23_ASD_SC40$membership)
colnames(L23_ASD_SC40_membership) <- c("membership")
L23_ASD_SC50_membership <- as.matrix(L23_ASD_SC50$membership)
colnames(L23_ASD_SC50_membership) <- c("membership")

# L23_Normal
L23_Normal_SC10_membership <- as.matrix(L23_Normal_SC10$membership)
colnames(L23_Normal_SC10_membership) <- c("membership")
L23_Normal_SC20_membership <- as.matrix(L23_Normal_SC20$membership)
colnames(L23_Normal_SC20_membership) <- c("membership")
L23_Normal_SC30_membership <- as.matrix(L23_Normal_SC30$membership)
colnames(L23_Normal_SC30_membership) <- c("membership")
L23_Normal_SC40_membership <- as.matrix(L23_Normal_SC40$membership)
colnames(L23_Normal_SC40_membership) <- c("membership")
L23_Normal_SC50_membership <- as.matrix(L23_Normal_SC50$membership)
colnames(L23_Normal_SC50_membership) <- c("membership")

# INPV_ASD
INPV_ASD_SC10_membership <- as.matrix(INPV_ASD_SC10$membership)
colnames(INPV_ASD_SC10_membership) <- c("membership")
INPV_ASD_SC20_membership <- as.matrix(INPV_ASD_SC20$membership)
colnames(INPV_ASD_SC20_membership) <- c("membership")
INPV_ASD_SC30_membership <- as.matrix(INPV_ASD_SC30$membership)
colnames(INPV_ASD_SC30_membership) <- c("membership")
INPV_ASD_SC40_membership <- as.matrix(INPV_ASD_SC40$membership)
colnames(INPV_ASD_SC40_membership) <- c("membership")
INPV_ASD_SC50_membership <- as.matrix(INPV_ASD_SC50$membership)
colnames(INPV_ASD_SC50_membership) <- c("membership")

# INPV_Normal
INPV_Normal_SC10_membership <- as.matrix(INPV_Normal_SC10$membership)
colnames(INPV_Normal_SC10_membership) <- c("membership")
INPV_Normal_SC20_membership <- as.matrix(INPV_Normal_SC20$membership)
colnames(INPV_Normal_SC20_membership) <- c("membership")
INPV_Normal_SC30_membership <- as.matrix(INPV_Normal_SC30$membership)
colnames(INPV_Normal_SC30_membership) <- c("membership")
INPV_Normal_SC40_membership <- as.matrix(INPV_Normal_SC40$membership)
colnames(INPV_Normal_SC40_membership) <- c("membership")
INPV_Normal_SC50_membership <- as.matrix(INPV_Normal_SC50$membership)
colnames(INPV_Normal_SC50_membership) <- c("membership")

# L4_ASD
L4_ASD_SC10_membership <- as.matrix(L4_ASD_SC10$membership)
colnames(L4_ASD_SC10_membership) <- c("membership")
L4_ASD_SC20_membership <- as.matrix(L4_ASD_SC20$membership)
colnames(L4_ASD_SC20_membership) <- c("membership")
L4_ASD_SC30_membership <- as.matrix(L4_ASD_SC30$membership)
colnames(L4_ASD_SC30_membership) <- c("membership")
L4_ASD_SC40_membership <- as.matrix(L4_ASD_SC40$membership)
colnames(L4_ASD_SC40_membership) <- c("membership")
L4_ASD_SC50_membership <- as.matrix(L4_ASD_SC50$membership)
colnames(L4_ASD_SC50_membership) <- c("membership")

# L4_Normal
L4_Normal_SC10_membership <- as.matrix(L4_Normal_SC10$membership)
colnames(L4_Normal_SC10_membership) <- c("membership")
L4_Normal_SC20_membership <- as.matrix(L4_Normal_SC20$membership)
colnames(L4_Normal_SC20_membership) <- c("membership")
L4_Normal_SC30_membership <- as.matrix(L4_Normal_SC30$membership)
colnames(L4_Normal_SC30_membership) <- c("membership")
L4_Normal_SC40_membership <- as.matrix(L4_Normal_SC40$membership)
colnames(L4_Normal_SC40_membership) <- c("membership")
L4_Normal_SC50_membership <- as.matrix(L4_Normal_SC50$membership)
colnames(L4_Normal_SC50_membership) <- c("membership")

# INSST_ASD
INSST_ASD_SC10_membership <- as.matrix(INSST_ASD_SC10$membership)
colnames(INSST_ASD_SC10_membership) <- c("membership")
INSST_ASD_SC20_membership <- as.matrix(INSST_ASD_SC20$membership)
colnames(INSST_ASD_SC20_membership) <- c("membership")
INSST_ASD_SC30_membership <- as.matrix(INSST_ASD_SC30$membership)
colnames(INSST_ASD_SC30_membership) <- c("membership")
INSST_ASD_SC40_membership <- as.matrix(INSST_ASD_SC40$membership)
colnames(INSST_ASD_SC40_membership) <- c("membership")
INSST_ASD_SC50_membership <- as.matrix(INSST_ASD_SC50$membership)
colnames(INSST_ASD_SC50_membership) <- c("membership")

# INSST_Normal
INSST_Normal_SC10_membership <- as.matrix(INSST_Normal_SC10$membership)
colnames(INSST_Normal_SC10_membership) <- c("membership")
INSST_Normal_SC20_membership <- as.matrix(INSST_Normal_SC20$membership)
colnames(INSST_Normal_SC20_membership) <- c("membership")
INSST_Normal_SC30_membership <- as.matrix(INSST_Normal_SC30$membership)
colnames(INSST_Normal_SC30_membership) <- c("membership")
INSST_Normal_SC40_membership <- as.matrix(INSST_Normal_SC40$membership)
colnames(INSST_Normal_SC40_membership) <- c("membership")
INSST_Normal_SC50_membership <- as.matrix(INSST_Normal_SC50$membership)
colnames(INSST_Normal_SC50_membership) <- c("membership")

# Neumat_ASD
Neumat_ASD_SC10_membership <- as.matrix(Neumat_ASD_SC10$membership)
colnames(Neumat_ASD_SC10_membership) <- c("membership")
Neumat_ASD_SC20_membership <- as.matrix(Neumat_ASD_SC20$membership)
colnames(Neumat_ASD_SC20_membership) <- c("membership")
Neumat_ASD_SC30_membership <- as.matrix(Neumat_ASD_SC30$membership)
colnames(Neumat_ASD_SC30_membership) <- c("membership")
Neumat_ASD_SC40_membership <- as.matrix(Neumat_ASD_SC40$membership)
colnames(Neumat_ASD_SC40_membership) <- c("membership")
Neumat_ASD_SC50_membership <- as.matrix(Neumat_ASD_SC50$membership)
colnames(Neumat_ASD_SC50_membership) <- c("membership")

# Neumat_Normal
Neumat_Normal_SC10_membership <- as.matrix(Neumat_Normal_SC10$membership)
colnames(Neumat_Normal_SC10_membership) <- c("membership")
Neumat_Normal_SC20_membership <- as.matrix(Neumat_Normal_SC20$membership)
colnames(Neumat_Normal_SC20_membership) <- c("membership")
Neumat_Normal_SC30_membership <- as.matrix(Neumat_Normal_SC30$membership)
colnames(Neumat_Normal_SC30_membership) <- c("membership")
Neumat_Normal_SC40_membership <- as.matrix(Neumat_Normal_SC40$membership)
colnames(Neumat_Normal_SC40_membership) <- c("membership")
Neumat_Normal_SC50_membership <- as.matrix(Neumat_Normal_SC50$membership)
colnames(Neumat_Normal_SC50_membership) <- c("membership")

# ASTPP_ASD
ASTPP_ASD_SC10_membership <- as.matrix(ASTPP_ASD_SC10$membership)
colnames(ASTPP_ASD_SC10_membership) <- c("membership")
ASTPP_ASD_SC20_membership <- as.matrix(ASTPP_ASD_SC20$membership)
colnames(ASTPP_ASD_SC20_membership) <- c("membership")
ASTPP_ASD_SC30_membership <- as.matrix(ASTPP_ASD_SC30$membership)
colnames(ASTPP_ASD_SC30_membership) <- c("membership")
ASTPP_ASD_SC40_membership <- as.matrix(ASTPP_ASD_SC40$membership)
colnames(ASTPP_ASD_SC40_membership) <- c("membership")
ASTPP_ASD_SC50_membership <- as.matrix(ASTPP_ASD_SC50$membership)
colnames(ASTPP_ASD_SC50_membership) <- c("membership")

# ASTPP_Normal
ASTPP_Normal_SC10_membership <- as.matrix(ASTPP_Normal_SC10$membership)
colnames(ASTPP_Normal_SC10_membership) <- c("membership")
ASTPP_Normal_SC20_membership <- as.matrix(ASTPP_Normal_SC20$membership)
colnames(ASTPP_Normal_SC20_membership) <- c("membership")
ASTPP_Normal_SC30_membership <- as.matrix(ASTPP_Normal_SC30$membership)
colnames(ASTPP_Normal_SC30_membership) <- c("membership")
ASTPP_Normal_SC40_membership <- as.matrix(ASTPP_Normal_SC40$membership)
colnames(ASTPP_Normal_SC40_membership) <- c("membership")
ASTPP_Normal_SC50_membership <- as.matrix(ASTPP_Normal_SC50$membership)
colnames(ASTPP_Normal_SC50_membership) <- c("membership")

ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), dims = 1:30)
Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), dims = 1:30)
ACC_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), dims = 1:30)
ACC_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), dims = 1:30)
PFC_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), dims = 1:30)
PFC_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), dims = 1:30)
Lower18_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), dims = 1:30)
Lower18_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), dims = 1:30)
Larger18_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), dims = 1:30)
Larger18_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), dims = 1:30)
Male_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), dims = 1:30)
Male_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), dims = 1:30)
Female_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), dims = 1:30)
Female_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), dims = 1:30)
NeuNRGNII_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), dims = 1:30)
NeuNRGNII_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), dims = 1:30)
L56_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), dims = 1:30)
L56_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), dims = 1:30)
Oligodendrocytes_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), dims = 1:30)
Oligodendrocytes_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), dims = 1:30)
OPC_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), dims = 1:30)
OPC_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), dims = 1:30)
ASTFB_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), dims = 1:30)
ASTFB_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), dims = 1:30)
Endothelial_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), dims = 1:30)
Endothelial_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), dims = 1:30)
Microglia_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), dims = 1:30)
Microglia_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), dims = 1:30)
NeuNRGNI_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), dims = 1:30)
NeuNRGNI_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), dims = 1:30)
INVIP_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), dims = 1:30)
INVIP_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), dims = 1:30)
L56CC_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), dims = 1:30)
L56CC_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), dims = 1:30)
INSV2C_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), dims = 1:30)
INSV2C_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), dims = 1:30)
L23_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), dims = 1:30)
L23_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), dims = 1:30)
INPV_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), dims = 1:30)
INPV_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), dims = 1:30)
L4_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), dims = 1:30)
L4_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), dims = 1:30)
INSST_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), dims = 1:30)
INSST_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), dims = 1:30)
Neumat_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), dims = 1:30)
Neumat_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), dims = 1:30)
ASTPP_ASD_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), dims = 1:30)
ASTPP_Normal_diffusion_comp <- get_diffusion_comp(sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), dims = 1:30)

# ASD
ASD_SC10_compactness <- mc_compactness(cell.membership = ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)),
                                      sc.reduction = ASD_diffusion_comp)
ASD_SC20_compactness <- mc_compactness(cell.membership = ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)),
                                      sc.reduction = ASD_diffusion_comp)
ASD_SC30_compactness <- mc_compactness(cell.membership = ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)),
                                      sc.reduction = ASD_diffusion_comp)
ASD_SC40_compactness <- mc_compactness(cell.membership = ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)),
                                      sc.reduction = ASD_diffusion_comp)
ASD_SC50_compactness <- mc_compactness(cell.membership = ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)),
                                      sc.reduction = ASD_diffusion_comp)

# Normal
Normal_SC10_compactness <- mc_compactness(cell.membership = Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)),
                                      sc.reduction = Normal_diffusion_comp)
Normal_SC20_compactness <- mc_compactness(cell.membership = Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)),
                                      sc.reduction = Normal_diffusion_comp)
Normal_SC30_compactness <- mc_compactness(cell.membership = Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)),
                                      sc.reduction = Normal_diffusion_comp)
Normal_SC40_compactness <- mc_compactness(cell.membership = Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)),
                                      sc.reduction = Normal_diffusion_comp)
Normal_SC50_compactness <- mc_compactness(cell.membership = Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)),
                                      sc.reduction = Normal_diffusion_comp)

# ACC_ASD
ACC_ASD_SC10_compactness <- mc_compactness(cell.membership = ACC_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)),
                                      sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC20_compactness <- mc_compactness(cell.membership = ACC_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)),
                                      sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC30_compactness <- mc_compactness(cell.membership = ACC_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)),
                                      sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC40_compactness <- mc_compactness(cell.membership = ACC_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)),
                                      sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC50_compactness <- mc_compactness(cell.membership = ACC_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)),
                                      sc.reduction = ACC_ASD_diffusion_comp)

# ACC_Normal
ACC_Normal_SC10_compactness <- mc_compactness(cell.membership = ACC_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)),
                                      sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC20_compactness <- mc_compactness(cell.membership = ACC_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)),
                                      sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC30_compactness <- mc_compactness(cell.membership = ACC_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)),
                                      sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC40_compactness <- mc_compactness(cell.membership = ACC_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)),
                                      sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC50_compactness <- mc_compactness(cell.membership = ACC_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)),
                                      sc.reduction = ACC_Normal_diffusion_comp)

# PFC_ASD
PFC_ASD_SC10_compactness <- mc_compactness(cell.membership = PFC_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)),
                                      sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC20_compactness <- mc_compactness(cell.membership = PFC_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)),
                                      sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC30_compactness <- mc_compactness(cell.membership = PFC_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)),
                                      sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC40_compactness <- mc_compactness(cell.membership = PFC_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)),
                                      sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC50_compactness <- mc_compactness(cell.membership = PFC_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)),
                                      sc.reduction = PFC_ASD_diffusion_comp)

# PFC_Normal
PFC_Normal_SC10_compactness <- mc_compactness(cell.membership = PFC_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)),
                                      sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC20_compactness <- mc_compactness(cell.membership = PFC_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)),
                                      sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC30_compactness <- mc_compactness(cell.membership = PFC_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)),
                                      sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC40_compactness <- mc_compactness(cell.membership = PFC_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)),
                                      sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC50_compactness <- mc_compactness(cell.membership = PFC_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)),
                                      sc.reduction = PFC_Normal_diffusion_comp)

# Lower18_ASD
Lower18_ASD_SC10_compactness <- mc_compactness(cell.membership = Lower18_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)),
                                      sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC20_compactness <- mc_compactness(cell.membership = Lower18_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)),
                                      sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC30_compactness <- mc_compactness(cell.membership = Lower18_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)),
                                      sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC40_compactness <- mc_compactness(cell.membership = Lower18_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)),
                                      sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC50_compactness <- mc_compactness(cell.membership = Lower18_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)),
                                      sc.reduction = Lower18_ASD_diffusion_comp)

# Lower18_Normal
Lower18_Normal_SC10_compactness <- mc_compactness(cell.membership = Lower18_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)),
                                      sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC20_compactness <- mc_compactness(cell.membership = Lower18_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)),
                                      sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC30_compactness <- mc_compactness(cell.membership = Lower18_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)),
                                      sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC40_compactness <- mc_compactness(cell.membership = Lower18_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)),
                                      sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC50_compactness <- mc_compactness(cell.membership = Lower18_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)),
                                      sc.reduction = Lower18_Normal_diffusion_comp)

# Larger18_ASD
Larger18_ASD_SC10_compactness <- mc_compactness(cell.membership = Larger18_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)),
                                      sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC20_compactness <- mc_compactness(cell.membership = Larger18_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)),
                                      sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC30_compactness <- mc_compactness(cell.membership = Larger18_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)),
                                      sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC40_compactness <- mc_compactness(cell.membership = Larger18_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)),
                                      sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC50_compactness <- mc_compactness(cell.membership = Larger18_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)),
                                      sc.reduction = Larger18_ASD_diffusion_comp)

# Larger18_Normal
Larger18_Normal_SC10_compactness <- mc_compactness(cell.membership = Larger18_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)),
                                      sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC20_compactness <- mc_compactness(cell.membership = Larger18_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)),
                                      sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC30_compactness <- mc_compactness(cell.membership = Larger18_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)),
                                      sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC40_compactness <- mc_compactness(cell.membership = Larger18_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)),
                                      sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC50_compactness <- mc_compactness(cell.membership = Larger18_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)),
                                      sc.reduction = Larger18_Normal_diffusion_comp)

# Male_ASD
Male_ASD_SC10_compactness <- mc_compactness(cell.membership = Male_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)),
                                      sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC20_compactness <- mc_compactness(cell.membership = Male_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)),
                                      sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC30_compactness <- mc_compactness(cell.membership = Male_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)),
                                      sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC40_compactness <- mc_compactness(cell.membership = Male_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)),
                                      sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC50_compactness <- mc_compactness(cell.membership = Male_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)),
                                      sc.reduction = Male_ASD_diffusion_comp)

# Male_Normal
Male_Normal_SC10_compactness <- mc_compactness(cell.membership = Male_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)),
                                      sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC20_compactness <- mc_compactness(cell.membership = Male_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)),
                                      sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC30_compactness <- mc_compactness(cell.membership = Male_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)),
                                      sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC40_compactness <- mc_compactness(cell.membership = Male_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)),
                                      sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC50_compactness <- mc_compactness(cell.membership = Male_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)),
                                      sc.reduction = Male_Normal_diffusion_comp)

# Female_ASD
Female_ASD_SC10_compactness <- mc_compactness(cell.membership = Female_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)),
                                      sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC20_compactness <- mc_compactness(cell.membership = Female_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)),
                                      sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC30_compactness <- mc_compactness(cell.membership = Female_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)),
                                      sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC40_compactness <- mc_compactness(cell.membership = Female_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)),
                                      sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC50_compactness <- mc_compactness(cell.membership = Female_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)),
                                      sc.reduction = Female_ASD_diffusion_comp)

# Female_Normal
Female_Normal_SC10_compactness <- mc_compactness(cell.membership = Female_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)),
                                      sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC20_compactness <- mc_compactness(cell.membership = Female_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)),
                                      sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC30_compactness <- mc_compactness(cell.membership = Female_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)),
                                      sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC40_compactness <- mc_compactness(cell.membership = Female_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)),
                                      sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC50_compactness <- mc_compactness(cell.membership = Female_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)),
                                      sc.reduction = Female_Normal_diffusion_comp)

# NeuNRGNII_ASD
NeuNRGNII_ASD_SC10_compactness <- mc_compactness(cell.membership = NeuNRGNII_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC20_compactness <- mc_compactness(cell.membership = NeuNRGNII_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC30_compactness <- mc_compactness(cell.membership = NeuNRGNII_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC40_compactness <- mc_compactness(cell.membership = NeuNRGNII_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC50_compactness <- mc_compactness(cell.membership = NeuNRGNII_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_ASD_diffusion_comp)

# NeuNRGNII_Normal
NeuNRGNII_Normal_SC10_compactness <- mc_compactness(cell.membership = NeuNRGNII_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC20_compactness <- mc_compactness(cell.membership = NeuNRGNII_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC30_compactness <- mc_compactness(cell.membership = NeuNRGNII_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC40_compactness <- mc_compactness(cell.membership = NeuNRGNII_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC50_compactness <- mc_compactness(cell.membership = NeuNRGNII_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNII_Normal_diffusion_comp)

# L56_ASD
L56_ASD_SC10_compactness <- mc_compactness(cell.membership = L56_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC20_compactness <- mc_compactness(cell.membership = L56_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC30_compactness <- mc_compactness(cell.membership = L56_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC40_compactness <- mc_compactness(cell.membership = L56_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC50_compactness <- mc_compactness(cell.membership = L56_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56_ASD_diffusion_comp)

# L56_Normal
L56_Normal_SC10_compactness <- mc_compactness(cell.membership = L56_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])),
                                      sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC20_compactness <- mc_compactness(cell.membership = L56_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])),
                                      sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC30_compactness <- mc_compactness(cell.membership = L56_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])),
                                      sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC40_compactness <- mc_compactness(cell.membership = L56_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])),
                                      sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC50_compactness <- mc_compactness(cell.membership = L56_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])),
                                      sc.reduction = L56_Normal_diffusion_comp)

# Oligodendrocytes_ASD
Oligodendrocytes_ASD_SC10_compactness <- mc_compactness(cell.membership = Oligodendrocytes_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC20_compactness <- mc_compactness(cell.membership = Oligodendrocytes_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC30_compactness <- mc_compactness(cell.membership = Oligodendrocytes_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC40_compactness <- mc_compactness(cell.membership = Oligodendrocytes_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC50_compactness <- mc_compactness(cell.membership = Oligodendrocytes_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_ASD_diffusion_comp)

# Oligodendrocytes_Normal
Oligodendrocytes_Normal_SC10_compactness <- mc_compactness(cell.membership = Oligodendrocytes_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC20_compactness <- mc_compactness(cell.membership = Oligodendrocytes_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC30_compactness <- mc_compactness(cell.membership = Oligodendrocytes_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC40_compactness <- mc_compactness(cell.membership = Oligodendrocytes_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC50_compactness <- mc_compactness(cell.membership = Oligodendrocytes_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])),
                                      sc.reduction = Oligodendrocytes_Normal_diffusion_comp)

# OPC_ASD
OPC_ASD_SC10_compactness <- mc_compactness(cell.membership = OPC_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])),
                                      sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC20_compactness <- mc_compactness(cell.membership = OPC_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])),
                                      sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC30_compactness <- mc_compactness(cell.membership = OPC_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])),
                                      sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC40_compactness <- mc_compactness(cell.membership = OPC_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])),
                                      sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC50_compactness <- mc_compactness(cell.membership = OPC_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])),
                                      sc.reduction = OPC_ASD_diffusion_comp)

# OPC_Normal
OPC_Normal_SC10_compactness <- mc_compactness(cell.membership = OPC_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])),
                                      sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC20_compactness <- mc_compactness(cell.membership = OPC_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])),
                                      sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC30_compactness <- mc_compactness(cell.membership = OPC_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])),
                                      sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC40_compactness <- mc_compactness(cell.membership = OPC_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])),
                                      sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC50_compactness <- mc_compactness(cell.membership = OPC_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])),
                                      sc.reduction = OPC_Normal_diffusion_comp)

# ASTFB_ASD
ASTFB_ASD_SC10_compactness <- mc_compactness(cell.membership = ASTFB_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC20_compactness <- mc_compactness(cell.membership = ASTFB_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC30_compactness <- mc_compactness(cell.membership = ASTFB_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC40_compactness <- mc_compactness(cell.membership = ASTFB_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC50_compactness <- mc_compactness(cell.membership = ASTFB_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_ASD_diffusion_comp)

# ASTFB_Normal
ASTFB_Normal_SC10_compactness <- mc_compactness(cell.membership = ASTFB_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC20_compactness <- mc_compactness(cell.membership = ASTFB_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC30_compactness <- mc_compactness(cell.membership = ASTFB_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC40_compactness <- mc_compactness(cell.membership = ASTFB_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC50_compactness <- mc_compactness(cell.membership = ASTFB_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])),
                                      sc.reduction = ASTFB_Normal_diffusion_comp)

# Endothelial_ASD
Endothelial_ASD_SC10_compactness <- mc_compactness(cell.membership = Endothelial_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC20_compactness <- mc_compactness(cell.membership = Endothelial_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC30_compactness <- mc_compactness(cell.membership = Endothelial_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC40_compactness <- mc_compactness(cell.membership = Endothelial_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC50_compactness <- mc_compactness(cell.membership = Endothelial_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_ASD_diffusion_comp)

# Endothelial_Normal
Endothelial_Normal_SC10_compactness <- mc_compactness(cell.membership = Endothelial_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC20_compactness <- mc_compactness(cell.membership = Endothelial_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC30_compactness <- mc_compactness(cell.membership = Endothelial_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC40_compactness <- mc_compactness(cell.membership = Endothelial_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC50_compactness <- mc_compactness(cell.membership = Endothelial_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])),
                                      sc.reduction = Endothelial_Normal_diffusion_comp)

# Microglia_ASD
Microglia_ASD_SC10_compactness <- mc_compactness(cell.membership = Microglia_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])),
                                      sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC20_compactness <- mc_compactness(cell.membership = Microglia_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])),
                                      sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC30_compactness <- mc_compactness(cell.membership = Microglia_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])),
                                      sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC40_compactness <- mc_compactness(cell.membership = Microglia_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])),
                                      sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC50_compactness <- mc_compactness(cell.membership = Microglia_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])),
                                      sc.reduction = Microglia_ASD_diffusion_comp)

# Microglia_Normal
Microglia_Normal_SC10_compactness <- mc_compactness(cell.membership = Microglia_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])),
                                      sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC20_compactness <- mc_compactness(cell.membership = Microglia_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])),
                                      sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC30_compactness <- mc_compactness(cell.membership = Microglia_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])),
                                      sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC40_compactness <- mc_compactness(cell.membership = Microglia_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])),
                                      sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC50_compactness <- mc_compactness(cell.membership = Microglia_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])),
                                      sc.reduction = Microglia_Normal_diffusion_comp)

# NeuNRGNI_ASD
NeuNRGNI_ASD_SC10_compactness <- mc_compactness(cell.membership = NeuNRGNI_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC20_compactness <- mc_compactness(cell.membership = NeuNRGNI_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC30_compactness <- mc_compactness(cell.membership = NeuNRGNI_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC40_compactness <- mc_compactness(cell.membership = NeuNRGNI_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC50_compactness <- mc_compactness(cell.membership = NeuNRGNI_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_ASD_diffusion_comp)

# NeuNRGNI_Normal
NeuNRGNI_Normal_SC10_compactness <- mc_compactness(cell.membership = NeuNRGNI_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC20_compactness <- mc_compactness(cell.membership = NeuNRGNI_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC30_compactness <- mc_compactness(cell.membership = NeuNRGNI_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC40_compactness <- mc_compactness(cell.membership = NeuNRGNI_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC50_compactness <- mc_compactness(cell.membership = NeuNRGNI_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])),
                                      sc.reduction = NeuNRGNI_Normal_diffusion_comp)

# INVIP_ASD
INVIP_ASD_SC10_compactness <- mc_compactness(cell.membership = INVIP_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])),
                                      sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC20_compactness <- mc_compactness(cell.membership = INVIP_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])),
                                      sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC30_compactness <- mc_compactness(cell.membership = INVIP_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])),
                                      sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC40_compactness <- mc_compactness(cell.membership = INVIP_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])),
                                      sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC50_compactness <- mc_compactness(cell.membership = INVIP_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])),
                                      sc.reduction = INVIP_ASD_diffusion_comp)

# INVIP_Normal
INVIP_Normal_SC10_compactness <- mc_compactness(cell.membership = INVIP_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])),
                                      sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC20_compactness <- mc_compactness(cell.membership = INVIP_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])),
                                      sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC30_compactness <- mc_compactness(cell.membership = INVIP_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])),
                                      sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC40_compactness <- mc_compactness(cell.membership = INVIP_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])),
                                      sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC50_compactness <- mc_compactness(cell.membership = INVIP_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])),
                                      sc.reduction = INVIP_Normal_diffusion_comp)

# L56CC_ASD
L56CC_ASD_SC10_compactness <- mc_compactness(cell.membership = L56CC_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC20_compactness <- mc_compactness(cell.membership = L56CC_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC30_compactness <- mc_compactness(cell.membership = L56CC_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC40_compactness <- mc_compactness(cell.membership = L56CC_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC50_compactness <- mc_compactness(cell.membership = L56CC_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])),
                                      sc.reduction = L56CC_ASD_diffusion_comp)

# L56CC_Normal
L56CC_Normal_SC10_compactness <- mc_compactness(cell.membership = L56CC_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])),
                                      sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC20_compactness <- mc_compactness(cell.membership = L56CC_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])),
                                      sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC30_compactness <- mc_compactness(cell.membership = L56CC_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])),
                                      sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC40_compactness <- mc_compactness(cell.membership = L56CC_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])),
                                      sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC50_compactness <- mc_compactness(cell.membership = L56CC_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])),
                                      sc.reduction = L56CC_Normal_diffusion_comp)

# INSV2C_ASD
INSV2C_ASD_SC10_compactness <- mc_compactness(cell.membership = INSV2C_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC20_compactness <- mc_compactness(cell.membership = INSV2C_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC30_compactness <- mc_compactness(cell.membership = INSV2C_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC40_compactness <- mc_compactness(cell.membership = INSV2C_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC50_compactness <- mc_compactness(cell.membership = INSV2C_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_ASD_diffusion_comp)

# INSV2C_Normal
INSV2C_Normal_SC10_compactness <- mc_compactness(cell.membership = INSV2C_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC20_compactness <- mc_compactness(cell.membership = INSV2C_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC30_compactness <- mc_compactness(cell.membership = INSV2C_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC40_compactness <- mc_compactness(cell.membership = INSV2C_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC50_compactness <- mc_compactness(cell.membership = INSV2C_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])),
                                      sc.reduction = INSV2C_Normal_diffusion_comp)

# L23_ASD
L23_ASD_SC10_compactness <- mc_compactness(cell.membership = L23_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])),
                                      sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC20_compactness <- mc_compactness(cell.membership = L23_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])),
                                      sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC30_compactness <- mc_compactness(cell.membership = L23_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])),
                                      sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC40_compactness <- mc_compactness(cell.membership = L23_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])),
                                      sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC50_compactness <- mc_compactness(cell.membership = L23_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])),
                                      sc.reduction = L23_ASD_diffusion_comp)

# L23_Normal
L23_Normal_SC10_compactness <- mc_compactness(cell.membership = L23_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])),
                                      sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC20_compactness <- mc_compactness(cell.membership = L23_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])),
                                      sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC30_compactness <- mc_compactness(cell.membership = L23_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])),
                                      sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC40_compactness <- mc_compactness(cell.membership = L23_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])),
                                      sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC50_compactness <- mc_compactness(cell.membership = L23_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])),
                                      sc.reduction = L23_Normal_diffusion_comp)

# INPV_ASD
INPV_ASD_SC10_compactness <- mc_compactness(cell.membership = INPV_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])),
                                      sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC20_compactness <- mc_compactness(cell.membership = INPV_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])),
                                      sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC30_compactness <- mc_compactness(cell.membership = INPV_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])),
                                      sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC40_compactness <- mc_compactness(cell.membership = INPV_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])),
                                      sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC50_compactness <- mc_compactness(cell.membership = INPV_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])),
                                      sc.reduction = INPV_ASD_diffusion_comp)

# INPV_Normal
INPV_Normal_SC10_compactness <- mc_compactness(cell.membership = INPV_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])),
                                      sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC20_compactness <- mc_compactness(cell.membership = INPV_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])),
                                      sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC30_compactness <- mc_compactness(cell.membership = INPV_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])),
                                      sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC40_compactness <- mc_compactness(cell.membership = INPV_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])),
                                      sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC50_compactness <- mc_compactness(cell.membership = INPV_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])),
                                      sc.reduction = INPV_Normal_diffusion_comp)

# L4_ASD
L4_ASD_SC10_compactness <- mc_compactness(cell.membership = L4_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])),
                                      sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC20_compactness <- mc_compactness(cell.membership = L4_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])),
                                      sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC30_compactness <- mc_compactness(cell.membership = L4_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])),
                                      sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC40_compactness <- mc_compactness(cell.membership = L4_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])),
                                      sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC50_compactness <- mc_compactness(cell.membership = L4_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])),
                                      sc.reduction = L4_ASD_diffusion_comp)

# L4_Normal
L4_Normal_SC10_compactness <- mc_compactness(cell.membership = L4_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])),
                                      sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC20_compactness <- mc_compactness(cell.membership = L4_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])),
                                      sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC30_compactness <- mc_compactness(cell.membership = L4_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])),
                                      sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC40_compactness <- mc_compactness(cell.membership = L4_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])),
                                      sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC50_compactness <- mc_compactness(cell.membership = L4_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])),
                                      sc.reduction = L4_Normal_diffusion_comp)

# INSST_ASD
INSST_ASD_SC10_compactness <- mc_compactness(cell.membership = INSST_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC20_compactness <- mc_compactness(cell.membership = INSST_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC30_compactness <- mc_compactness(cell.membership = INSST_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC40_compactness <- mc_compactness(cell.membership = INSST_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC50_compactness <- mc_compactness(cell.membership = INSST_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])),
                                      sc.reduction = INSST_ASD_diffusion_comp)

# INSST_Normal
INSST_Normal_SC10_compactness <- mc_compactness(cell.membership = INSST_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])),
                                      sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC20_compactness <- mc_compactness(cell.membership = INSST_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])),
                                      sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC30_compactness <- mc_compactness(cell.membership = INSST_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])),
                                      sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC40_compactness <- mc_compactness(cell.membership = INSST_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])),
                                      sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC50_compactness <- mc_compactness(cell.membership = INSST_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])),
                                      sc.reduction = INSST_Normal_diffusion_comp)

# Neumat_ASD
Neumat_ASD_SC10_compactness <- mc_compactness(cell.membership = Neumat_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])),
                                      sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC20_compactness <- mc_compactness(cell.membership = Neumat_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])),
                                      sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC30_compactness <- mc_compactness(cell.membership = Neumat_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])),
                                      sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC40_compactness <- mc_compactness(cell.membership = Neumat_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])),
                                      sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC50_compactness <- mc_compactness(cell.membership = Neumat_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])),
                                      sc.reduction = Neumat_ASD_diffusion_comp)

# Neumat_Normal
Neumat_Normal_SC10_compactness <- mc_compactness(cell.membership = Neumat_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])),
                                      sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC20_compactness <- mc_compactness(cell.membership = Neumat_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])),
                                      sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC30_compactness <- mc_compactness(cell.membership = Neumat_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])),
                                      sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC40_compactness <- mc_compactness(cell.membership = Neumat_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])),
                                      sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC50_compactness <- mc_compactness(cell.membership = Neumat_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])),
                                      sc.reduction = Neumat_Normal_diffusion_comp)

# ASTPP_ASD
ASTPP_ASD_SC10_compactness <- mc_compactness(cell.membership = ASTPP_ASD_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC20_compactness <- mc_compactness(cell.membership = ASTPP_ASD_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC30_compactness <- mc_compactness(cell.membership = ASTPP_ASD_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC40_compactness <- mc_compactness(cell.membership = ASTPP_ASD_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC50_compactness <- mc_compactness(cell.membership = ASTPP_ASD_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_ASD_diffusion_comp)

# ASTPP_Normal
ASTPP_Normal_SC10_compactness <- mc_compactness(cell.membership = ASTPP_Normal_SC10_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC20_compactness <- mc_compactness(cell.membership = ASTPP_Normal_SC20_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC30_compactness <- mc_compactness(cell.membership = ASTPP_Normal_SC30_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC40_compactness <- mc_compactness(cell.membership = ASTPP_Normal_SC40_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC50_compactness <- mc_compactness(cell.membership = ASTPP_Normal_SC50_membership,
                                      sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])),
                                      sc.reduction = ASTPP_Normal_diffusion_comp)

# Separation is the distance to the closest metacell. The higher the separation value the better.
# ASD
ASD_SC10_separation <- mc_separation(cell.membership = ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), 
                                    sc.reduction = ASD_diffusion_comp)
ASD_SC20_separation <- mc_separation(cell.membership = ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), 
                                    sc.reduction = ASD_diffusion_comp)
ASD_SC30_separation <- mc_separation(cell.membership = ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), 
                                    sc.reduction = ASD_diffusion_comp)
ASD_SC40_separation <- mc_separation(cell.membership = ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), 
                                    sc.reduction = ASD_diffusion_comp)
ASD_SC50_separation <- mc_separation(cell.membership = ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), 
                                    sc.reduction = ASD_diffusion_comp)

# Normal
Normal_SC10_separation <- mc_separation(cell.membership = Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), 
                                    sc.reduction = Normal_diffusion_comp)
Normal_SC20_separation <- mc_separation(cell.membership = Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), 
                                    sc.reduction = Normal_diffusion_comp)
Normal_SC30_separation <- mc_separation(cell.membership = Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), 
                                    sc.reduction = Normal_diffusion_comp)
Normal_SC40_separation <- mc_separation(cell.membership = Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), 
                                    sc.reduction = Normal_diffusion_comp)
Normal_SC50_separation <- mc_separation(cell.membership = Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), 
                                    sc.reduction = Normal_diffusion_comp)

# ACC_ASD
ACC_ASD_SC10_separation <- mc_separation(cell.membership = ACC_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), 
                                    sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC20_separation <- mc_separation(cell.membership = ACC_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), 
                                    sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC30_separation <- mc_separation(cell.membership = ACC_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), 
                                    sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC40_separation <- mc_separation(cell.membership = ACC_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), 
                                    sc.reduction = ACC_ASD_diffusion_comp)
ACC_ASD_SC50_separation <- mc_separation(cell.membership = ACC_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), 
                                    sc.reduction = ACC_ASD_diffusion_comp)

# ACC_Normal
ACC_Normal_SC10_separation <- mc_separation(cell.membership = ACC_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), 
                                    sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC20_separation <- mc_separation(cell.membership = ACC_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), 
                                    sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC30_separation <- mc_separation(cell.membership = ACC_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), 
                                    sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC40_separation <- mc_separation(cell.membership = ACC_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), 
                                    sc.reduction = ACC_Normal_diffusion_comp)
ACC_Normal_SC50_separation <- mc_separation(cell.membership = ACC_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), 
                                    sc.reduction = ACC_Normal_diffusion_comp)

# PFC_ASD
PFC_ASD_SC10_separation <- mc_separation(cell.membership = PFC_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), 
                                    sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC20_separation <- mc_separation(cell.membership = PFC_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), 
                                    sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC30_separation <- mc_separation(cell.membership = PFC_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), 
                                    sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC40_separation <- mc_separation(cell.membership = PFC_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), 
                                    sc.reduction = PFC_ASD_diffusion_comp)
PFC_ASD_SC50_separation <- mc_separation(cell.membership = PFC_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), 
                                    sc.reduction = PFC_ASD_diffusion_comp)

# PFC_Normal
PFC_Normal_SC10_separation <- mc_separation(cell.membership = PFC_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), 
                                    sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC20_separation <- mc_separation(cell.membership = PFC_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), 
                                    sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC30_separation <- mc_separation(cell.membership = PFC_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), 
                                    sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC40_separation <- mc_separation(cell.membership = PFC_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), 
                                    sc.reduction = PFC_Normal_diffusion_comp)
PFC_Normal_SC50_separation <- mc_separation(cell.membership = PFC_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), 
                                    sc.reduction = PFC_Normal_diffusion_comp)

# Lower18_ASD
Lower18_ASD_SC10_separation <- mc_separation(cell.membership = Lower18_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), 
                                    sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC20_separation <- mc_separation(cell.membership = Lower18_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), 
                                    sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC30_separation <- mc_separation(cell.membership = Lower18_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), 
                                    sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC40_separation <- mc_separation(cell.membership = Lower18_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), 
                                    sc.reduction = Lower18_ASD_diffusion_comp)
Lower18_ASD_SC50_separation <- mc_separation(cell.membership = Lower18_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), 
                                    sc.reduction = Lower18_ASD_diffusion_comp)

# Lower18_Normal
Lower18_Normal_SC10_separation <- mc_separation(cell.membership = Lower18_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), 
                                    sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC20_separation <- mc_separation(cell.membership = Lower18_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), 
                                    sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC30_separation <- mc_separation(cell.membership = Lower18_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), 
                                    sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC40_separation <- mc_separation(cell.membership = Lower18_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), 
                                    sc.reduction = Lower18_Normal_diffusion_comp)
Lower18_Normal_SC50_separation <- mc_separation(cell.membership = Lower18_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), 
                                    sc.reduction = Lower18_Normal_diffusion_comp)

# Larger18_ASD
Larger18_ASD_SC10_separation <- mc_separation(cell.membership = Larger18_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), 
                                    sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC20_separation <- mc_separation(cell.membership = Larger18_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), 
                                    sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC30_separation <- mc_separation(cell.membership = Larger18_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), 
                                    sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC40_separation <- mc_separation(cell.membership = Larger18_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), 
                                    sc.reduction = Larger18_ASD_diffusion_comp)
Larger18_ASD_SC50_separation <- mc_separation(cell.membership = Larger18_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), 
                                    sc.reduction = Larger18_ASD_diffusion_comp)

# Larger18_Normal
Larger18_Normal_SC10_separation <- mc_separation(cell.membership = Larger18_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), 
                                    sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC20_separation <- mc_separation(cell.membership = Larger18_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), 
                                    sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC30_separation <- mc_separation(cell.membership = Larger18_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), 
                                    sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC40_separation <- mc_separation(cell.membership = Larger18_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), 
                                    sc.reduction = Larger18_Normal_diffusion_comp)
Larger18_Normal_SC50_separation <- mc_separation(cell.membership = Larger18_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), 
                                    sc.reduction = Larger18_Normal_diffusion_comp)

# Male_ASD
Male_ASD_SC10_separation <- mc_separation(cell.membership = Male_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), 
                                    sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC20_separation <- mc_separation(cell.membership = Male_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), 
                                    sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC30_separation <- mc_separation(cell.membership = Male_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), 
                                    sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC40_separation <- mc_separation(cell.membership = Male_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), 
                                    sc.reduction = Male_ASD_diffusion_comp)
Male_ASD_SC50_separation <- mc_separation(cell.membership = Male_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), 
                                    sc.reduction = Male_ASD_diffusion_comp)

# Male_Normal
Male_Normal_SC10_separation <- mc_separation(cell.membership = Male_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), 
                                    sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC20_separation <- mc_separation(cell.membership = Male_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), 
                                    sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC30_separation <- mc_separation(cell.membership = Male_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), 
                                    sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC40_separation <- mc_separation(cell.membership = Male_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), 
                                    sc.reduction = Male_Normal_diffusion_comp)
Male_Normal_SC50_separation <- mc_separation(cell.membership = Male_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), 
                                    sc.reduction = Male_Normal_diffusion_comp)

# Female_ASD
Female_ASD_SC10_separation <- mc_separation(cell.membership = Female_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), 
                                    sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC20_separation <- mc_separation(cell.membership = Female_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), 
                                    sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC30_separation <- mc_separation(cell.membership = Female_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), 
                                    sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC40_separation <- mc_separation(cell.membership = Female_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), 
                                    sc.reduction = Female_ASD_diffusion_comp)
Female_ASD_SC50_separation <- mc_separation(cell.membership = Female_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), 
                                    sc.reduction = Female_ASD_diffusion_comp)

# Female_Normal
Female_Normal_SC10_separation <- mc_separation(cell.membership = Female_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), 
                                    sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC20_separation <- mc_separation(cell.membership = Female_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), 
                                    sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC30_separation <- mc_separation(cell.membership = Female_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), 
                                    sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC40_separation <- mc_separation(cell.membership = Female_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), 
                                    sc.reduction = Female_Normal_diffusion_comp)
Female_Normal_SC50_separation <- mc_separation(cell.membership = Female_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), 
                                    sc.reduction = Female_Normal_diffusion_comp)

# NeuNRGNII_ASD
NeuNRGNII_ASD_SC10_separation <- mc_separation(cell.membership = NeuNRGNII_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC20_separation <- mc_separation(cell.membership = NeuNRGNII_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC30_separation <- mc_separation(cell.membership = NeuNRGNII_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC40_separation <- mc_separation(cell.membership = NeuNRGNII_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_ASD_diffusion_comp)
NeuNRGNII_ASD_SC50_separation <- mc_separation(cell.membership = NeuNRGNII_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_ASD_diffusion_comp)

# NeuNRGNII_Normal
NeuNRGNII_Normal_SC10_separation <- mc_separation(cell.membership = NeuNRGNII_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC20_separation <- mc_separation(cell.membership = NeuNRGNII_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC30_separation <- mc_separation(cell.membership = NeuNRGNII_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC40_separation <- mc_separation(cell.membership = NeuNRGNII_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_Normal_diffusion_comp)
NeuNRGNII_Normal_SC50_separation <- mc_separation(cell.membership = NeuNRGNII_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNII_Normal_diffusion_comp)

# L56_ASD
L56_ASD_SC10_separation <- mc_separation(cell.membership = L56_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC20_separation <- mc_separation(cell.membership = L56_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC30_separation <- mc_separation(cell.membership = L56_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC40_separation <- mc_separation(cell.membership = L56_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56_ASD_diffusion_comp)
L56_ASD_SC50_separation <- mc_separation(cell.membership = L56_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56_ASD_diffusion_comp)

# L56_Normal
L56_Normal_SC10_separation <- mc_separation(cell.membership = L56_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), 
                                    sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC20_separation <- mc_separation(cell.membership = L56_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), 
                                    sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC30_separation <- mc_separation(cell.membership = L56_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), 
                                    sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC40_separation <- mc_separation(cell.membership = L56_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), 
                                    sc.reduction = L56_Normal_diffusion_comp)
L56_Normal_SC50_separation <- mc_separation(cell.membership = L56_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), 
                                    sc.reduction = L56_Normal_diffusion_comp)

# Oligodendrocytes_ASD
Oligodendrocytes_ASD_SC10_separation <- mc_separation(cell.membership = Oligodendrocytes_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC20_separation <- mc_separation(cell.membership = Oligodendrocytes_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC30_separation <- mc_separation(cell.membership = Oligodendrocytes_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC40_separation <- mc_separation(cell.membership = Oligodendrocytes_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_ASD_diffusion_comp)
Oligodendrocytes_ASD_SC50_separation <- mc_separation(cell.membership = Oligodendrocytes_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_ASD_diffusion_comp)

# Oligodendrocytes_Normal
Oligodendrocytes_Normal_SC10_separation <- mc_separation(cell.membership = Oligodendrocytes_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC20_separation <- mc_separation(cell.membership = Oligodendrocytes_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC30_separation <- mc_separation(cell.membership = Oligodendrocytes_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC40_separation <- mc_separation(cell.membership = Oligodendrocytes_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_Normal_diffusion_comp)
Oligodendrocytes_Normal_SC50_separation <- mc_separation(cell.membership = Oligodendrocytes_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), 
                                    sc.reduction = Oligodendrocytes_Normal_diffusion_comp)

# OPC_ASD
OPC_ASD_SC10_separation <- mc_separation(cell.membership = OPC_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC20_separation <- mc_separation(cell.membership = OPC_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC30_separation <- mc_separation(cell.membership = OPC_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC40_separation <- mc_separation(cell.membership = OPC_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = OPC_ASD_diffusion_comp)
OPC_ASD_SC50_separation <- mc_separation(cell.membership = OPC_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = OPC_ASD_diffusion_comp)

# OPC_Normal
OPC_Normal_SC10_separation <- mc_separation(cell.membership = OPC_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), 
                                    sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC20_separation <- mc_separation(cell.membership = OPC_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), 
                                    sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC30_separation <- mc_separation(cell.membership = OPC_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), 
                                    sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC40_separation <- mc_separation(cell.membership = OPC_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), 
                                    sc.reduction = OPC_Normal_diffusion_comp)
OPC_Normal_SC50_separation <- mc_separation(cell.membership = OPC_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), 
                                    sc.reduction = OPC_Normal_diffusion_comp)

# ASTFB_ASD
ASTFB_ASD_SC10_separation <- mc_separation(cell.membership = ASTFB_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC20_separation <- mc_separation(cell.membership = ASTFB_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC30_separation <- mc_separation(cell.membership = ASTFB_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC40_separation <- mc_separation(cell.membership = ASTFB_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_ASD_diffusion_comp)
ASTFB_ASD_SC50_separation <- mc_separation(cell.membership = ASTFB_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_ASD_diffusion_comp)

# ASTFB_Normal
ASTFB_Normal_SC10_separation <- mc_separation(cell.membership = ASTFB_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC20_separation <- mc_separation(cell.membership = ASTFB_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC30_separation <- mc_separation(cell.membership = ASTFB_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC40_separation <- mc_separation(cell.membership = ASTFB_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_Normal_diffusion_comp)
ASTFB_Normal_SC50_separation <- mc_separation(cell.membership = ASTFB_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTFB_Normal_diffusion_comp)

# Endothelial_ASD
Endothelial_ASD_SC10_separation <- mc_separation(cell.membership = Endothelial_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC20_separation <- mc_separation(cell.membership = Endothelial_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC30_separation <- mc_separation(cell.membership = Endothelial_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC40_separation <- mc_separation(cell.membership = Endothelial_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_ASD_diffusion_comp)
Endothelial_ASD_SC50_separation <- mc_separation(cell.membership = Endothelial_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_ASD_diffusion_comp)

# Endothelial_Normal
Endothelial_Normal_SC10_separation <- mc_separation(cell.membership = Endothelial_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC20_separation <- mc_separation(cell.membership = Endothelial_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC30_separation <- mc_separation(cell.membership = Endothelial_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC40_separation <- mc_separation(cell.membership = Endothelial_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_Normal_diffusion_comp)
Endothelial_Normal_SC50_separation <- mc_separation(cell.membership = Endothelial_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), 
                                    sc.reduction = Endothelial_Normal_diffusion_comp)

# Microglia_ASD
Microglia_ASD_SC10_separation <- mc_separation(cell.membership = Microglia_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC20_separation <- mc_separation(cell.membership = Microglia_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC30_separation <- mc_separation(cell.membership = Microglia_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC40_separation <- mc_separation(cell.membership = Microglia_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_ASD_diffusion_comp)
Microglia_ASD_SC50_separation <- mc_separation(cell.membership = Microglia_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_ASD_diffusion_comp)

# Microglia_Normal
Microglia_Normal_SC10_separation <- mc_separation(cell.membership = Microglia_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC20_separation <- mc_separation(cell.membership = Microglia_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC30_separation <- mc_separation(cell.membership = Microglia_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC40_separation <- mc_separation(cell.membership = Microglia_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_Normal_diffusion_comp)
Microglia_Normal_SC50_separation <- mc_separation(cell.membership = Microglia_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), 
                                    sc.reduction = Microglia_Normal_diffusion_comp)

# NeuNRGNI_ASD
NeuNRGNI_ASD_SC10_separation <- mc_separation(cell.membership = NeuNRGNI_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC20_separation <- mc_separation(cell.membership = NeuNRGNI_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC30_separation <- mc_separation(cell.membership = NeuNRGNI_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC40_separation <- mc_separation(cell.membership = NeuNRGNI_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_ASD_diffusion_comp)
NeuNRGNI_ASD_SC50_separation <- mc_separation(cell.membership = NeuNRGNI_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_ASD_diffusion_comp)

# NeuNRGNI_Normal
NeuNRGNI_Normal_SC10_separation <- mc_separation(cell.membership = NeuNRGNI_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC20_separation <- mc_separation(cell.membership = NeuNRGNI_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC30_separation <- mc_separation(cell.membership = NeuNRGNI_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC40_separation <- mc_separation(cell.membership = NeuNRGNI_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_Normal_diffusion_comp)
NeuNRGNI_Normal_SC50_separation <- mc_separation(cell.membership = NeuNRGNI_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), 
                                    sc.reduction = NeuNRGNI_Normal_diffusion_comp)

# INVIP_ASD
INVIP_ASD_SC10_separation <- mc_separation(cell.membership = INVIP_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC20_separation <- mc_separation(cell.membership = INVIP_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC30_separation <- mc_separation(cell.membership = INVIP_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC40_separation <- mc_separation(cell.membership = INVIP_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_ASD_diffusion_comp)
INVIP_ASD_SC50_separation <- mc_separation(cell.membership = INVIP_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_ASD_diffusion_comp)

# INVIP_Normal
INVIP_Normal_SC10_separation <- mc_separation(cell.membership = INVIP_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC20_separation <- mc_separation(cell.membership = INVIP_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC30_separation <- mc_separation(cell.membership = INVIP_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC40_separation <- mc_separation(cell.membership = INVIP_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_Normal_diffusion_comp)
INVIP_Normal_SC50_separation <- mc_separation(cell.membership = INVIP_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), 
                                    sc.reduction = INVIP_Normal_diffusion_comp)

# L56CC_ASD
L56CC_ASD_SC10_separation <- mc_separation(cell.membership = L56CC_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC20_separation <- mc_separation(cell.membership = L56CC_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC30_separation <- mc_separation(cell.membership = L56CC_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC40_separation <- mc_separation(cell.membership = L56CC_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_ASD_diffusion_comp)
L56CC_ASD_SC50_separation <- mc_separation(cell.membership = L56CC_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_ASD_diffusion_comp)

# L56CC_Normal
L56CC_Normal_SC10_separation <- mc_separation(cell.membership = L56CC_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC20_separation <- mc_separation(cell.membership = L56CC_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC30_separation <- mc_separation(cell.membership = L56CC_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC40_separation <- mc_separation(cell.membership = L56CC_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_Normal_diffusion_comp)
L56CC_Normal_SC50_separation <- mc_separation(cell.membership = L56CC_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), 
                                    sc.reduction = L56CC_Normal_diffusion_comp)

# INSV2C_ASD
INSV2C_ASD_SC10_separation <- mc_separation(cell.membership = INSV2C_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC20_separation <- mc_separation(cell.membership = INSV2C_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC30_separation <- mc_separation(cell.membership = INSV2C_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC40_separation <- mc_separation(cell.membership = INSV2C_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_ASD_diffusion_comp)
INSV2C_ASD_SC50_separation <- mc_separation(cell.membership = INSV2C_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_ASD_diffusion_comp)

# INSV2C_Normal
INSV2C_Normal_SC10_separation <- mc_separation(cell.membership = INSV2C_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC20_separation <- mc_separation(cell.membership = INSV2C_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC30_separation <- mc_separation(cell.membership = INSV2C_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC40_separation <- mc_separation(cell.membership = INSV2C_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_Normal_diffusion_comp)
INSV2C_Normal_SC50_separation <- mc_separation(cell.membership = INSV2C_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), 
                                    sc.reduction = INSV2C_Normal_diffusion_comp)

# L23_ASD
L23_ASD_SC10_separation <- mc_separation(cell.membership = L23_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC20_separation <- mc_separation(cell.membership = L23_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC30_separation <- mc_separation(cell.membership = L23_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC40_separation <- mc_separation(cell.membership = L23_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L23_ASD_diffusion_comp)
L23_ASD_SC50_separation <- mc_separation(cell.membership = L23_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L23_ASD_diffusion_comp)

# L23_Normal
L23_Normal_SC10_separation <- mc_separation(cell.membership = L23_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), 
                                    sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC20_separation <- mc_separation(cell.membership = L23_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), 
                                    sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC30_separation <- mc_separation(cell.membership = L23_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), 
                                    sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC40_separation <- mc_separation(cell.membership = L23_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), 
                                    sc.reduction = L23_Normal_diffusion_comp)
L23_Normal_SC50_separation <- mc_separation(cell.membership = L23_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), 
                                    sc.reduction = L23_Normal_diffusion_comp)

# INPV_ASD
INPV_ASD_SC10_separation <- mc_separation(cell.membership = INPV_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC20_separation <- mc_separation(cell.membership = INPV_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC30_separation <- mc_separation(cell.membership = INPV_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC40_separation <- mc_separation(cell.membership = INPV_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INPV_ASD_diffusion_comp)
INPV_ASD_SC50_separation <- mc_separation(cell.membership = INPV_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INPV_ASD_diffusion_comp)

# INPV_Normal
INPV_Normal_SC10_separation <- mc_separation(cell.membership = INPV_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), 
                                    sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC20_separation <- mc_separation(cell.membership = INPV_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), 
                                    sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC30_separation <- mc_separation(cell.membership = INPV_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), 
                                    sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC40_separation <- mc_separation(cell.membership = INPV_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), 
                                    sc.reduction = INPV_Normal_diffusion_comp)
INPV_Normal_SC50_separation <- mc_separation(cell.membership = INPV_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), 
                                    sc.reduction = INPV_Normal_diffusion_comp)

# L4_ASD
L4_ASD_SC10_separation <- mc_separation(cell.membership = L4_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC20_separation <- mc_separation(cell.membership = L4_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC30_separation <- mc_separation(cell.membership = L4_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC40_separation <- mc_separation(cell.membership = L4_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L4_ASD_diffusion_comp)
L4_ASD_SC50_separation <- mc_separation(cell.membership = L4_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), 
                                    sc.reduction = L4_ASD_diffusion_comp)

# L4_Normal
L4_Normal_SC10_separation <- mc_separation(cell.membership = L4_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), 
                                    sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC20_separation <- mc_separation(cell.membership = L4_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), 
                                    sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC30_separation <- mc_separation(cell.membership = L4_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), 
                                    sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC40_separation <- mc_separation(cell.membership = L4_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), 
                                    sc.reduction = L4_Normal_diffusion_comp)
L4_Normal_SC50_separation <- mc_separation(cell.membership = L4_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), 
                                    sc.reduction = L4_Normal_diffusion_comp)

# INSST_ASD
INSST_ASD_SC10_separation <- mc_separation(cell.membership = INSST_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC20_separation <- mc_separation(cell.membership = INSST_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC30_separation <- mc_separation(cell.membership = INSST_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC40_separation <- mc_separation(cell.membership = INSST_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSST_ASD_diffusion_comp)
INSST_ASD_SC50_separation <- mc_separation(cell.membership = INSST_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), 
                                    sc.reduction = INSST_ASD_diffusion_comp)

# INSST_Normal
INSST_Normal_SC10_separation <- mc_separation(cell.membership = INSST_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), 
                                    sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC20_separation <- mc_separation(cell.membership = INSST_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), 
                                    sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC30_separation <- mc_separation(cell.membership = INSST_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), 
                                    sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC40_separation <- mc_separation(cell.membership = INSST_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), 
                                    sc.reduction = INSST_Normal_diffusion_comp)
INSST_Normal_SC50_separation <- mc_separation(cell.membership = INSST_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), 
                                    sc.reduction = INSST_Normal_diffusion_comp)

# Neumat_ASD
Neumat_ASD_SC10_separation <- mc_separation(cell.membership = Neumat_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC20_separation <- mc_separation(cell.membership = Neumat_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC30_separation <- mc_separation(cell.membership = Neumat_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC40_separation <- mc_separation(cell.membership = Neumat_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_ASD_diffusion_comp)
Neumat_ASD_SC50_separation <- mc_separation(cell.membership = Neumat_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_ASD_diffusion_comp)

# Neumat_Normal
Neumat_Normal_SC10_separation <- mc_separation(cell.membership = Neumat_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC20_separation <- mc_separation(cell.membership = Neumat_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC30_separation <- mc_separation(cell.membership = Neumat_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC40_separation <- mc_separation(cell.membership = Neumat_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_Normal_diffusion_comp)
Neumat_Normal_SC50_separation <- mc_separation(cell.membership = Neumat_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), 
                                    sc.reduction = Neumat_Normal_diffusion_comp)

# ASTPP_ASD
ASTPP_ASD_SC10_separation <- mc_separation(cell.membership = ASTPP_ASD_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC20_separation <- mc_separation(cell.membership = ASTPP_ASD_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC30_separation <- mc_separation(cell.membership = ASTPP_ASD_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC40_separation <- mc_separation(cell.membership = ASTPP_ASD_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_ASD_diffusion_comp)
ASTPP_ASD_SC50_separation <- mc_separation(cell.membership = ASTPP_ASD_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_ASD_diffusion_comp)

# ASTPP_Normal
ASTPP_Normal_SC10_separation <- mc_separation(cell.membership = ASTPP_Normal_SC10_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC20_separation <- mc_separation(cell.membership = ASTPP_Normal_SC20_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC30_separation <- mc_separation(cell.membership = ASTPP_Normal_SC30_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC40_separation <- mc_separation(cell.membership = ASTPP_Normal_SC40_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_Normal_diffusion_comp)
ASTPP_Normal_SC50_separation <- mc_separation(cell.membership = ASTPP_Normal_SC50_membership, 
                                    sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), 
                                    sc.reduction = ASTPP_Normal_diffusion_comp)

# INV is the mean-normalized variance of gene expression within the metacell. The lower the INV value the better. Note that it is the only metric that is latent-space independent.
# ASD
ASD_SC10_INV <- mc_INV(cell.membership = ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), group.label = "membership")
ASD_SC20_INV <- mc_INV(cell.membership = ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), group.label = "membership")
ASD_SC30_INV <- mc_INV(cell.membership = ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), group.label = "membership")
ASD_SC40_INV <- mc_INV(cell.membership = ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), group.label = "membership")
ASD_SC50_INV <- mc_INV(cell.membership = ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), group.label = "membership")

# Normal
Normal_SC10_INV <- mc_INV(cell.membership = Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), group.label = "membership")
Normal_SC20_INV <- mc_INV(cell.membership = Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), group.label = "membership")
Normal_SC30_INV <- mc_INV(cell.membership = Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), group.label = "membership")
Normal_SC40_INV <- mc_INV(cell.membership = Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), group.label = "membership")
Normal_SC50_INV <- mc_INV(cell.membership = Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_mRNAs_data, Control_lncRNAs_data)), group.label = "membership")

# ACC_ASD
ACC_ASD_SC10_INV <- mc_INV(cell.membership = ACC_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), group.label = "membership")
ACC_ASD_SC20_INV <- mc_INV(cell.membership = ACC_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), group.label = "membership")
ACC_ASD_SC30_INV <- mc_INV(cell.membership = ACC_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), group.label = "membership")
ACC_ASD_SC40_INV <- mc_INV(cell.membership = ACC_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), group.label = "membership")
ACC_ASD_SC50_INV <- mc_INV(cell.membership = ACC_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ACC_ASD_mRNAs_data, ACC_ASD_lncRNAs_data)), group.label = "membership")

# ACC_Normal
ACC_Normal_SC10_INV <- mc_INV(cell.membership = ACC_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), group.label = "membership")
ACC_Normal_SC20_INV <- mc_INV(cell.membership = ACC_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), group.label = "membership")
ACC_Normal_SC30_INV <- mc_INV(cell.membership = ACC_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), group.label = "membership")
ACC_Normal_SC40_INV <- mc_INV(cell.membership = ACC_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), group.label = "membership")
ACC_Normal_SC50_INV <- mc_INV(cell.membership = ACC_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(ACC_control_mRNAs_data, ACC_control_lncRNAs_data)), group.label = "membership")

# PFC_ASD
PFC_ASD_SC10_INV <- mc_INV(cell.membership = PFC_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), group.label = "membership")
PFC_ASD_SC20_INV <- mc_INV(cell.membership = PFC_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), group.label = "membership")
PFC_ASD_SC30_INV <- mc_INV(cell.membership = PFC_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), group.label = "membership")
PFC_ASD_SC40_INV <- mc_INV(cell.membership = PFC_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), group.label = "membership")
PFC_ASD_SC50_INV <- mc_INV(cell.membership = PFC_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(PFC_ASD_mRNAs_data, PFC_ASD_lncRNAs_data)), group.label = "membership")

# PFC_Normal
PFC_Normal_SC10_INV <- mc_INV(cell.membership = PFC_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), group.label = "membership")
PFC_Normal_SC20_INV <- mc_INV(cell.membership = PFC_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), group.label = "membership")
PFC_Normal_SC30_INV <- mc_INV(cell.membership = PFC_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), group.label = "membership")
PFC_Normal_SC40_INV <- mc_INV(cell.membership = PFC_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), group.label = "membership")
PFC_Normal_SC50_INV <- mc_INV(cell.membership = PFC_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(PFC_control_mRNAs_data, PFC_control_lncRNAs_data)), group.label = "membership")

# Lower18_ASD
Lower18_ASD_SC10_INV <- mc_INV(cell.membership = Lower18_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), group.label = "membership")
Lower18_ASD_SC20_INV <- mc_INV(cell.membership = Lower18_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), group.label = "membership")
Lower18_ASD_SC30_INV <- mc_INV(cell.membership = Lower18_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), group.label = "membership")
Lower18_ASD_SC40_INV <- mc_INV(cell.membership = Lower18_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), group.label = "membership")
Lower18_ASD_SC50_INV <- mc_INV(cell.membership = Lower18_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(Lower18_ASD_mRNAs_data, Lower18_ASD_lncRNAs_data)), group.label = "membership")

# Lower18_Normal
Lower18_Normal_SC10_INV <- mc_INV(cell.membership = Lower18_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), group.label = "membership")
Lower18_Normal_SC20_INV <- mc_INV(cell.membership = Lower18_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), group.label = "membership")
Lower18_Normal_SC30_INV <- mc_INV(cell.membership = Lower18_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), group.label = "membership")
Lower18_Normal_SC40_INV <- mc_INV(cell.membership = Lower18_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), group.label = "membership")
Lower18_Normal_SC50_INV <- mc_INV(cell.membership = Lower18_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Lower18_control_mRNAs_data, Lower18_control_lncRNAs_data)), group.label = "membership")

# Larger18_ASD
Larger18_ASD_SC10_INV <- mc_INV(cell.membership = Larger18_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), group.label = "membership")
Larger18_ASD_SC20_INV <- mc_INV(cell.membership = Larger18_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), group.label = "membership")
Larger18_ASD_SC30_INV <- mc_INV(cell.membership = Larger18_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), group.label = "membership")
Larger18_ASD_SC40_INV <- mc_INV(cell.membership = Larger18_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), group.label = "membership")
Larger18_ASD_SC50_INV <- mc_INV(cell.membership = Larger18_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(Larger18_ASD_mRNAs_data, Larger18_ASD_lncRNAs_data)), group.label = "membership")

# Larger18_Normal
Larger18_Normal_SC10_INV <- mc_INV(cell.membership = Larger18_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), group.label = "membership")
Larger18_Normal_SC20_INV <- mc_INV(cell.membership = Larger18_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), group.label = "membership")
Larger18_Normal_SC30_INV <- mc_INV(cell.membership = Larger18_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), group.label = "membership")
Larger18_Normal_SC40_INV <- mc_INV(cell.membership = Larger18_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), group.label = "membership")
Larger18_Normal_SC50_INV <- mc_INV(cell.membership = Larger18_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Larger18_control_mRNAs_data, Larger18_control_lncRNAs_data)), group.label = "membership")

# Male_ASD
Male_ASD_SC10_INV <- mc_INV(cell.membership = Male_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), group.label = "membership")
Male_ASD_SC20_INV <- mc_INV(cell.membership = Male_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), group.label = "membership")
Male_ASD_SC30_INV <- mc_INV(cell.membership = Male_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), group.label = "membership")
Male_ASD_SC40_INV <- mc_INV(cell.membership = Male_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), group.label = "membership")
Male_ASD_SC50_INV <- mc_INV(cell.membership = Male_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(Male_ASD_mRNAs_data, Male_ASD_lncRNAs_data)), group.label = "membership")

# Male_Normal
Male_Normal_SC10_INV <- mc_INV(cell.membership = Male_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), group.label = "membership")
Male_Normal_SC20_INV <- mc_INV(cell.membership = Male_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), group.label = "membership")
Male_Normal_SC30_INV <- mc_INV(cell.membership = Male_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), group.label = "membership")
Male_Normal_SC40_INV <- mc_INV(cell.membership = Male_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), group.label = "membership")
Male_Normal_SC50_INV <- mc_INV(cell.membership = Male_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Male_control_mRNAs_data, Male_control_lncRNAs_data)), group.label = "membership")

# Female_ASD
Female_ASD_SC10_INV <- mc_INV(cell.membership = Female_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), group.label = "membership")
Female_ASD_SC20_INV <- mc_INV(cell.membership = Female_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), group.label = "membership")
Female_ASD_SC30_INV <- mc_INV(cell.membership = Female_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), group.label = "membership")
Female_ASD_SC40_INV <- mc_INV(cell.membership = Female_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), group.label = "membership")
Female_ASD_SC50_INV <- mc_INV(cell.membership = Female_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(Female_ASD_mRNAs_data, Female_ASD_lncRNAs_data)), group.label = "membership")

# Female_Normal
Female_Normal_SC10_INV <- mc_INV(cell.membership = Female_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), group.label = "membership")
Female_Normal_SC20_INV <- mc_INV(cell.membership = Female_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), group.label = "membership")
Female_Normal_SC30_INV <- mc_INV(cell.membership = Female_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), group.label = "membership")
Female_Normal_SC40_INV <- mc_INV(cell.membership = Female_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), group.label = "membership")
Female_Normal_SC50_INV <- mc_INV(cell.membership = Female_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Female_control_mRNAs_data, Female_control_lncRNAs_data)), group.label = "membership")

# NeuNRGNII_ASD
NeuNRGNII_ASD_SC10_INV <- mc_INV(cell.membership = NeuNRGNII_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_ASD_SC20_INV <- mc_INV(cell.membership = NeuNRGNII_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_ASD_SC30_INV <- mc_INV(cell.membership = NeuNRGNII_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_ASD_SC40_INV <- mc_INV(cell.membership = NeuNRGNII_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_ASD_SC50_INV <- mc_INV(cell.membership = NeuNRGNII_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNII_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNII_ASD_lncRNAs_data"]])), group.label = "membership")

# NeuNRGNII_Normal
NeuNRGNII_Normal_SC10_INV <- mc_INV(cell.membership = NeuNRGNII_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_Normal_SC20_INV <- mc_INV(cell.membership = NeuNRGNII_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_Normal_SC30_INV <- mc_INV(cell.membership = NeuNRGNII_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_Normal_SC40_INV <- mc_INV(cell.membership = NeuNRGNII_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNII_Normal_SC50_INV <- mc_INV(cell.membership = NeuNRGNII_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNII_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNII_control_lncRNAs_data"]])), group.label = "membership")

# L56_ASD
L56_ASD_SC10_INV <- mc_INV(cell.membership = L56_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), group.label = "membership")
L56_ASD_SC20_INV <- mc_INV(cell.membership = L56_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), group.label = "membership")
L56_ASD_SC30_INV <- mc_INV(cell.membership = L56_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), group.label = "membership")
L56_ASD_SC40_INV <- mc_INV(cell.membership = L56_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), group.label = "membership")
L56_ASD_SC50_INV <- mc_INV(cell.membership = L56_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56_ASD_lncRNAs_data"]])), group.label = "membership")

# L56_Normal
L56_Normal_SC10_INV <- mc_INV(cell.membership = L56_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), group.label = "membership")
L56_Normal_SC20_INV <- mc_INV(cell.membership = L56_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), group.label = "membership")
L56_Normal_SC30_INV <- mc_INV(cell.membership = L56_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), group.label = "membership")
L56_Normal_SC40_INV <- mc_INV(cell.membership = L56_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), group.label = "membership")
L56_Normal_SC50_INV <- mc_INV(cell.membership = L56_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56_control_lncRNAs_data"]])), group.label = "membership")

# Oligodendrocytes_ASD
Oligodendrocytes_ASD_SC10_INV <- mc_INV(cell.membership = Oligodendrocytes_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_ASD_SC20_INV <- mc_INV(cell.membership = Oligodendrocytes_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_ASD_SC30_INV <- mc_INV(cell.membership = Oligodendrocytes_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_ASD_SC40_INV <- mc_INV(cell.membership = Oligodendrocytes_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_ASD_SC50_INV <- mc_INV(cell.membership = Oligodendrocytes_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Oligodendrocytes_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Oligodendrocytes_ASD_lncRNAs_data"]])), group.label = "membership")

# Oligodendrocytes_Normal
Oligodendrocytes_Normal_SC10_INV <- mc_INV(cell.membership = Oligodendrocytes_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_Normal_SC20_INV <- mc_INV(cell.membership = Oligodendrocytes_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_Normal_SC30_INV <- mc_INV(cell.membership = Oligodendrocytes_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_Normal_SC40_INV <- mc_INV(cell.membership = Oligodendrocytes_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), group.label = "membership")
Oligodendrocytes_Normal_SC50_INV <- mc_INV(cell.membership = Oligodendrocytes_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Oligodendrocytes_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Oligodendrocytes_control_lncRNAs_data"]])), group.label = "membership")

# OPC_ASD
OPC_ASD_SC10_INV <- mc_INV(cell.membership = OPC_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), group.label = "membership")
OPC_ASD_SC20_INV <- mc_INV(cell.membership = OPC_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), group.label = "membership")
OPC_ASD_SC30_INV <- mc_INV(cell.membership = OPC_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), group.label = "membership")
OPC_ASD_SC40_INV <- mc_INV(cell.membership = OPC_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), group.label = "membership")
OPC_ASD_SC50_INV <- mc_INV(cell.membership = OPC_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["OPC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["OPC_ASD_lncRNAs_data"]])), group.label = "membership")

# OPC_Normal
OPC_Normal_SC10_INV <- mc_INV(cell.membership = OPC_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), group.label = "membership")
OPC_Normal_SC20_INV <- mc_INV(cell.membership = OPC_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), group.label = "membership")
OPC_Normal_SC30_INV <- mc_INV(cell.membership = OPC_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), group.label = "membership")
OPC_Normal_SC40_INV <- mc_INV(cell.membership = OPC_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), group.label = "membership")
OPC_Normal_SC50_INV <- mc_INV(cell.membership = OPC_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["OPC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["OPC_control_lncRNAs_data"]])), group.label = "membership")

# ASTFB_ASD
ASTFB_ASD_SC10_INV <- mc_INV(cell.membership = ASTFB_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), group.label = "membership")
ASTFB_ASD_SC20_INV <- mc_INV(cell.membership = ASTFB_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), group.label = "membership")
ASTFB_ASD_SC30_INV <- mc_INV(cell.membership = ASTFB_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), group.label = "membership")
ASTFB_ASD_SC40_INV <- mc_INV(cell.membership = ASTFB_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), group.label = "membership")
ASTFB_ASD_SC50_INV <- mc_INV(cell.membership = ASTFB_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTFB_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTFB_ASD_lncRNAs_data"]])), group.label = "membership")

# ASTFB_Normal
ASTFB_Normal_SC10_INV <- mc_INV(cell.membership = ASTFB_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), group.label = "membership")
ASTFB_Normal_SC20_INV <- mc_INV(cell.membership = ASTFB_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), group.label = "membership")
ASTFB_Normal_SC30_INV <- mc_INV(cell.membership = ASTFB_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), group.label = "membership")
ASTFB_Normal_SC40_INV <- mc_INV(cell.membership = ASTFB_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), group.label = "membership")
ASTFB_Normal_SC50_INV <- mc_INV(cell.membership = ASTFB_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTFB_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTFB_control_lncRNAs_data"]])), group.label = "membership")

# Endothelial_ASD
Endothelial_ASD_SC10_INV <- mc_INV(cell.membership = Endothelial_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), group.label = "membership")
Endothelial_ASD_SC20_INV <- mc_INV(cell.membership = Endothelial_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), group.label = "membership")
Endothelial_ASD_SC30_INV <- mc_INV(cell.membership = Endothelial_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), group.label = "membership")
Endothelial_ASD_SC40_INV <- mc_INV(cell.membership = Endothelial_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), group.label = "membership")
Endothelial_ASD_SC50_INV <- mc_INV(cell.membership = Endothelial_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Endothelial_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Endothelial_ASD_lncRNAs_data"]])), group.label = "membership")

# Endothelial_Normal
Endothelial_Normal_SC10_INV <- mc_INV(cell.membership = Endothelial_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), group.label = "membership")
Endothelial_Normal_SC20_INV <- mc_INV(cell.membership = Endothelial_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), group.label = "membership")
Endothelial_Normal_SC30_INV <- mc_INV(cell.membership = Endothelial_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), group.label = "membership")
Endothelial_Normal_SC40_INV <- mc_INV(cell.membership = Endothelial_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), group.label = "membership")
Endothelial_Normal_SC50_INV <- mc_INV(cell.membership = Endothelial_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Endothelial_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Endothelial_control_lncRNAs_data"]])), group.label = "membership")

# Microglia_ASD
Microglia_ASD_SC10_INV <- mc_INV(cell.membership = Microglia_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), group.label = "membership")
Microglia_ASD_SC20_INV <- mc_INV(cell.membership = Microglia_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), group.label = "membership")
Microglia_ASD_SC30_INV <- mc_INV(cell.membership = Microglia_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), group.label = "membership")
Microglia_ASD_SC40_INV <- mc_INV(cell.membership = Microglia_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), group.label = "membership")
Microglia_ASD_SC50_INV <- mc_INV(cell.membership = Microglia_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Microglia_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Microglia_ASD_lncRNAs_data"]])), group.label = "membership")

# Microglia_Normal
Microglia_Normal_SC10_INV <- mc_INV(cell.membership = Microglia_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), group.label = "membership")
Microglia_Normal_SC20_INV <- mc_INV(cell.membership = Microglia_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), group.label = "membership")
Microglia_Normal_SC30_INV <- mc_INV(cell.membership = Microglia_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), group.label = "membership")
Microglia_Normal_SC40_INV <- mc_INV(cell.membership = Microglia_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), group.label = "membership")
Microglia_Normal_SC50_INV <- mc_INV(cell.membership = Microglia_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Microglia_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Microglia_control_lncRNAs_data"]])), group.label = "membership")

# NeuNRGNI_ASD
NeuNRGNI_ASD_SC10_INV <- mc_INV(cell.membership = NeuNRGNI_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_ASD_SC20_INV <- mc_INV(cell.membership = NeuNRGNI_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_ASD_SC30_INV <- mc_INV(cell.membership = NeuNRGNI_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_ASD_SC40_INV <- mc_INV(cell.membership = NeuNRGNI_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_ASD_SC50_INV <- mc_INV(cell.membership = NeuNRGNI_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["NeuNRGNI_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["NeuNRGNI_ASD_lncRNAs_data"]])), group.label = "membership")

# NeuNRGNI_Normal
NeuNRGNI_Normal_SC10_INV <- mc_INV(cell.membership = NeuNRGNI_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_Normal_SC20_INV <- mc_INV(cell.membership = NeuNRGNI_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_Normal_SC30_INV <- mc_INV(cell.membership = NeuNRGNI_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_Normal_SC40_INV <- mc_INV(cell.membership = NeuNRGNI_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), group.label = "membership")
NeuNRGNI_Normal_SC50_INV <- mc_INV(cell.membership = NeuNRGNI_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["NeuNRGNI_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["NeuNRGNI_control_lncRNAs_data"]])), group.label = "membership")

# INVIP_ASD
INVIP_ASD_SC10_INV <- mc_INV(cell.membership = INVIP_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), group.label = "membership")
INVIP_ASD_SC20_INV <- mc_INV(cell.membership = INVIP_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), group.label = "membership")
INVIP_ASD_SC30_INV <- mc_INV(cell.membership = INVIP_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), group.label = "membership")
INVIP_ASD_SC40_INV <- mc_INV(cell.membership = INVIP_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), group.label = "membership")
INVIP_ASD_SC50_INV <- mc_INV(cell.membership = INVIP_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INVIP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INVIP_ASD_lncRNAs_data"]])), group.label = "membership")

# INVIP_Normal
INVIP_Normal_SC10_INV <- mc_INV(cell.membership = INVIP_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), group.label = "membership")
INVIP_Normal_SC20_INV <- mc_INV(cell.membership = INVIP_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), group.label = "membership")
INVIP_Normal_SC30_INV <- mc_INV(cell.membership = INVIP_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), group.label = "membership")
INVIP_Normal_SC40_INV <- mc_INV(cell.membership = INVIP_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), group.label = "membership")
INVIP_Normal_SC50_INV <- mc_INV(cell.membership = INVIP_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INVIP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INVIP_control_lncRNAs_data"]])), group.label = "membership")

# L56CC_ASD
L56CC_ASD_SC10_INV <- mc_INV(cell.membership = L56CC_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), group.label = "membership")
L56CC_ASD_SC20_INV <- mc_INV(cell.membership = L56CC_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), group.label = "membership")
L56CC_ASD_SC30_INV <- mc_INV(cell.membership = L56CC_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), group.label = "membership")
L56CC_ASD_SC40_INV <- mc_INV(cell.membership = L56CC_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), group.label = "membership")
L56CC_ASD_SC50_INV <- mc_INV(cell.membership = L56CC_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L56CC_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L56CC_ASD_lncRNAs_data"]])), group.label = "membership")

# L56CC_Normal
L56CC_Normal_SC10_INV <- mc_INV(cell.membership = L56CC_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), group.label = "membership")
L56CC_Normal_SC20_INV <- mc_INV(cell.membership = L56CC_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), group.label = "membership")
L56CC_Normal_SC30_INV <- mc_INV(cell.membership = L56CC_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), group.label = "membership")
L56CC_Normal_SC40_INV <- mc_INV(cell.membership = L56CC_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), group.label = "membership")
L56CC_Normal_SC50_INV <- mc_INV(cell.membership = L56CC_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L56CC_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L56CC_control_lncRNAs_data"]])), group.label = "membership")

# INSV2C_ASD
INSV2C_ASD_SC10_INV <- mc_INV(cell.membership = INSV2C_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), group.label = "membership")
INSV2C_ASD_SC20_INV <- mc_INV(cell.membership = INSV2C_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), group.label = "membership")
INSV2C_ASD_SC30_INV <- mc_INV(cell.membership = INSV2C_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), group.label = "membership")
INSV2C_ASD_SC40_INV <- mc_INV(cell.membership = INSV2C_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), group.label = "membership")
INSV2C_ASD_SC50_INV <- mc_INV(cell.membership = INSV2C_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSV2C_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSV2C_ASD_lncRNAs_data"]])), group.label = "membership")

# INSV2C_Normal
INSV2C_Normal_SC10_INV <- mc_INV(cell.membership = INSV2C_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), group.label = "membership")
INSV2C_Normal_SC20_INV <- mc_INV(cell.membership = INSV2C_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), group.label = "membership")
INSV2C_Normal_SC30_INV <- mc_INV(cell.membership = INSV2C_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), group.label = "membership")
INSV2C_Normal_SC40_INV <- mc_INV(cell.membership = INSV2C_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), group.label = "membership")
INSV2C_Normal_SC50_INV <- mc_INV(cell.membership = INSV2C_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSV2C_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSV2C_control_lncRNAs_data"]])), group.label = "membership")

# L23_ASD
L23_ASD_SC10_INV <- mc_INV(cell.membership = L23_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), group.label = "membership")
L23_ASD_SC20_INV <- mc_INV(cell.membership = L23_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), group.label = "membership")
L23_ASD_SC30_INV <- mc_INV(cell.membership = L23_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), group.label = "membership")
L23_ASD_SC40_INV <- mc_INV(cell.membership = L23_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), group.label = "membership")
L23_ASD_SC50_INV <- mc_INV(cell.membership = L23_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L23_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L23_ASD_lncRNAs_data"]])), group.label = "membership")

# L23_Normal
L23_Normal_SC10_INV <- mc_INV(cell.membership = L23_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), group.label = "membership")
L23_Normal_SC20_INV <- mc_INV(cell.membership = L23_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), group.label = "membership")
L23_Normal_SC30_INV <- mc_INV(cell.membership = L23_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), group.label = "membership")
L23_Normal_SC40_INV <- mc_INV(cell.membership = L23_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), group.label = "membership")
L23_Normal_SC50_INV <- mc_INV(cell.membership = L23_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L23_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L23_control_lncRNAs_data"]])), group.label = "membership")

# INPV_ASD
INPV_ASD_SC10_INV <- mc_INV(cell.membership = INPV_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), group.label = "membership")
INPV_ASD_SC20_INV <- mc_INV(cell.membership = INPV_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), group.label = "membership")
INPV_ASD_SC30_INV <- mc_INV(cell.membership = INPV_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), group.label = "membership")
INPV_ASD_SC40_INV <- mc_INV(cell.membership = INPV_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), group.label = "membership")
INPV_ASD_SC50_INV <- mc_INV(cell.membership = INPV_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INPV_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INPV_ASD_lncRNAs_data"]])), group.label = "membership")

# INPV_Normal
INPV_Normal_SC10_INV <- mc_INV(cell.membership = INPV_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), group.label = "membership")
INPV_Normal_SC20_INV <- mc_INV(cell.membership = INPV_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), group.label = "membership")
INPV_Normal_SC30_INV <- mc_INV(cell.membership = INPV_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), group.label = "membership")
INPV_Normal_SC40_INV <- mc_INV(cell.membership = INPV_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), group.label = "membership")
INPV_Normal_SC50_INV <- mc_INV(cell.membership = INPV_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INPV_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INPV_control_lncRNAs_data"]])), group.label = "membership")

# L4_ASD
L4_ASD_SC10_INV <- mc_INV(cell.membership = L4_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), group.label = "membership")
L4_ASD_SC20_INV <- mc_INV(cell.membership = L4_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), group.label = "membership")
L4_ASD_SC30_INV <- mc_INV(cell.membership = L4_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), group.label = "membership")
L4_ASD_SC40_INV <- mc_INV(cell.membership = L4_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), group.label = "membership")
L4_ASD_SC50_INV <- mc_INV(cell.membership = L4_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["L4_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["L4_ASD_lncRNAs_data"]])), group.label = "membership")

# L4_Normal
L4_Normal_SC10_INV <- mc_INV(cell.membership = L4_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), group.label = "membership")
L4_Normal_SC20_INV <- mc_INV(cell.membership = L4_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), group.label = "membership")
L4_Normal_SC30_INV <- mc_INV(cell.membership = L4_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), group.label = "membership")
L4_Normal_SC40_INV <- mc_INV(cell.membership = L4_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), group.label = "membership")
L4_Normal_SC50_INV <- mc_INV(cell.membership = L4_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["L4_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["L4_control_lncRNAs_data"]])), group.label = "membership")

# INSST_ASD
INSST_ASD_SC10_INV <- mc_INV(cell.membership = INSST_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), group.label = "membership")
INSST_ASD_SC20_INV <- mc_INV(cell.membership = INSST_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), group.label = "membership")
INSST_ASD_SC30_INV <- mc_INV(cell.membership = INSST_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), group.label = "membership")
INSST_ASD_SC40_INV <- mc_INV(cell.membership = INSST_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), group.label = "membership")
INSST_ASD_SC50_INV <- mc_INV(cell.membership = INSST_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["INSST_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["INSST_ASD_lncRNAs_data"]])), group.label = "membership")

# INSST_Normal
INSST_Normal_SC10_INV <- mc_INV(cell.membership = INSST_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), group.label = "membership")
INSST_Normal_SC20_INV <- mc_INV(cell.membership = INSST_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), group.label = "membership")
INSST_Normal_SC30_INV <- mc_INV(cell.membership = INSST_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), group.label = "membership")
INSST_Normal_SC40_INV <- mc_INV(cell.membership = INSST_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), group.label = "membership")
INSST_Normal_SC50_INV <- mc_INV(cell.membership = INSST_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["INSST_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["INSST_control_lncRNAs_data"]])), group.label = "membership")

# Neumat_ASD
Neumat_ASD_SC10_INV <- mc_INV(cell.membership = Neumat_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), group.label = "membership")
Neumat_ASD_SC20_INV <- mc_INV(cell.membership = Neumat_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), group.label = "membership")
Neumat_ASD_SC30_INV <- mc_INV(cell.membership = Neumat_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), group.label = "membership")
Neumat_ASD_SC40_INV <- mc_INV(cell.membership = Neumat_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), group.label = "membership")
Neumat_ASD_SC50_INV <- mc_INV(cell.membership = Neumat_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["Neumat_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["Neumat_ASD_lncRNAs_data"]])), group.label = "membership")

# Neumat_Normal
Neumat_Normal_SC10_INV <- mc_INV(cell.membership = Neumat_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), group.label = "membership")
Neumat_Normal_SC20_INV <- mc_INV(cell.membership = Neumat_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), group.label = "membership")
Neumat_Normal_SC30_INV <- mc_INV(cell.membership = Neumat_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), group.label = "membership")
Neumat_Normal_SC40_INV <- mc_INV(cell.membership = Neumat_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), group.label = "membership")
Neumat_Normal_SC50_INV <- mc_INV(cell.membership = Neumat_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["Neumat_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["Neumat_control_lncRNAs_data"]])), group.label = "membership")

# ASTPP_ASD
ASTPP_ASD_SC10_INV <- mc_INV(cell.membership = ASTPP_ASD_SC10_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), group.label = "membership")
ASTPP_ASD_SC20_INV <- mc_INV(cell.membership = ASTPP_ASD_SC20_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), group.label = "membership")
ASTPP_ASD_SC30_INV <- mc_INV(cell.membership = ASTPP_ASD_SC30_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), group.label = "membership")
ASTPP_ASD_SC40_INV <- mc_INV(cell.membership = ASTPP_ASD_SC40_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), group.label = "membership")
ASTPP_ASD_SC50_INV <- mc_INV(cell.membership = ASTPP_ASD_SC50_membership, sc.obj = CreateSeuratObject(rbind(ASD_celltype_mRNAs_list[["ASTPP_ASD_mRNAs_data"]], ASD_celltype_lncRNAs_list[["ASTPP_ASD_lncRNAs_data"]])), group.label = "membership")

# ASTPP_Normal
ASTPP_Normal_SC10_INV <- mc_INV(cell.membership = ASTPP_Normal_SC10_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), group.label = "membership")
ASTPP_Normal_SC20_INV <- mc_INV(cell.membership = ASTPP_Normal_SC20_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), group.label = "membership")
ASTPP_Normal_SC30_INV <- mc_INV(cell.membership = ASTPP_Normal_SC30_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), group.label = "membership")
ASTPP_Normal_SC40_INV <- mc_INV(cell.membership = ASTPP_Normal_SC40_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), group.label = "membership")
ASTPP_Normal_SC50_INV <- mc_INV(cell.membership = ASTPP_Normal_SC50_membership, sc.obj = CreateSeuratObject(rbind(Control_celltype_mRNAs_list[["ASTPP_control_mRNAs_data"]], Control_celltype_lncRNAs_list[["ASTPP_control_lncRNAs_data"]])), group.label = "membership")

############################################################################## Cell trajectory inference ##############################################################################
## Cell trajectory inference using SCORPIUS R package (https://CRAN.R-project.org/package=SCORPIUS)
# ASD and Normal
ASD_space <- reduce_dimensionality(t(as.matrix(ASD_SC30_GE)), ndim =2)
Normal_space <- reduce_dimensionality(t(as.matrix(Normal_SC30_GE)), ndim = 2)

ASD_model <- infer_trajectory(ASD_space)
Normal_model <- infer_trajectory(Normal_space)

ASD_model_data <- t(as.matrix(ASD_SC30_GE))[match(names(sort(ASD_model$time)), colnames(ASD_SC30_GE)), ]
Normal_model_data <- t(as.matrix(Normal_SC30_GE))[match(names(sort(Normal_model$time)), colnames(Normal_SC30_GE)), ]

# ACC_ASD, ACC_Normal, PFC_ASD, PFC_Normal
ACC_ASD_space <- reduce_dimensionality(t(as.matrix(ACC_ASD_SC30_GE)), ndim =2)
ACC_Normal_space <- reduce_dimensionality(t(as.matrix(ACC_Normal_SC30_GE)), ndim = 2)
PFC_ASD_space <- reduce_dimensionality(t(as.matrix(PFC_ASD_SC30_GE)), ndim =2)
PFC_Normal_space <- reduce_dimensionality(t(as.matrix(PFC_Normal_SC30_GE)), ndim = 2)

ACC_ASD_model <- infer_trajectory(ACC_ASD_space)
ACC_Normal_model <- infer_trajectory(ACC_Normal_space)
PFC_ASD_model <- infer_trajectory(PFC_ASD_space)
PFC_Normal_model <- infer_trajectory(PFC_Normal_space)

ACC_ASD_model_data <- t(as.matrix(ACC_ASD_SC30_GE))[match(names(sort(ACC_ASD_model$time)), colnames(ACC_ASD_SC30_GE)), ]
ACC_Normal_model_data <- t(as.matrix(ACC_Normal_SC30_GE))[match(names(sort(ACC_Normal_model$time)), colnames(ACC_Normal_SC30_GE)), ]
PFC_ASD_model_data <- t(as.matrix(PFC_ASD_SC30_GE))[match(names(sort(PFC_ASD_model$time)), colnames(PFC_ASD_SC30_GE)), ]
PFC_Normal_model_data <- t(as.matrix(PFC_Normal_SC30_GE))[match(names(sort(PFC_Normal_model$time)), colnames(PFC_Normal_SC30_GE)), ]

# Lower18_ASD, Lower18_Normal, Larger18_ASD, Larger18_Normal
Lower18_ASD_space <- reduce_dimensionality(t(as.matrix(Lower18_ASD_SC30_GE)), ndim =2)
Lower18_Normal_space <- reduce_dimensionality(t(as.matrix(Lower18_Normal_SC30_GE)), ndim = 2)
Larger18_ASD_space <- reduce_dimensionality(t(as.matrix(Larger18_ASD_SC30_GE)), ndim =2)
Larger18_Normal_space <- reduce_dimensionality(t(as.matrix(Larger18_Normal_SC30_GE)), ndim = 2)

Lower18_ASD_model <- infer_trajectory(Lower18_ASD_space)
Lower18_Normal_model <- infer_trajectory(Lower18_Normal_space)
Larger18_ASD_model <- infer_trajectory(Larger18_ASD_space)
Larger18_Normal_model <- infer_trajectory(Larger18_Normal_space)

Lower18_ASD_model_data <- t(as.matrix(Lower18_ASD_SC30_GE))[match(names(sort(Lower18_ASD_model$time)), colnames(Lower18_ASD_SC30_GE)), ]
Lower18_Normal_model_data <- t(as.matrix(Lower18_Normal_SC30_GE))[match(names(sort(Lower18_Normal_model$time)), colnames(Lower18_Normal_SC30_GE)), ]
Larger18_ASD_model_data <- t(as.matrix(Larger18_ASD_SC30_GE))[match(names(sort(Larger18_ASD_model$time)), colnames(Larger18_ASD_SC30_GE)), ]
Larger18_Normal_model_data <- t(as.matrix(Larger18_Normal_SC30_GE))[match(names(sort(Larger18_Normal_model$time)), colnames(Larger18_Normal_SC30_GE)), ]

# Male_ASD, Male_Normal, Female_ASD, Female_Normal
Male_ASD_space <- reduce_dimensionality(t(as.matrix(Male_ASD_SC30_GE)), ndim =2)
Male_Normal_space <- reduce_dimensionality(t(as.matrix(Male_Normal_SC30_GE)), ndim = 2)
Female_ASD_space <- reduce_dimensionality(t(as.matrix(Female_ASD_SC30_GE)), ndim =2)
Female_Normal_space <- reduce_dimensionality(t(as.matrix(Female_Normal_SC30_GE)), ndim = 2)

Male_ASD_model <- infer_trajectory(Male_ASD_space)
Male_Normal_model <- infer_trajectory(Male_Normal_space)
Female_ASD_model <- infer_trajectory(Female_ASD_space)
Female_Normal_model <- infer_trajectory(Female_Normal_space)

Male_ASD_model_data <- t(as.matrix(Male_ASD_SC30_GE))[match(names(sort(Male_ASD_model$time)), colnames(Male_ASD_SC30_GE)), ]
Male_Normal_model_data <- t(as.matrix(Male_Normal_SC30_GE))[match(names(sort(Male_Normal_model$time)), colnames(Male_Normal_SC30_GE)), ]
Female_ASD_model_data <- t(as.matrix(Female_ASD_SC30_GE))[match(names(sort(Female_ASD_model$time)), colnames(Female_ASD_SC30_GE)), ]
Female_Normal_model_data <- t(as.matrix(Female_Normal_SC30_GE))[match(names(sort(Female_Normal_model$time)), colnames(Female_Normal_SC30_GE)), ]

# NeuNRGNII_ASD, NeuNRGNII_Normal
NeuNRGNII_ASD_space <- reduce_dimensionality(t(as.matrix(NeuNRGNII_ASD_SC30_GE)), ndim =2)
NeuNRGNII_Normal_space <- reduce_dimensionality(t(as.matrix(NeuNRGNII_Normal_SC30_GE)), ndim = 2)

NeuNRGNII_ASD_model <- infer_trajectory(NeuNRGNII_ASD_space)
NeuNRGNII_Normal_model <- infer_trajectory(NeuNRGNII_Normal_space)

NeuNRGNII_ASD_model_data <- t(as.matrix(NeuNRGNII_ASD_SC30_GE))[match(names(sort(NeuNRGNII_ASD_model$time)), colnames(NeuNRGNII_ASD_SC30_GE)), ]
NeuNRGNII_Normal_model_data <- t(as.matrix(NeuNRGNII_Normal_SC30_GE))[match(names(sort(NeuNRGNII_Normal_model$time)), colnames(NeuNRGNII_Normal_SC30_GE)), ]

# L56_ASD, L56_Normal
L56_ASD_space <- reduce_dimensionality(t(as.matrix(L56_ASD_SC30_GE)), ndim =2)
L56_Normal_space <- reduce_dimensionality(t(as.matrix(L56_Normal_SC30_GE)), ndim = 2)

L56_ASD_model <- infer_trajectory(L56_ASD_space)
L56_Normal_model <- infer_trajectory(L56_Normal_space)

L56_ASD_model_data <- t(as.matrix(L56_ASD_SC30_GE))[match(names(sort(L56_ASD_model$time)), colnames(L56_ASD_SC30_GE)), ]
L56_Normal_model_data <- t(as.matrix(L56_Normal_SC30_GE))[match(names(sort(L56_Normal_model$time)), colnames(L56_Normal_SC30_GE)), ]

# Oligodendrocytes_ASD, Oligodendrocytes_Normal
Oligodendrocytes_ASD_space <- reduce_dimensionality(t(as.matrix(Oligodendrocytes_ASD_SC30_GE)), ndim =2)
Oligodendrocytes_Normal_space <- reduce_dimensionality(t(as.matrix(Oligodendrocytes_Normal_SC30_GE)), ndim = 2)

Oligodendrocytes_ASD_model <- infer_trajectory(Oligodendrocytes_ASD_space)
Oligodendrocytes_Normal_model <- infer_trajectory(Oligodendrocytes_Normal_space)

Oligodendrocytes_ASD_model_data <- t(as.matrix(Oligodendrocytes_ASD_SC30_GE))[match(names(sort(Oligodendrocytes_ASD_model$time)), colnames(Oligodendrocytes_ASD_SC30_GE)), ]
Oligodendrocytes_Normal_model_data <- t(as.matrix(Oligodendrocytes_Normal_SC30_GE))[match(names(sort(Oligodendrocytes_Normal_model$time)), colnames(Oligodendrocytes_Normal_SC30_GE)), ]

# OPC_ASD, OPC_Normal
OPC_ASD_space <- reduce_dimensionality(t(as.matrix(OPC_ASD_SC30_GE)), ndim =2)
OPC_Normal_space <- reduce_dimensionality(t(as.matrix(OPC_Normal_SC30_GE)), ndim = 2)

OPC_ASD_model <- infer_trajectory(OPC_ASD_space)
OPC_Normal_model <- infer_trajectory(OPC_Normal_space)

OPC_ASD_model_data <- t(as.matrix(OPC_ASD_SC30_GE))[match(names(sort(OPC_ASD_model$time)), colnames(OPC_ASD_SC30_GE)), ]
OPC_Normal_model_data <- t(as.matrix(OPC_Normal_SC30_GE))[match(names(sort(OPC_Normal_model$time)), colnames(OPC_Normal_SC30_GE)), ]

# ASTFB_ASD, ASTFB_Normal
ASTFB_ASD_space <- reduce_dimensionality(t(as.matrix(ASTFB_ASD_SC30_GE)), ndim =2)
ASTFB_Normal_space <- reduce_dimensionality(t(as.matrix(ASTFB_Normal_SC30_GE)), ndim = 2)

ASTFB_ASD_model <- infer_trajectory(ASTFB_ASD_space)
ASTFB_Normal_model <- infer_trajectory(ASTFB_Normal_space)

ASTFB_ASD_model_data <- t(as.matrix(ASTFB_ASD_SC30_GE))[match(names(sort(ASTFB_ASD_model$time)), colnames(ASTFB_ASD_SC30_GE)), ]
ASTFB_Normal_model_data <- t(as.matrix(ASTFB_Normal_SC30_GE))[match(names(sort(ASTFB_Normal_model$time)), colnames(ASTFB_Normal_SC30_GE)), ]

# Endothelial_ASD, Endothelial_Normal
Endothelial_ASD_space <- reduce_dimensionality(t(as.matrix(Endothelial_ASD_SC30_GE)), ndim =2)
Endothelial_Normal_space <- reduce_dimensionality(t(as.matrix(Endothelial_Normal_SC30_GE)), ndim = 2)

Endothelial_ASD_model <- infer_trajectory(Endothelial_ASD_space)
Endothelial_Normal_model <- infer_trajectory(Endothelial_Normal_space)

Endothelial_ASD_model_data <- t(as.matrix(Endothelial_ASD_SC30_GE))[match(names(sort(Endothelial_ASD_model$time)), colnames(Endothelial_ASD_SC30_GE)), ]
Endothelial_Normal_model_data <- t(as.matrix(Endothelial_Normal_SC30_GE))[match(names(sort(Endothelial_Normal_model$time)), colnames(Endothelial_Normal_SC30_GE)), ]

# Microglia_ASD, Microglia_Normal
Microglia_ASD_space <- reduce_dimensionality(t(as.matrix(Microglia_ASD_SC30_GE)), ndim =2)
Microglia_Normal_space <- reduce_dimensionality(t(as.matrix(Microglia_Normal_SC30_GE)), ndim = 2)

Microglia_ASD_model <- infer_trajectory(Microglia_ASD_space)
Microglia_Normal_model <- infer_trajectory(Microglia_Normal_space)

Microglia_ASD_model_data <- t(as.matrix(Microglia_ASD_SC30_GE))[match(names(sort(Microglia_ASD_model$time)), colnames(Microglia_ASD_SC30_GE)), ]
Microglia_Normal_model_data <- t(as.matrix(Microglia_Normal_SC30_GE))[match(names(sort(Microglia_Normal_model$time)), colnames(Microglia_Normal_SC30_GE)), ]

# NeuNRGNI_ASD, NeuNRGNI_Normal
NeuNRGNI_ASD_space <- reduce_dimensionality(t(as.matrix(NeuNRGNI_ASD_SC30_GE)), ndim =2)
NeuNRGNI_Normal_space <- reduce_dimensionality(t(as.matrix(NeuNRGNI_Normal_SC30_GE)), ndim = 2)

NeuNRGNI_ASD_model <- infer_trajectory(NeuNRGNI_ASD_space)
NeuNRGNI_Normal_model <- infer_trajectory(NeuNRGNI_Normal_space)

NeuNRGNI_ASD_model_data <- t(as.matrix(NeuNRGNI_ASD_SC30_GE))[match(names(sort(NeuNRGNI_ASD_model$time)), colnames(NeuNRGNI_ASD_SC30_GE)), ]
NeuNRGNI_Normal_model_data <- t(as.matrix(NeuNRGNI_Normal_SC30_GE))[match(names(sort(NeuNRGNI_Normal_model$time)), colnames(NeuNRGNI_Normal_SC30_GE)), ]

# INVIP_ASD, INVIP_Normal
INVIP_ASD_space <- reduce_dimensionality(t(as.matrix(INVIP_ASD_SC30_GE)), ndim =2)
INVIP_Normal_space <- reduce_dimensionality(t(as.matrix(INVIP_Normal_SC30_GE)), ndim = 2)

INVIP_ASD_model <- infer_trajectory(INVIP_ASD_space)
INVIP_Normal_model <- infer_trajectory(INVIP_Normal_space)

INVIP_ASD_model_data <- t(as.matrix(INVIP_ASD_SC30_GE))[match(names(sort(INVIP_ASD_model$time)), colnames(INVIP_ASD_SC30_GE)), ]
INVIP_Normal_model_data <- t(as.matrix(INVIP_Normal_SC30_GE))[match(names(sort(INVIP_Normal_model$time)), colnames(INVIP_Normal_SC30_GE)), ]

# L56CC_ASD, L56CC_Normal
L56CC_ASD_space <- reduce_dimensionality(t(as.matrix(L56CC_ASD_SC30_GE)), ndim =2)
L56CC_Normal_space <- reduce_dimensionality(t(as.matrix(L56CC_Normal_SC30_GE)), ndim = 2)

L56CC_ASD_model <- infer_trajectory(L56CC_ASD_space)
L56CC_Normal_model <- infer_trajectory(L56CC_Normal_space)

L56CC_ASD_model_data <- t(as.matrix(L56CC_ASD_SC30_GE))[match(names(sort(L56CC_ASD_model$time)), colnames(L56CC_ASD_SC30_GE)), ]
L56CC_Normal_model_data <- t(as.matrix(L56CC_Normal_SC30_GE))[match(names(sort(L56CC_Normal_model$time)), colnames(L56CC_Normal_SC30_GE)), ]

# INSV2C_ASD, INSV2C_Normal
INSV2C_ASD_space <- reduce_dimensionality(t(as.matrix(INSV2C_ASD_SC30_GE)), ndim =2)
INSV2C_Normal_space <- reduce_dimensionality(t(as.matrix(INSV2C_Normal_SC30_GE)), ndim = 2)

INSV2C_ASD_model <- infer_trajectory(INSV2C_ASD_space)
INSV2C_Normal_model <- infer_trajectory(INSV2C_Normal_space)

INSV2C_ASD_model_data <- t(as.matrix(INSV2C_ASD_SC30_GE))[match(names(sort(INSV2C_ASD_model$time)), colnames(INSV2C_ASD_SC30_GE)), ]
INSV2C_Normal_model_data <- t(as.matrix(INSV2C_Normal_SC30_GE))[match(names(sort(INSV2C_Normal_model$time)), colnames(INSV2C_Normal_SC30_GE)), ]

# L23_ASD, L23_Normal
L23_ASD_space <- reduce_dimensionality(t(as.matrix(L23_ASD_SC30_GE)), ndim =2)
L23_Normal_space <- reduce_dimensionality(t(as.matrix(L23_Normal_SC30_GE)), ndim = 2)

L23_ASD_model <- infer_trajectory(L23_ASD_space)
L23_Normal_model <- infer_trajectory(L23_Normal_space)

L23_ASD_model_data <- t(as.matrix(L23_ASD_SC30_GE))[match(names(sort(L23_ASD_model$time)), colnames(L23_ASD_SC30_GE)), ]
L23_Normal_model_data <- t(as.matrix(L23_Normal_SC30_GE))[match(names(sort(L23_Normal_model$time)), colnames(L23_Normal_SC30_GE)), ]

# INPV_ASD, INPV_Normal
INPV_ASD_space <- reduce_dimensionality(t(as.matrix(INPV_ASD_SC30_GE)), ndim =2)
INPV_Normal_space <- reduce_dimensionality(t(as.matrix(INPV_Normal_SC30_GE)), ndim = 2)

INPV_ASD_model <- infer_trajectory(INPV_ASD_space)
INPV_Normal_model <- infer_trajectory(INPV_Normal_space)

INPV_ASD_model_data <- t(as.matrix(INPV_ASD_SC30_GE))[match(names(sort(INPV_ASD_model$time)), colnames(INPV_ASD_SC30_GE)), ]
INPV_Normal_model_data <- t(as.matrix(INPV_Normal_SC30_GE))[match(names(sort(INPV_Normal_model$time)), colnames(INPV_Normal_SC30_GE)), ]

# L4_ASD, L4_Normal
L4_ASD_space <- reduce_dimensionality(t(as.matrix(L4_ASD_SC30_GE)), ndim =2)
L4_Normal_space <- reduce_dimensionality(t(as.matrix(L4_Normal_SC30_GE)), ndim = 2)

L4_ASD_model <- infer_trajectory(L4_ASD_space)
L4_Normal_model <- infer_trajectory(L4_Normal_space)

L4_ASD_model_data <- t(as.matrix(L4_ASD_SC30_GE))[match(names(sort(L4_ASD_model$time)), colnames(L4_ASD_SC30_GE)), ]
L4_Normal_model_data <- t(as.matrix(L4_Normal_SC30_GE))[match(names(sort(L4_Normal_model$time)), colnames(L4_Normal_SC30_GE)), ]

# INSST_ASD, INSST_Normal
INSST_ASD_space <- reduce_dimensionality(t(as.matrix(INSST_ASD_SC30_GE)), ndim =2)
INSST_Normal_space <- reduce_dimensionality(t(as.matrix(INSST_Normal_SC30_GE)), ndim = 2)

INSST_ASD_model <- infer_trajectory(INSST_ASD_space)
INSST_Normal_model <- infer_trajectory(INSST_Normal_space)

INSST_ASD_model_data <- t(as.matrix(INSST_ASD_SC30_GE))[match(names(sort(INSST_ASD_model$time)), colnames(INSST_ASD_SC30_GE)), ]
INSST_Normal_model_data <- t(as.matrix(INSST_Normal_SC30_GE))[match(names(sort(INSST_Normal_model$time)), colnames(INSST_Normal_SC30_GE)), ]

# Neumat_ASD, Neumat_Normal
Neumat_ASD_space <- reduce_dimensionality(t(as.matrix(Neumat_ASD_SC30_GE)), ndim =2)
Neumat_Normal_space <- reduce_dimensionality(t(as.matrix(Neumat_Normal_SC30_GE)), ndim = 2)

Neumat_ASD_model <- infer_trajectory(Neumat_ASD_space)
Neumat_Normal_model <- infer_trajectory(Neumat_Normal_space)

Neumat_ASD_model_data <- t(as.matrix(Neumat_ASD_SC30_GE))[match(names(sort(Neumat_ASD_model$time)), colnames(Neumat_ASD_SC30_GE)), ]
Neumat_Normal_model_data <- t(as.matrix(Neumat_Normal_SC30_GE))[match(names(sort(Neumat_Normal_model$time)), colnames(Neumat_Normal_SC30_GE)), ]

# ASTPP_ASD, ASTPP_Normal
ASTPP_ASD_space <- reduce_dimensionality(t(as.matrix(ASTPP_ASD_SC30_GE)), ndim =2)
ASTPP_Normal_space <- reduce_dimensionality(t(as.matrix(ASTPP_Normal_SC30_GE)), ndim = 2)

ASTPP_ASD_model <- infer_trajectory(ASTPP_ASD_space)
ASTPP_Normal_model <- infer_trajectory(ASTPP_Normal_space)

ASTPP_ASD_model_data <- t(as.matrix(ASTPP_ASD_SC30_GE))[match(names(sort(ASTPP_ASD_model$time)), colnames(ASTPP_ASD_SC30_GE)), ]
ASTPP_Normal_model_data <- t(as.matrix(ASTPP_Normal_SC30_GE))[match(names(sort(ASTPP_Normal_model$time)), colnames(ASTPP_Normal_SC30_GE)), ]

################################################################################  Clean application ############################################################################
# ASD and Normal
ASD_res_darkcausality <- darkcausality_parallel(ASD_model_data, cause, effect, num.cores = 48)
Normal_res_darkcausality <- darkcausality_parallel(Normal_model_data, cause, effect, num.cores = 48)

# ACC_ASD, ACC_Normal, PFC_ASD, PFC_Normal
ACC_ASD_res_darkcausality <- darkcausality_parallel(ACC_ASD_model_data, cause, effect, num.cores = 48)
ACC_Normal_res_darkcausality <- darkcausality_parallel(ACC_Normal_model_data, cause, effect, num.cores = 48)
PFC_ASD_res_darkcausality <- darkcausality_parallel(PFC_ASD_model_data, cause, effect, num.cores = 48)
PFC_Normal_res_darkcausality <- darkcausality_parallel(PFC_Normal_model_data, cause, effect, num.cores = 48)

# Lower18_ASD, Lower18_Normal, Larger18_ASD, Larger18_Normal
Lower18_ASD_res_darkcausality <- darkcausality_parallel(Lower18_ASD_model_data, cause, effect, num.cores = 48)
Lower18_Normal_res_darkcausality <- darkcausality_parallel(Lower18_Normal_model_data, cause, effect, num.cores = 48)
Larger18_ASD_res_darkcausality <- darkcausality_parallel(Larger18_ASD_model_data, cause, effect, num.cores = 48)
Larger18_Normal_res_darkcausality <- darkcausality_parallel(Larger18_Normal_model_data, cause, effect, num.cores = 48)

# Male_ASD, Male_Normal, Female_ASD, Female_Normal
Male_ASD_res_darkcausality <- darkcausality_parallel(Male_ASD_model_data, cause, effect, num.cores = 48)
Male_Normal_res_darkcausality <- darkcausality_parallel(Male_Normal_model_data, cause, effect, num.cores = 48)
Female_ASD_res_darkcausality <- darkcausality_parallel(Female_ASD_model_data, cause, effect, num.cores = 48)
Female_Normal_res_darkcausality <- darkcausality_parallel(Female_Normal_model_data, cause, effect, num.cores = 48)

# NeuNRGNII_ASD, NeuNRGNII_Normal
NeuNRGNII_ASD_res_darkcausality <- darkcausality_parallel(NeuNRGNII_ASD_model_data, cause, effect, num.cores = 48)
NeuNRGNII_Normal_res_darkcausality <- darkcausality_parallel(NeuNRGNII_Normal_model_data, cause, effect, num.cores = 48)

# L56_ASD, L56_Normal
L56_ASD_res_darkcausality <- darkcausality_parallel(L56_ASD_model_data, cause, effect, num.cores = 48)
L56_Normal_res_darkcausality <- darkcausality_parallel(L56_Normal_model_data, cause, effect, num.cores = 48)

# Oligodendrocytes_ASD, Oligodendrocytes_Normal
Oligodendrocytes_ASD_res_darkcausality <- darkcausality_parallel(Oligodendrocytes_ASD_model_data, cause, effect, num.cores = 48)
Oligodendrocytes_Normal_res_darkcausality <- darkcausality_parallel(Oligodendrocytes_Normal_model_data, cause, effect, num.cores = 48)

# OPC_ASD, OPC_Normal
OPC_ASD_res_darkcausality <- darkcausality_parallel(OPC_ASD_model_data, cause, effect, num.cores = 48)
OPC_Normal_res_darkcausality <- darkcausality_parallel(OPC_Normal_model_data, cause, effect, num.cores = 48)
 
# ASTFB_ASD, ASTFB_Normal
ASTFB_ASD_res_darkcausality <- darkcausality_parallel(ASTFB_ASD_model_data, cause, effect, num.cores = 48)
ASTFB_Normal_res_darkcausality <- darkcausality_parallel(ASTFB_Normal_model_data, cause, effect, num.cores = 48)

# Endothelial_ASD, Endothelial_Normal
Endothelial_ASD_res_darkcausality <- darkcausality_parallel(Endothelial_ASD_model_data, cause, effect, num.cores = 48)
Endothelial_Normal_res_darkcausality <- darkcausality_parallel(Endothelial_Normal_model_data, cause, effect, num.cores = 48)

# Microglia_ASD, Microglia_Normal
Microglia_ASD_res_darkcausality <- darkcausality_parallel(Microglia_ASD_model_data, cause, effect, num.cores = 48)
Microglia_Normal_res_darkcausality <- darkcausality_parallel(Microglia_Normal_model_data, cause, effect, num.cores = 48)

# NeuNRGNI_ASD, NeuNRGNI_Normal
NeuNRGNI_ASD_res_darkcausality <- darkcausality_parallel(NeuNRGNI_ASD_model_data, cause, effect, num.cores = 48)
NeuNRGNI_Normal_res_darkcausality <- darkcausality_parallel(NeuNRGNI_Normal_model_data, cause, effect, num.cores = 48)

# INVIP_ASD, INVIP_Normal
INVIP_ASD_res_darkcausality <- darkcausality_parallel(INVIP_ASD_model_data, cause, effect, num.cores = 48)
INVIP_Normal_res_darkcausality <- darkcausality_parallel(INVIP_Normal_model_data, cause, effect, num.cores = 48)

# L56CC_ASD, L56CC_Normal
L56CC_ASD_res_darkcausality <- darkcausality_parallel(L56CC_ASD_model_data, cause, effect, num.cores = 48)
L56CC_Normal_res_darkcausality <- darkcausality_parallel(L56CC_Normal_model_data, cause, effect, num.cores = 48)

# INSV2C_ASD, INSV2C_Normal
INSV2C_ASD_res_darkcausality <- darkcausality_parallel(INSV2C_ASD_model_data, cause, effect, num.cores = 48)
INSV2C_Normal_res_darkcausality <- darkcausality_parallel(INSV2C_Normal_model_data, cause, effect, num.cores = 48)

# L23_ASD, L23_Normal
L23_ASD_res_darkcausality <- darkcausality_parallel(L23_ASD_model_data, cause, effect, num.cores = 48)
L23_Normal_res_darkcausality <- darkcausality_parallel(L23_Normal_model_data, cause, effect, num.cores = 48)

# INPV_ASD, INPV_Normal
INPV_ASD_res_darkcausality <- darkcausality_parallel(INPV_ASD_model_data, cause, effect, num.cores = 48)
INPV_Normal_res_darkcausality <- darkcausality_parallel(INPV_Normal_model_data, cause, effect, num.cores = 48)

# L4_ASD, L4_Normal
L4_ASD_res_darkcausality <- darkcausality_parallel(L4_ASD_model_data, cause, effect, num.cores = 48)
L4_Normal_res_darkcausality <- darkcausality_parallel(L4_Normal_model_data, cause, effect, num.cores = 48)

# INSST_ASD, INSST_Normal
INSST_ASD_res_darkcausality <- darkcausality_parallel(INSST_ASD_model_data, cause, effect, num.cores = 48)
INSST_Normal_res_darkcausality <- darkcausality_parallel(INSST_Normal_model_data, cause, effect, num.cores = 48)

# Neumat_ASD, Neumat_Normal
Neumat_ASD_res_darkcausality <- darkcausality_parallel(Neumat_ASD_model_data, cause, effect, num.cores = 48)
Neumat_Normal_res_darkcausality <- darkcausality_parallel(Neumat_Normal_model_data, cause, effect, num.cores = 48)

# ASTPP_ASD, ASTPP_Normal
ASTPP_ASD_res_darkcausality <- darkcausality_parallel(ASTPP_ASD_model_data, cause, effect, num.cores = 48)
ASTPP_Normal_res_darkcausality <- darkcausality_parallel(ASTPP_Normal_model_data, cause, effect, num.cores = 48)

## ASD and Normal
ASD_res_darkcausality[[2]] <- ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
ASD_res_darkcausality[[3]] <- ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
ASD_res_darkcausality[[4]] <- ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
ASD_res_darkcausality_sig <- sigPC(ASD_res_darkcausality[[2]], ASD_res_darkcausality[[3]], ASD_res_darkcausality[[4]])
Normal_res_darkcausality[[2]] <- Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Normal_res_darkcausality[[3]] <- Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Normal_res_darkcausality[[4]] <- Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Normal_res_darkcausality_sig <- sigPC(Normal_res_darkcausality[[2]], Normal_res_darkcausality[[3]], Normal_res_darkcausality[[4]])
ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(ASD_res_darkcausality_sig[[1]] + ASD_res_darkcausality_sig[[2]] + ASD_res_darkcausality_sig[[3]]))
Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Normal_res_darkcausality_sig[[1]] + Normal_res_darkcausality_sig[[2]] + Normal_res_darkcausality_sig[[3]]))

# ACC_ASD, ACC_Normal, PFC_ASD, PFC_Normal
ACC_ASD_res_darkcausality[[2]] <- ACC_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
ACC_ASD_res_darkcausality[[3]] <- ACC_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
ACC_ASD_res_darkcausality[[4]] <- ACC_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
ACC_ASD_res_darkcausality_sig <- sigPC(ACC_ASD_res_darkcausality[[2]], ACC_ASD_res_darkcausality[[3]], ACC_ASD_res_darkcausality[[4]])
ACC_Normal_res_darkcausality[[2]] <- ACC_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
ACC_Normal_res_darkcausality[[3]] <- ACC_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
ACC_Normal_res_darkcausality[[4]] <- ACC_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
ACC_Normal_res_darkcausality_sig <- sigPC(ACC_Normal_res_darkcausality[[2]], ACC_Normal_res_darkcausality[[3]], ACC_Normal_res_darkcausality[[4]])
PFC_ASD_res_darkcausality[[2]] <- PFC_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
PFC_ASD_res_darkcausality[[3]] <- PFC_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
PFC_ASD_res_darkcausality[[4]] <- PFC_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
PFC_ASD_res_darkcausality_sig <- sigPC(PFC_ASD_res_darkcausality[[2]], PFC_ASD_res_darkcausality[[3]], PFC_ASD_res_darkcausality[[4]])
PFC_Normal_res_darkcausality[[2]] <- PFC_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
PFC_Normal_res_darkcausality[[3]] <- PFC_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
PFC_Normal_res_darkcausality[[4]] <- PFC_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
PFC_Normal_res_darkcausality_sig <- sigPC(PFC_Normal_res_darkcausality[[2]], PFC_Normal_res_darkcausality[[3]], PFC_Normal_res_darkcausality[[4]])
ACC_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(ACC_ASD_res_darkcausality_sig[[1]] + ACC_ASD_res_darkcausality_sig[[2]] + ACC_ASD_res_darkcausality_sig[[3]]))
ACC_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(ACC_Normal_res_darkcausality_sig[[1]] + ACC_Normal_res_darkcausality_sig[[2]] + ACC_Normal_res_darkcausality_sig[[3]]))
PFC_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(PFC_ASD_res_darkcausality_sig[[1]] + PFC_ASD_res_darkcausality_sig[[2]] + PFC_ASD_res_darkcausality_sig[[3]]))
PFC_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(PFC_Normal_res_darkcausality_sig[[1]] + PFC_Normal_res_darkcausality_sig[[2]] + PFC_Normal_res_darkcausality_sig[[3]]))

# Lower18_ASD, Lower18_Normal, Larger18_ASD, Larger18_Normal
Lower18_ASD_res_darkcausality[[2]] <- Lower18_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Lower18_ASD_res_darkcausality[[3]] <- Lower18_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Lower18_ASD_res_darkcausality[[4]] <- Lower18_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Lower18_ASD_res_darkcausality_sig <- sigPC(Lower18_ASD_res_darkcausality[[2]], Lower18_ASD_res_darkcausality[[3]], Lower18_ASD_res_darkcausality[[4]])
Lower18_Normal_res_darkcausality[[2]] <- Lower18_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Lower18_Normal_res_darkcausality[[3]] <- Lower18_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Lower18_Normal_res_darkcausality[[4]] <- Lower18_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Lower18_Normal_res_darkcausality_sig <- sigPC(Lower18_Normal_res_darkcausality[[2]], Lower18_Normal_res_darkcausality[[3]], Lower18_Normal_res_darkcausality[[4]])
Larger18_ASD_res_darkcausality[[2]] <- Larger18_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Larger18_ASD_res_darkcausality[[3]] <- Larger18_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Larger18_ASD_res_darkcausality[[4]] <- Larger18_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Larger18_ASD_res_darkcausality_sig <- sigPC(Larger18_ASD_res_darkcausality[[2]], Larger18_ASD_res_darkcausality[[3]], Larger18_ASD_res_darkcausality[[4]])
Larger18_Normal_res_darkcausality[[2]] <- Larger18_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Larger18_Normal_res_darkcausality[[3]] <- Larger18_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Larger18_Normal_res_darkcausality[[4]] <- Larger18_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Larger18_Normal_res_darkcausality_sig <- sigPC(Larger18_Normal_res_darkcausality[[2]], Larger18_Normal_res_darkcausality[[3]], Larger18_Normal_res_darkcausality[[4]])
Lower18_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Lower18_ASD_res_darkcausality_sig[[1]] + Lower18_ASD_res_darkcausality_sig[[2]] + Lower18_ASD_res_darkcausality_sig[[3]]))
Lower18_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Lower18_Normal_res_darkcausality_sig[[1]] + Lower18_Normal_res_darkcausality_sig[[2]] + Lower18_Normal_res_darkcausality_sig[[3]]))
Larger18_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Larger18_ASD_res_darkcausality_sig[[1]] + Larger18_ASD_res_darkcausality_sig[[2]] + Larger18_ASD_res_darkcausality_sig[[3]]))
Larger18_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Larger18_Normal_res_darkcausality_sig[[1]] + Larger18_Normal_res_darkcausality_sig[[2]] + Larger18_Normal_res_darkcausality_sig[[3]]))

# Male_ASD, Male_Normal, Female_ASD, Female_Normal
Male_ASD_res_darkcausality[[2]] <- Male_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Male_ASD_res_darkcausality[[3]] <- Male_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Male_ASD_res_darkcausality[[4]] <- Male_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Male_ASD_res_darkcausality_sig <- sigPC(Male_ASD_res_darkcausality[[2]], Male_ASD_res_darkcausality[[3]], Male_ASD_res_darkcausality[[4]])
Male_Normal_res_darkcausality[[2]] <- Male_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Male_Normal_res_darkcausality[[3]] <- Male_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Male_Normal_res_darkcausality[[4]] <- Male_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Male_Normal_res_darkcausality_sig <- sigPC(Male_Normal_res_darkcausality[[2]], Male_Normal_res_darkcausality[[3]], Male_Normal_res_darkcausality[[4]])
Female_ASD_res_darkcausality[[2]] <- Female_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Female_ASD_res_darkcausality[[3]] <- Female_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Female_ASD_res_darkcausality[[4]] <- Female_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Female_ASD_res_darkcausality_sig <- sigPC(Female_ASD_res_darkcausality[[2]], Female_ASD_res_darkcausality[[3]], Female_ASD_res_darkcausality[[4]])
Female_Normal_res_darkcausality[[2]] <- Female_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Female_Normal_res_darkcausality[[3]] <- Female_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Female_Normal_res_darkcausality[[4]] <- Female_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Female_Normal_res_darkcausality_sig <- sigPC(Female_Normal_res_darkcausality[[2]], Female_Normal_res_darkcausality[[3]], Female_Normal_res_darkcausality[[4]])
Male_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Male_ASD_res_darkcausality_sig[[1]] + Male_ASD_res_darkcausality_sig[[2]] + Male_ASD_res_darkcausality_sig[[3]]))
Male_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Male_Normal_res_darkcausality_sig[[1]] + Male_Normal_res_darkcausality_sig[[2]] + Male_Normal_res_darkcausality_sig[[3]]))
Female_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Female_ASD_res_darkcausality_sig[[1]] + Female_ASD_res_darkcausality_sig[[2]] + Female_ASD_res_darkcausality_sig[[3]]))
Female_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Female_Normal_res_darkcausality_sig[[1]] + Female_Normal_res_darkcausality_sig[[2]] + Female_Normal_res_darkcausality_sig[[3]]))

# NeuNRGNII_ASD, NeuNRGNII_Normal
NeuNRGNII_ASD_res_darkcausality[[2]] <- NeuNRGNII_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
NeuNRGNII_ASD_res_darkcausality[[3]] <- NeuNRGNII_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
NeuNRGNII_ASD_res_darkcausality[[4]] <- NeuNRGNII_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
NeuNRGNII_ASD_res_darkcausality_sig <- sigPC(NeuNRGNII_ASD_res_darkcausality[[2]], NeuNRGNII_ASD_res_darkcausality[[3]], NeuNRGNII_ASD_res_darkcausality[[4]])
NeuNRGNII_Normal_res_darkcausality[[2]] <- NeuNRGNII_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
NeuNRGNII_Normal_res_darkcausality[[3]] <- NeuNRGNII_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
NeuNRGNII_Normal_res_darkcausality[[4]] <- NeuNRGNII_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
NeuNRGNII_Normal_res_darkcausality_sig <- sigPC(NeuNRGNII_Normal_res_darkcausality[[2]], NeuNRGNII_Normal_res_darkcausality[[3]], NeuNRGNII_Normal_res_darkcausality[[4]])
NeuNRGNII_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(NeuNRGNII_ASD_res_darkcausality_sig[[1]] + NeuNRGNII_ASD_res_darkcausality_sig[[2]] + NeuNRGNII_ASD_res_darkcausality_sig[[3]]))
NeuNRGNII_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(NeuNRGNII_Normal_res_darkcausality_sig[[1]] + NeuNRGNII_Normal_res_darkcausality_sig[[2]] + NeuNRGNII_Normal_res_darkcausality_sig[[3]]))

# L56_ASD, L56_Normal
L56_ASD_res_darkcausality[[2]] <- L56_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L56_ASD_res_darkcausality[[3]] <- L56_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L56_ASD_res_darkcausality[[4]] <- L56_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L56_ASD_res_darkcausality_sig <- sigPC(L56_ASD_res_darkcausality[[2]], L56_ASD_res_darkcausality[[3]], L56_ASD_res_darkcausality[[4]])
L56_Normal_res_darkcausality[[2]] <- L56_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L56_Normal_res_darkcausality[[3]] <- L56_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L56_Normal_res_darkcausality[[4]] <- L56_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L56_Normal_res_darkcausality_sig <- sigPC(L56_Normal_res_darkcausality[[2]], L56_Normal_res_darkcausality[[3]], L56_Normal_res_darkcausality[[4]])
L56_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L56_ASD_res_darkcausality_sig[[1]] + L56_ASD_res_darkcausality_sig[[2]] + L56_ASD_res_darkcausality_sig[[3]]))
L56_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L56_Normal_res_darkcausality_sig[[1]] + L56_Normal_res_darkcausality_sig[[2]] + L56_Normal_res_darkcausality_sig[[3]]))

# Oligodendrocytes_ASD, Oligodendrocytes_Normal
Oligodendrocytes_ASD_res_darkcausality[[2]] <- Oligodendrocytes_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Oligodendrocytes_ASD_res_darkcausality[[3]] <- Oligodendrocytes_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Oligodendrocytes_ASD_res_darkcausality[[4]] <- Oligodendrocytes_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Oligodendrocytes_ASD_res_darkcausality_sig <- sigPC(Oligodendrocytes_ASD_res_darkcausality[[2]], Oligodendrocytes_ASD_res_darkcausality[[3]], Oligodendrocytes_ASD_res_darkcausality[[4]])
Oligodendrocytes_Normal_res_darkcausality[[2]] <- Oligodendrocytes_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Oligodendrocytes_Normal_res_darkcausality[[3]] <- Oligodendrocytes_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Oligodendrocytes_Normal_res_darkcausality[[4]] <- Oligodendrocytes_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Oligodendrocytes_Normal_res_darkcausality_sig <- sigPC(Oligodendrocytes_Normal_res_darkcausality[[2]], Oligodendrocytes_Normal_res_darkcausality[[3]], Oligodendrocytes_Normal_res_darkcausality[[4]])
Oligodendrocytes_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Oligodendrocytes_ASD_res_darkcausality_sig[[1]] + Oligodendrocytes_ASD_res_darkcausality_sig[[2]] + Oligodendrocytes_ASD_res_darkcausality_sig[[3]]))
Oligodendrocytes_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Oligodendrocytes_Normal_res_darkcausality_sig[[1]] + Oligodendrocytes_Normal_res_darkcausality_sig[[2]] + Oligodendrocytes_Normal_res_darkcausality_sig[[3]]))

# OPC_ASD, OPC_Normal
OPC_ASD_res_darkcausality[[2]] <- OPC_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
OPC_ASD_res_darkcausality[[3]] <- OPC_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
OPC_ASD_res_darkcausality[[4]] <- OPC_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
OPC_ASD_res_darkcausality_sig <- sigPC(OPC_ASD_res_darkcausality[[2]], OPC_ASD_res_darkcausality[[3]], OPC_ASD_res_darkcausality[[4]])
OPC_Normal_res_darkcausality[[2]] <- OPC_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
OPC_Normal_res_darkcausality[[3]] <- OPC_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
OPC_Normal_res_darkcausality[[4]] <- OPC_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
OPC_Normal_res_darkcausality_sig <- sigPC(OPC_Normal_res_darkcausality[[2]], OPC_Normal_res_darkcausality[[3]], OPC_Normal_res_darkcausality[[4]])
OPC_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(OPC_ASD_res_darkcausality_sig[[1]] + OPC_ASD_res_darkcausality_sig[[2]] + OPC_ASD_res_darkcausality_sig[[3]]))
OPC_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(OPC_Normal_res_darkcausality_sig[[1]] + OPC_Normal_res_darkcausality_sig[[2]] + OPC_Normal_res_darkcausality_sig[[3]]))

# ASTFB_ASD, ASTFB_Normal
ASTFB_ASD_res_darkcausality[[2]] <- ASTFB_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
ASTFB_ASD_res_darkcausality[[3]] <- ASTFB_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
ASTFB_ASD_res_darkcausality[[4]] <- ASTFB_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
ASTFB_ASD_res_darkcausality_sig <- sigPC(ASTFB_ASD_res_darkcausality[[2]], ASTFB_ASD_res_darkcausality[[3]], ASTFB_ASD_res_darkcausality[[4]])
ASTFB_Normal_res_darkcausality[[2]] <- ASTFB_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
ASTFB_Normal_res_darkcausality[[3]] <- ASTFB_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
ASTFB_Normal_res_darkcausality[[4]] <- ASTFB_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
ASTFB_Normal_res_darkcausality_sig <- sigPC(ASTFB_Normal_res_darkcausality[[2]], ASTFB_Normal_res_darkcausality[[3]], ASTFB_Normal_res_darkcausality[[4]])
ASTFB_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(ASTFB_ASD_res_darkcausality_sig[[1]] + ASTFB_ASD_res_darkcausality_sig[[2]] + ASTFB_ASD_res_darkcausality_sig[[3]]))
ASTFB_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(ASTFB_Normal_res_darkcausality_sig[[1]] + ASTFB_Normal_res_darkcausality_sig[[2]] + ASTFB_Normal_res_darkcausality_sig[[3]]))

# Endothelial_ASD, Endothelial_Normal
Endothelial_ASD_res_darkcausality[[2]] <- Endothelial_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Endothelial_ASD_res_darkcausality[[3]] <- Endothelial_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Endothelial_ASD_res_darkcausality[[4]] <- Endothelial_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Endothelial_ASD_res_darkcausality_sig <- sigPC(Endothelial_ASD_res_darkcausality[[2]], Endothelial_ASD_res_darkcausality[[3]], Endothelial_ASD_res_darkcausality[[4]])
Endothelial_Normal_res_darkcausality[[2]] <- Endothelial_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Endothelial_Normal_res_darkcausality[[3]] <- Endothelial_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Endothelial_Normal_res_darkcausality[[4]] <- Endothelial_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Endothelial_Normal_res_darkcausality_sig <- sigPC(Endothelial_Normal_res_darkcausality[[2]], Endothelial_Normal_res_darkcausality[[3]], Endothelial_Normal_res_darkcausality[[4]])
Endothelial_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Endothelial_ASD_res_darkcausality_sig[[1]] + Endothelial_ASD_res_darkcausality_sig[[2]] + Endothelial_ASD_res_darkcausality_sig[[3]]))
Endothelial_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Endothelial_Normal_res_darkcausality_sig[[1]] + Endothelial_Normal_res_darkcausality_sig[[2]] + Endothelial_Normal_res_darkcausality_sig[[3]]))

# Microglia_ASD, Microglia_Normal
Microglia_ASD_res_darkcausality[[2]] <- Microglia_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Microglia_ASD_res_darkcausality[[3]] <- Microglia_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Microglia_ASD_res_darkcausality[[4]] <- Microglia_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Microglia_ASD_res_darkcausality_sig <- sigPC(Microglia_ASD_res_darkcausality[[2]], Microglia_ASD_res_darkcausality[[3]], Microglia_ASD_res_darkcausality[[4]])
Microglia_Normal_res_darkcausality[[2]] <- Microglia_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Microglia_Normal_res_darkcausality[[3]] <- Microglia_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Microglia_Normal_res_darkcausality[[4]] <- Microglia_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Microglia_Normal_res_darkcausality_sig <- sigPC(Microglia_Normal_res_darkcausality[[2]], Microglia_Normal_res_darkcausality[[3]], Microglia_Normal_res_darkcausality[[4]])
Microglia_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Microglia_ASD_res_darkcausality_sig[[1]] + Microglia_ASD_res_darkcausality_sig[[2]] + Microglia_ASD_res_darkcausality_sig[[3]]))
Microglia_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Microglia_Normal_res_darkcausality_sig[[1]] + Microglia_Normal_res_darkcausality_sig[[2]] + Microglia_Normal_res_darkcausality_sig[[3]]))

# NeuNRGNI_ASD, NeuNRGNI_Normal
NeuNRGNI_ASD_res_darkcausality[[2]] <- NeuNRGNI_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
NeuNRGNI_ASD_res_darkcausality[[3]] <- NeuNRGNI_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
NeuNRGNI_ASD_res_darkcausality[[4]] <- NeuNRGNI_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
NeuNRGNI_ASD_res_darkcausality_sig <- sigPC(NeuNRGNI_ASD_res_darkcausality[[2]], NeuNRGNI_ASD_res_darkcausality[[3]], NeuNRGNI_ASD_res_darkcausality[[4]])
NeuNRGNI_Normal_res_darkcausality[[2]] <- NeuNRGNI_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
NeuNRGNI_Normal_res_darkcausality[[3]] <- NeuNRGNI_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
NeuNRGNI_Normal_res_darkcausality[[4]] <- NeuNRGNI_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
NeuNRGNI_Normal_res_darkcausality_sig <- sigPC(NeuNRGNI_Normal_res_darkcausality[[2]], NeuNRGNI_Normal_res_darkcausality[[3]], NeuNRGNI_Normal_res_darkcausality[[4]])
NeuNRGNI_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(NeuNRGNI_ASD_res_darkcausality_sig[[1]] + NeuNRGNI_ASD_res_darkcausality_sig[[2]] + NeuNRGNI_ASD_res_darkcausality_sig[[3]]))
NeuNRGNI_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(NeuNRGNI_Normal_res_darkcausality_sig[[1]] + NeuNRGNI_Normal_res_darkcausality_sig[[2]] + NeuNRGNI_Normal_res_darkcausality_sig[[3]]))

# INVIP_ASD, INVIP_Normal
INVIP_ASD_res_darkcausality[[2]] <- INVIP_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INVIP_ASD_res_darkcausality[[3]] <- INVIP_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INVIP_ASD_res_darkcausality[[4]] <- INVIP_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INVIP_ASD_res_darkcausality_sig <- sigPC(INVIP_ASD_res_darkcausality[[2]], INVIP_ASD_res_darkcausality[[3]], INVIP_ASD_res_darkcausality[[4]])
INVIP_Normal_res_darkcausality[[2]] <- INVIP_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INVIP_Normal_res_darkcausality[[3]] <- INVIP_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INVIP_Normal_res_darkcausality[[4]] <- INVIP_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INVIP_Normal_res_darkcausality_sig <- sigPC(INVIP_Normal_res_darkcausality[[2]], INVIP_Normal_res_darkcausality[[3]], INVIP_Normal_res_darkcausality[[4]])
INVIP_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INVIP_ASD_res_darkcausality_sig[[1]] + INVIP_ASD_res_darkcausality_sig[[2]] + INVIP_ASD_res_darkcausality_sig[[3]]))
INVIP_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INVIP_Normal_res_darkcausality_sig[[1]] + INVIP_Normal_res_darkcausality_sig[[2]] + INVIP_Normal_res_darkcausality_sig[[3]]))

# L56CC_ASD, L56CC_Normal
L56CC_ASD_res_darkcausality[[2]] <- L56CC_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L56CC_ASD_res_darkcausality[[3]] <- L56CC_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L56CC_ASD_res_darkcausality[[4]] <- L56CC_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L56CC_ASD_res_darkcausality_sig <- sigPC(L56CC_ASD_res_darkcausality[[2]], L56CC_ASD_res_darkcausality[[3]], L56CC_ASD_res_darkcausality[[4]])
L56CC_Normal_res_darkcausality[[2]] <- L56CC_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L56CC_Normal_res_darkcausality[[3]] <- L56CC_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L56CC_Normal_res_darkcausality[[4]] <- L56CC_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L56CC_Normal_res_darkcausality_sig <- sigPC(L56CC_Normal_res_darkcausality[[2]], L56CC_Normal_res_darkcausality[[3]], L56CC_Normal_res_darkcausality[[4]])
L56CC_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L56CC_ASD_res_darkcausality_sig[[1]] + L56CC_ASD_res_darkcausality_sig[[2]] + L56CC_ASD_res_darkcausality_sig[[3]]))
L56CC_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L56CC_Normal_res_darkcausality_sig[[1]] + L56CC_Normal_res_darkcausality_sig[[2]] + L56CC_Normal_res_darkcausality_sig[[3]]))

# INSV2C_ASD, INSV2C_Normal
INSV2C_ASD_res_darkcausality[[2]] <- INSV2C_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INSV2C_ASD_res_darkcausality[[3]] <- INSV2C_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INSV2C_ASD_res_darkcausality[[4]] <- INSV2C_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INSV2C_ASD_res_darkcausality_sig <- sigPC(INSV2C_ASD_res_darkcausality[[2]], INSV2C_ASD_res_darkcausality[[3]], INSV2C_ASD_res_darkcausality[[4]])
INSV2C_Normal_res_darkcausality[[2]] <- INSV2C_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INSV2C_Normal_res_darkcausality[[3]] <- INSV2C_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INSV2C_Normal_res_darkcausality[[4]] <- INSV2C_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INSV2C_Normal_res_darkcausality_sig <- sigPC(INSV2C_Normal_res_darkcausality[[2]], INSV2C_Normal_res_darkcausality[[3]], INSV2C_Normal_res_darkcausality[[4]])
INSV2C_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INSV2C_ASD_res_darkcausality_sig[[1]] + INSV2C_ASD_res_darkcausality_sig[[2]] + INSV2C_ASD_res_darkcausality_sig[[3]]))
INSV2C_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INSV2C_Normal_res_darkcausality_sig[[1]] + INSV2C_Normal_res_darkcausality_sig[[2]] + INSV2C_Normal_res_darkcausality_sig[[3]]))

# L23_ASD, L23_Normal
L23_ASD_res_darkcausality[[2]] <- L23_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L23_ASD_res_darkcausality[[3]] <- L23_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L23_ASD_res_darkcausality[[4]] <- L23_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L23_ASD_res_darkcausality_sig <- sigPC(L23_ASD_res_darkcausality[[2]], L23_ASD_res_darkcausality[[3]], L23_ASD_res_darkcausality[[4]])
L23_Normal_res_darkcausality[[2]] <- L23_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L23_Normal_res_darkcausality[[3]] <- L23_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L23_Normal_res_darkcausality[[4]] <- L23_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L23_Normal_res_darkcausality_sig <- sigPC(L23_Normal_res_darkcausality[[2]], L23_Normal_res_darkcausality[[3]], L23_Normal_res_darkcausality[[4]])
L23_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L23_ASD_res_darkcausality_sig[[1]] + L23_ASD_res_darkcausality_sig[[2]] + L23_ASD_res_darkcausality_sig[[3]]))
L23_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L23_Normal_res_darkcausality_sig[[1]] + L23_Normal_res_darkcausality_sig[[2]] + L23_Normal_res_darkcausality_sig[[3]]))

# INPV_ASD, INPV_Normal
INPV_ASD_res_darkcausality[[2]] <- INPV_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INPV_ASD_res_darkcausality[[3]] <- INPV_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INPV_ASD_res_darkcausality[[4]] <- INPV_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INPV_ASD_res_darkcausality_sig <- sigPC(INPV_ASD_res_darkcausality[[2]], INPV_ASD_res_darkcausality[[3]], INPV_ASD_res_darkcausality[[4]])
INPV_Normal_res_darkcausality[[2]] <- INPV_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INPV_Normal_res_darkcausality[[3]] <- INPV_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INPV_Normal_res_darkcausality[[4]] <- INPV_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INPV_Normal_res_darkcausality_sig <- sigPC(INPV_Normal_res_darkcausality[[2]], INPV_Normal_res_darkcausality[[3]], INPV_Normal_res_darkcausality[[4]])
INPV_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INPV_ASD_res_darkcausality_sig[[1]] + INPV_ASD_res_darkcausality_sig[[2]] + INPV_ASD_res_darkcausality_sig[[3]]))
INPV_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INPV_Normal_res_darkcausality_sig[[1]] + INPV_Normal_res_darkcausality_sig[[2]] + INPV_Normal_res_darkcausality_sig[[3]]))

# L4_ASD, L4_Normal
L4_ASD_res_darkcausality[[2]] <- L4_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L4_ASD_res_darkcausality[[3]] <- L4_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L4_ASD_res_darkcausality[[4]] <- L4_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L4_ASD_res_darkcausality_sig <- sigPC(L4_ASD_res_darkcausality[[2]], L4_ASD_res_darkcausality[[3]], L4_ASD_res_darkcausality[[4]])
L4_Normal_res_darkcausality[[2]] <- L4_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
L4_Normal_res_darkcausality[[3]] <- L4_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
L4_Normal_res_darkcausality[[4]] <- L4_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
L4_Normal_res_darkcausality_sig <- sigPC(L4_Normal_res_darkcausality[[2]], L4_Normal_res_darkcausality[[3]], L4_Normal_res_darkcausality[[4]])
L4_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L4_ASD_res_darkcausality_sig[[1]] + L4_ASD_res_darkcausality_sig[[2]] + L4_ASD_res_darkcausality_sig[[3]]))
L4_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(L4_Normal_res_darkcausality_sig[[1]] + L4_Normal_res_darkcausality_sig[[2]] + L4_Normal_res_darkcausality_sig[[3]]))

# INSST_ASD, INSST_Normal
INSST_ASD_res_darkcausality[[2]] <- INSST_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INSST_ASD_res_darkcausality[[3]] <- INSST_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INSST_ASD_res_darkcausality[[4]] <- INSST_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INSST_ASD_res_darkcausality_sig <- sigPC(INSST_ASD_res_darkcausality[[2]], INSST_ASD_res_darkcausality[[3]], INSST_ASD_res_darkcausality[[4]])
INSST_Normal_res_darkcausality[[2]] <- INSST_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
INSST_Normal_res_darkcausality[[3]] <- INSST_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
INSST_Normal_res_darkcausality[[4]] <- INSST_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
INSST_Normal_res_darkcausality_sig <- sigPC(INSST_Normal_res_darkcausality[[2]], INSST_Normal_res_darkcausality[[3]], INSST_Normal_res_darkcausality[[4]])
INSST_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INSST_ASD_res_darkcausality_sig[[1]] + INSST_ASD_res_darkcausality_sig[[2]] + INSST_ASD_res_darkcausality_sig[[3]]))
INSST_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(INSST_Normal_res_darkcausality_sig[[1]] + INSST_Normal_res_darkcausality_sig[[2]] + INSST_Normal_res_darkcausality_sig[[3]]))

# Neumat_ASD, Neumat_Normal
Neumat_ASD_res_darkcausality[[2]] <- Neumat_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Neumat_ASD_res_darkcausality[[3]] <- Neumat_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Neumat_ASD_res_darkcausality[[4]] <- Neumat_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Neumat_ASD_res_darkcausality_sig <- sigPC(Neumat_ASD_res_darkcausality[[2]], Neumat_ASD_res_darkcausality[[3]], Neumat_ASD_res_darkcausality[[4]])
Neumat_Normal_res_darkcausality[[2]] <- Neumat_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
Neumat_Normal_res_darkcausality[[3]] <- Neumat_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
Neumat_Normal_res_darkcausality[[4]] <- Neumat_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
Neumat_Normal_res_darkcausality_sig <- sigPC(Neumat_Normal_res_darkcausality[[2]], Neumat_Normal_res_darkcausality[[3]], Neumat_Normal_res_darkcausality[[4]])
Neumat_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Neumat_ASD_res_darkcausality_sig[[1]] + Neumat_ASD_res_darkcausality_sig[[2]] + Neumat_ASD_res_darkcausality_sig[[3]]))
Neumat_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(Neumat_Normal_res_darkcausality_sig[[1]] + Neumat_Normal_res_darkcausality_sig[[2]] + Neumat_Normal_res_darkcausality_sig[[3]]))

# ASTPP_ASD, ASTPP_Normal
ASTPP_ASD_res_darkcausality[[2]] <- ASTPP_ASD_res_darkcausality[[2]] %>% replace(is.na(.), 0)
ASTPP_ASD_res_darkcausality[[3]] <- ASTPP_ASD_res_darkcausality[[3]] %>% replace(is.na(.), 0)
ASTPP_ASD_res_darkcausality[[4]] <- ASTPP_ASD_res_darkcausality[[4]] %>% replace(is.na(.), 0)
ASTPP_ASD_res_darkcausality_sig <- sigPC(ASTPP_ASD_res_darkcausality[[2]], ASTPP_ASD_res_darkcausality[[3]], ASTPP_ASD_res_darkcausality[[4]])
ASTPP_Normal_res_darkcausality[[2]] <- ASTPP_Normal_res_darkcausality[[2]] %>% replace(is.na(.), 0)
ASTPP_Normal_res_darkcausality[[3]] <- ASTPP_Normal_res_darkcausality[[3]] %>% replace(is.na(.), 0)
ASTPP_Normal_res_darkcausality[[4]] <- ASTPP_Normal_res_darkcausality[[4]] %>% replace(is.na(.), 0)
ASTPP_Normal_res_darkcausality_sig <- sigPC(ASTPP_Normal_res_darkcausality[[2]], ASTPP_Normal_res_darkcausality[[3]], ASTPP_Normal_res_darkcausality[[4]])
ASTPP_ASD_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(ASTPP_ASD_res_darkcausality_sig[[1]] + ASTPP_ASD_res_darkcausality_sig[[2]] + ASTPP_ASD_res_darkcausality_sig[[3]]))
ASTPP_Normal_res_darkcausality_graph <- graph_from_biadjacency_matrix(t(ASTPP_Normal_res_darkcausality_sig[[1]] + ASTPP_Normal_res_darkcausality_sig[[2]] + ASTPP_Normal_res_darkcausality_sig[[3]]))

################################################################################  Granger application ############################################################################
# ASD and Normal
ASD_res_granger <- granger_parallel(ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Normal_res_granger <- granger_parallel(Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(ASD_res_granger))
Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Normal_res_granger))

# ACC_ASD, ACC_Normal, PFC_ASD, PFC_Normal
ACC_ASD_res_granger <- granger_parallel(ACC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
ACC_Normal_res_granger <- granger_parallel(ACC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
PFC_ASD_res_granger <- granger_parallel(PFC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
PFC_Normal_res_granger <- granger_parallel(PFC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
ACC_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(ACC_ASD_res_granger))
ACC_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(ACC_Normal_res_granger))
PFC_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(PFC_ASD_res_granger))
PFC_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(PFC_Normal_res_granger))

# Lower18_ASD, Lower18_Normal, Larger18_ASD, Larger18_Normal
Lower18_ASD_res_granger <- granger_parallel(Lower18_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Lower18_Normal_res_granger <- granger_parallel(Lower18_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Larger18_ASD_res_granger <- granger_parallel(Larger18_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Larger18_Normal_res_granger <- granger_parallel(Larger18_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Lower18_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Lower18_ASD_res_granger))
Lower18_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Lower18_Normal_res_granger))
Larger18_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Larger18_ASD_res_granger))
Larger18_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Larger18_Normal_res_granger))

# Male_ASD, Male_Normal, Female_ASD, Female_Normal
Male_ASD_res_granger <- granger_parallel(Male_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Male_Normal_res_granger <- granger_parallel(Male_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Female_ASD_res_granger <- granger_parallel(Female_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Female_Normal_res_granger <- granger_parallel(Female_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Male_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Male_ASD_res_granger))
Male_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Male_Normal_res_granger))
Female_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Female_ASD_res_granger))
Female_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Female_Normal_res_granger))

# NeuNRGNII_ASD, NeuNRGNII_Normal
NeuNRGNII_ASD_res_granger <- granger_parallel(NeuNRGNII_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
NeuNRGNII_Normal_res_granger <- granger_parallel(NeuNRGNII_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
NeuNRGNII_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(NeuNRGNII_ASD_res_granger))
NeuNRGNII_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(NeuNRGNII_Normal_res_granger))

# L56_ASD, L56_Normal
L56_ASD_res_granger <- granger_parallel(L56_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L56_Normal_res_granger <- granger_parallel(L56_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L56_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(L56_ASD_res_granger))
L56_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(L56_Normal_res_granger))

# Oligodendrocytes_ASD, Oligodendrocytes_Normal
Oligodendrocytes_ASD_res_granger <- granger_parallel(Oligodendrocytes_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Oligodendrocytes_Normal_res_granger <- granger_parallel(Oligodendrocytes_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Oligodendrocytes_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Oligodendrocytes_ASD_res_granger))
Oligodendrocytes_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Oligodendrocytes_Normal_res_granger))

# OPC_ASD, OPC_Normal
OPC_ASD_res_granger <- granger_parallel(OPC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
OPC_Normal_res_granger <- granger_parallel(OPC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
OPC_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(OPC_ASD_res_granger))
OPC_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(OPC_Normal_res_granger))

# ASTFB_ASD, ASTFB_Normal
ASTFB_ASD_res_granger <- granger_parallel(ASTFB_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
ASTFB_Normal_res_granger <- granger_parallel(ASTFB_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
ASTFB_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(ASTFB_ASD_res_granger))
ASTFB_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(ASTFB_Normal_res_granger))

# Endothelial_ASD, Endothelial_Normal
Endothelial_ASD_res_granger <- granger_parallel(Endothelial_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Endothelial_Normal_res_granger <- granger_parallel(Endothelial_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Endothelial_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Endothelial_ASD_res_granger))
Endothelial_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Endothelial_Normal_res_granger))

# Microglia_ASD, Microglia_Normal
Microglia_ASD_res_granger <- granger_parallel(Microglia_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Microglia_Normal_res_granger <- granger_parallel(Microglia_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Microglia_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Microglia_ASD_res_granger))
Microglia_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Microglia_Normal_res_granger))

# NeuNRGNI_ASD, NeuNRGNI_Normal
NeuNRGNI_ASD_res_granger <- granger_parallel(NeuNRGNI_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
NeuNRGNI_Normal_res_granger <- granger_parallel(NeuNRGNI_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
NeuNRGNI_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(NeuNRGNI_ASD_res_granger))
NeuNRGNI_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(NeuNRGNI_Normal_res_granger))

# INVIP_ASD, INVIP_Normal
INVIP_ASD_res_granger <- granger_parallel(INVIP_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INVIP_Normal_res_granger <- granger_parallel(INVIP_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INVIP_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(INVIP_ASD_res_granger))
INVIP_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(INVIP_Normal_res_granger))

# L56CC_ASD, L56CC_Normal
L56CC_ASD_res_granger <- granger_parallel(L56CC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L56CC_Normal_res_granger <- granger_parallel(L56CC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L56CC_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(L56CC_ASD_res_granger))
L56CC_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(L56CC_Normal_res_granger))

# INSV2C_ASD, INSV2C_Normal
INSV2C_ASD_res_granger <- granger_parallel(INSV2C_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INSV2C_Normal_res_granger <- granger_parallel(INSV2C_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INSV2C_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(INSV2C_ASD_res_granger))
INSV2C_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(INSV2C_Normal_res_granger))

# L23_ASD, L23_Normal
L23_ASD_res_granger <- granger_parallel(L23_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L23_Normal_res_granger <- granger_parallel(L23_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L23_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(L23_ASD_res_granger))
L23_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(L23_Normal_res_granger))

# INPV_ASD, INPV_Normal
INPV_ASD_res_granger <- granger_parallel(INPV_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INPV_Normal_res_granger <- granger_parallel(INPV_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INPV_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(INPV_ASD_res_granger))
INPV_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(INPV_Normal_res_granger))

# L4_ASD, L4_Normal
L4_ASD_res_granger <- granger_parallel(L4_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L4_Normal_res_granger <- granger_parallel(L4_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
L4_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(L4_ASD_res_granger))
L4_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(L4_Normal_res_granger))

# INSST_ASD, INSST_Normal
INSST_ASD_res_granger <- granger_parallel(INSST_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INSST_Normal_res_granger <- granger_parallel(INSST_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
INSST_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(INSST_ASD_res_granger))
INSST_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(INSST_Normal_res_granger))

# Neumat_ASD, Neumat_Normal
Neumat_ASD_res_granger <- granger_parallel(Neumat_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Neumat_Normal_res_granger <- granger_parallel(Neumat_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
Neumat_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(Neumat_ASD_res_granger))
Neumat_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(Neumat_Normal_res_granger))

# ASTPP_ASD, ASTPP_Normal
ASTPP_ASD_res_granger <- granger_parallel(ASTPP_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
ASTPP_Normal_res_granger <- granger_parallel(ASTPP_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, num.cores = 48)
ASTPP_ASD_res_granger_graph <- graph_from_biadjacency_matrix(t(ASTPP_ASD_res_granger))
ASTPP_Normal_res_granger_graph <- graph_from_biadjacency_matrix(t(ASTPP_Normal_res_granger))

###########################################################################  seqICP application ################################################################################################
# ASD and Normal
ASD_res_seqICP <- seqICP_parallel(ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Normal_res_seqICP <- seqICP_parallel(Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(ASD_res_seqICP))
Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Normal_res_seqICP))

# ACC_ASD, ACC_Normal, PFC_ASD, PFC_Normal
ACC_ASD_res_seqICP <- seqICP_parallel(ACC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
ACC_Normal_res_seqICP <- seqICP_parallel(ACC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
PFC_ASD_res_seqICP <- seqICP_parallel(PFC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
PFC_Normal_res_seqICP <- seqICP_parallel(PFC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
ACC_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(ACC_ASD_res_seqICP))
ACC_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(ACC_Normal_res_seqICP))
PFC_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(PFC_ASD_res_seqICP))
PFC_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(PFC_Normal_res_seqICP))

# Lower18_ASD, Lower18_Normal, Larger18_ASD, Larger18_Normal
Lower18_ASD_res_seqICP <- seqICP_parallel(Lower18_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Lower18_Normal_res_seqICP <- seqICP_parallel(Lower18_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Larger18_ASD_res_seqICP <- seqICP_parallel(Larger18_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Larger18_Normal_res_seqICP <- seqICP_parallel(Larger18_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Lower18_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Lower18_ASD_res_seqICP))
Lower18_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Lower18_Normal_res_seqICP))
Larger18_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Larger18_ASD_res_seqICP))
Larger18_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Larger18_Normal_res_seqICP))

# Male_ASD, Male_Normal, Female_ASD, Female_Normal
Male_ASD_res_seqICP <- seqICP_parallel(Male_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Male_Normal_res_seqICP <- seqICP_parallel(Male_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Female_ASD_res_seqICP <- seqICP_parallel(Female_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Female_Normal_res_seqICP <- seqICP_parallel(Female_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Male_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Male_ASD_res_seqICP))
Male_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Male_Normal_res_seqICP))
Female_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Female_ASD_res_seqICP))
Female_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Female_Normal_res_seqICP))

# NeuNRGNII_ASD, NeuNRGNII_Normal
NeuNRGNII_ASD_res_seqICP <- seqICP_parallel(NeuNRGNII_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
NeuNRGNII_Normal_res_seqICP <- seqICP_parallel(NeuNRGNII_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
NeuNRGNII_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(NeuNRGNII_ASD_res_seqICP))
NeuNRGNII_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(NeuNRGNII_Normal_res_seqICP))

# L56_ASD, L56_Normal
L56_ASD_res_seqICP <- seqICP_parallel(L56_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L56_Normal_res_seqICP <- seqICP_parallel(L56_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L56_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L56_ASD_res_seqICP))
L56_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L56_Normal_res_seqICP))

# Oligodendrocytes_ASD, Oligodendrocytes_Normal
Oligodendrocytes_ASD_res_seqICP <- seqICP_parallel(Oligodendrocytes_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Oligodendrocytes_Normal_res_seqICP <- seqICP_parallel(Oligodendrocytes_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Oligodendrocytes_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Oligodendrocytes_ASD_res_seqICP))
Oligodendrocytes_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Oligodendrocytes_Normal_res_seqICP))

# OPC_ASD, OPC_Normal
OPC_ASD_res_seqICP <- seqICP_parallel(OPC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
OPC_Normal_res_seqICP <- seqICP_parallel(OPC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
OPC_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(OPC_ASD_res_seqICP))
OPC_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(OPC_Normal_res_seqICP))

# ASTFB_ASD, ASTFB_Normal
ASTFB_ASD_res_seqICP <- seqICP_parallel(ASTFB_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
ASTFB_Normal_res_seqICP <- seqICP_parallel(ASTFB_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
ASTFB_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(ASTFB_ASD_res_seqICP))
ASTFB_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(ASTFB_Normal_res_seqICP))

# Endothelial_ASD, Endothelial_Normal
Endothelial_ASD_res_seqICP <- seqICP_parallel(Endothelial_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Endothelial_Normal_res_seqICP <- seqICP_parallel(Endothelial_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Endothelial_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Endothelial_ASD_res_seqICP))
Endothelial_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Endothelial_Normal_res_seqICP))

# Microglia_ASD, Microglia_Normal
Microglia_ASD_res_seqICP <- seqICP_parallel(Microglia_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Microglia_Normal_res_seqICP <- seqICP_parallel(Microglia_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Microglia_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Microglia_ASD_res_seqICP))
Microglia_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Microglia_Normal_res_seqICP))

# NeuNRGNI_ASD, NeuNRGNI_Normal
NeuNRGNI_ASD_res_seqICP <- seqICP_parallel(NeuNRGNI_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
NeuNRGNI_Normal_res_seqICP <- seqICP_parallel(NeuNRGNI_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
NeuNRGNI_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(NeuNRGNI_ASD_res_seqICP))
NeuNRGNI_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(NeuNRGNI_Normal_res_seqICP))

# INVIP_ASD, INVIP_Normal
INVIP_ASD_res_seqICP <- seqICP_parallel(INVIP_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INVIP_Normal_res_seqICP <- seqICP_parallel(INVIP_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INVIP_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INVIP_ASD_res_seqICP))
INVIP_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INVIP_Normal_res_seqICP))

# L56CC_ASD, L56CC_Normal
L56CC_ASD_res_seqICP <- seqICP_parallel(L56CC_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L56CC_Normal_res_seqICP <- seqICP_parallel(L56CC_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L56CC_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L56CC_ASD_res_seqICP))
L56CC_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L56CC_Normal_res_seqICP))

# INSV2C_ASD, INSV2C_Normal
INSV2C_ASD_res_seqICP <- seqICP_parallel(INSV2C_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INSV2C_Normal_res_seqICP <- seqICP_parallel(INSV2C_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INSV2C_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INSV2C_ASD_res_seqICP))
INSV2C_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INSV2C_Normal_res_seqICP))

# L23_ASD, L23_Normal
L23_ASD_res_seqICP <- seqICP_parallel(L23_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L23_Normal_res_seqICP <- seqICP_parallel(L23_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L23_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L23_ASD_res_seqICP))
L23_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L23_Normal_res_seqICP))

# INPV_ASD, INPV_Normal
INPV_ASD_res_seqICP <- seqICP_parallel(INPV_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INPV_Normal_res_seqICP <- seqICP_parallel(INPV_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INPV_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INPV_ASD_res_seqICP))
INPV_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INPV_Normal_res_seqICP))

# L4_ASD, L4_Normal
L4_ASD_res_seqICP <- seqICP_parallel(L4_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L4_Normal_res_seqICP <- seqICP_parallel(L4_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
L4_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L4_ASD_res_seqICP))
L4_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(L4_Normal_res_seqICP))

# INSST_ASD, INSST_Normal
INSST_ASD_res_seqICP <- seqICP_parallel(INSST_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INSST_Normal_res_seqICP <- seqICP_parallel(INSST_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
INSST_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INSST_ASD_res_seqICP))
INSST_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(INSST_Normal_res_seqICP))

# Neumat_ASD, Neumat_Normal
Neumat_ASD_res_seqICP <- seqICP_parallel(Neumat_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Neumat_Normal_res_seqICP <- seqICP_parallel(Neumat_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
Neumat_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Neumat_ASD_res_seqICP))
Neumat_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(Neumat_Normal_res_seqICP))

# ASTPP_ASD, ASTPP_Normal
ASTPP_ASD_res_seqICP <- seqICP_parallel(ASTPP_ASD_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
ASTPP_Normal_res_seqICP <- seqICP_parallel(ASTPP_Normal_model_data, cause, effect, pvalue_cutoff = 0.01, method = "seqICPnl", num.cores = 48)
ASTPP_ASD_res_seqICP_graph <- graph_from_biadjacency_matrix(t(ASTPP_ASD_res_seqICP))
ASTPP_Normal_res_seqICP_graph <- graph_from_biadjacency_matrix(t(ASTPP_Normal_res_seqICP))

############################################################################ Incorporating priori information of lncRNA targets #########################################################################
lncRTarget_priori <- as.matrix(read.csv("LncTar.csv", header = TRUE, sep=","))
lncRTarget_priori_graph <- make_graph(c(t(lncRTarget_priori[, 1:2])), directed = FALSE)

# Clean
ASD_res_darkcausality_priori_graph <- ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Normal_res_darkcausality_priori_graph <- Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

ACC_ASD_res_darkcausality_priori_graph <- ACC_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
ACC_Normal_res_darkcausality_priori_graph <- ACC_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph
PFC_ASD_res_darkcausality_priori_graph <- PFC_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
PFC_Normal_res_darkcausality_priori_graph <- PFC_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

Lower18_ASD_res_darkcausality_priori_graph <- Lower18_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Lower18_Normal_res_darkcausality_priori_graph <- Lower18_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph
Larger18_ASD_res_darkcausality_priori_graph <- Larger18_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Larger18_Normal_res_darkcausality_priori_graph <- Larger18_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

Male_ASD_res_darkcausality_priori_graph <- Male_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Male_Normal_res_darkcausality_priori_graph <- Male_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph
Female_ASD_res_darkcausality_priori_graph <- Female_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Female_Normal_res_darkcausality_priori_graph <- Female_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

NeuNRGNII_ASD_res_darkcausality_priori_graph <- NeuNRGNII_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
NeuNRGNII_Normal_res_darkcausality_priori_graph <- NeuNRGNII_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

L56_ASD_res_darkcausality_priori_graph <- L56_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
L56_Normal_res_darkcausality_priori_graph <- L56_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

Oligodendrocytes_ASD_res_darkcausality_priori_graph <- Oligodendrocytes_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Oligodendrocytes_Normal_res_darkcausality_priori_graph <- Oligodendrocytes_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

OPC_ASD_res_darkcausality_priori_graph <- OPC_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
OPC_Normal_res_darkcausality_priori_graph <- OPC_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

ASTFB_ASD_res_darkcausality_priori_graph <- ASTFB_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
ASTFB_Normal_res_darkcausality_priori_graph <- ASTFB_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

Endothelial_ASD_res_darkcausality_priori_graph <- Endothelial_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Endothelial_Normal_res_darkcausality_priori_graph <- Endothelial_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

Microglia_ASD_res_darkcausality_priori_graph <- Microglia_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Microglia_Normal_res_darkcausality_priori_graph <- Microglia_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

NeuNRGNI_ASD_res_darkcausality_priori_graph <- NeuNRGNI_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
NeuNRGNI_Normal_res_darkcausality_priori_graph <- NeuNRGNI_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

INVIP_ASD_res_darkcausality_priori_graph <- INVIP_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
INVIP_Normal_res_darkcausality_priori_graph <- INVIP_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

L56CC_ASD_res_darkcausality_priori_graph <- L56CC_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
L56CC_Normal_res_darkcausality_priori_graph <- L56CC_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

INSV2C_ASD_res_darkcausality_priori_graph <- INSV2C_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
INSV2C_Normal_res_darkcausality_priori_graph <- INSV2C_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

L23_ASD_res_darkcausality_priori_graph <- L23_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
L23_Normal_res_darkcausality_priori_graph <- L23_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

INPV_ASD_res_darkcausality_priori_graph <- INPV_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
INPV_Normal_res_darkcausality_priori_graph <- INPV_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

L4_ASD_res_darkcausality_priori_graph <- L4_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
L4_Normal_res_darkcausality_priori_graph <- L4_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

INSST_ASD_res_darkcausality_priori_graph <- INSST_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
INSST_Normal_res_darkcausality_priori_graph <- INSST_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

Neumat_ASD_res_darkcausality_priori_graph <- Neumat_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
Neumat_Normal_res_darkcausality_priori_graph <- Neumat_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

ASTPP_ASD_res_darkcausality_priori_graph <- ASTPP_ASD_res_darkcausality_graph %s% lncRTarget_priori_graph
ASTPP_Normal_res_darkcausality_priori_graph <- ASTPP_Normal_res_darkcausality_graph %s% lncRTarget_priori_graph

# Granger
ASD_res_granger_priori_graph <- ASD_res_granger_graph %s% lncRTarget_priori_graph
Normal_res_granger_priori_graph <- Normal_res_granger_graph %s% lncRTarget_priori_graph

ACC_ASD_res_granger_priori_graph <- ACC_ASD_res_granger_graph %s% lncRTarget_priori_graph
ACC_Normal_res_granger_priori_graph <- ACC_Normal_res_granger_graph %s% lncRTarget_priori_graph
PFC_ASD_res_granger_priori_graph <- PFC_ASD_res_granger_graph %s% lncRTarget_priori_graph
PFC_Normal_res_granger_priori_graph <- PFC_Normal_res_granger_graph %s% lncRTarget_priori_graph

Lower18_ASD_res_granger_priori_graph <- Lower18_ASD_res_granger_graph %s% lncRTarget_priori_graph
Lower18_Normal_res_granger_priori_graph <- Lower18_Normal_res_granger_graph %s% lncRTarget_priori_graph
Larger18_ASD_res_granger_priori_graph <- Larger18_ASD_res_granger_graph %s% lncRTarget_priori_graph
Larger18_Normal_res_granger_priori_graph <- Larger18_Normal_res_granger_graph %s% lncRTarget_priori_graph

Male_ASD_res_granger_priori_graph <- Male_ASD_res_granger_graph %s% lncRTarget_priori_graph
Male_Normal_res_granger_priori_graph <- Male_Normal_res_granger_graph %s% lncRTarget_priori_graph
Female_ASD_res_granger_priori_graph <- Female_ASD_res_granger_graph %s% lncRTarget_priori_graph
Female_Normal_res_granger_priori_graph <- Female_Normal_res_granger_graph %s% lncRTarget_priori_graph

NeuNRGNII_ASD_res_granger_priori_graph <- NeuNRGNII_ASD_res_granger_graph %s% lncRTarget_priori_graph
NeuNRGNII_Normal_res_granger_priori_graph <- NeuNRGNII_Normal_res_granger_graph %s% lncRTarget_priori_graph

L56_ASD_res_granger_priori_graph <- L56_ASD_res_granger_graph %s% lncRTarget_priori_graph
L56_Normal_res_granger_priori_graph <- L56_Normal_res_granger_graph %s% lncRTarget_priori_graph

Oligodendrocytes_ASD_res_granger_priori_graph <- Oligodendrocytes_ASD_res_granger_graph %s% lncRTarget_priori_graph
Oligodendrocytes_Normal_res_granger_priori_graph <- Oligodendrocytes_Normal_res_granger_graph %s% lncRTarget_priori_graph

OPC_ASD_res_granger_priori_graph <- OPC_ASD_res_granger_graph %s% lncRTarget_priori_graph
OPC_Normal_res_granger_priori_graph <- OPC_Normal_res_granger_graph %s% lncRTarget_priori_graph

ASTFB_ASD_res_granger_priori_graph <- ASTFB_ASD_res_granger_graph %s% lncRTarget_priori_graph
ASTFB_Normal_res_granger_priori_graph <- ASTFB_Normal_res_granger_graph %s% lncRTarget_priori_graph

Endothelial_ASD_res_granger_priori_graph <- Endothelial_ASD_res_granger_graph %s% lncRTarget_priori_graph
Endothelial_Normal_res_granger_priori_graph <- Endothelial_Normal_res_granger_graph %s% lncRTarget_priori_graph

Microglia_ASD_res_granger_priori_graph <- Microglia_ASD_res_granger_graph %s% lncRTarget_priori_graph
Microglia_Normal_res_granger_priori_graph <- Microglia_Normal_res_granger_graph %s% lncRTarget_priori_graph

NeuNRGNI_ASD_res_granger_priori_graph <- NeuNRGNI_ASD_res_granger_graph %s% lncRTarget_priori_graph
NeuNRGNI_Normal_res_granger_priori_graph <- NeuNRGNI_Normal_res_granger_graph %s% lncRTarget_priori_graph

INVIP_ASD_res_granger_priori_graph <- INVIP_ASD_res_granger_graph %s% lncRTarget_priori_graph
INVIP_Normal_res_granger_priori_graph <- INVIP_Normal_res_granger_graph %s% lncRTarget_priori_graph

L56CC_ASD_res_granger_priori_graph <- L56CC_ASD_res_granger_graph %s% lncRTarget_priori_graph
L56CC_Normal_res_granger_priori_graph <- L56CC_Normal_res_granger_graph %s% lncRTarget_priori_graph

INSV2C_ASD_res_granger_priori_graph <- INSV2C_ASD_res_granger_graph %s% lncRTarget_priori_graph
INSV2C_Normal_res_granger_priori_graph <- INSV2C_Normal_res_granger_graph %s% lncRTarget_priori_graph

L23_ASD_res_granger_priori_graph <- L23_ASD_res_granger_graph %s% lncRTarget_priori_graph
L23_Normal_res_granger_priori_graph <- L23_Normal_res_granger_graph %s% lncRTarget_priori_graph

INPV_ASD_res_granger_priori_graph <- INPV_ASD_res_granger_graph %s% lncRTarget_priori_graph
INPV_Normal_res_granger_priori_graph <- INPV_Normal_res_granger_graph %s% lncRTarget_priori_graph

L4_ASD_res_granger_priori_graph <- L4_ASD_res_granger_graph %s% lncRTarget_priori_graph
L4_Normal_res_granger_priori_graph <- L4_Normal_res_granger_graph %s% lncRTarget_priori_graph

INSST_ASD_res_granger_priori_graph <- INSST_ASD_res_granger_graph %s% lncRTarget_priori_graph
INSST_Normal_res_granger_priori_graph <- INSST_Normal_res_granger_graph %s% lncRTarget_priori_graph

Neumat_ASD_res_granger_priori_graph <- Neumat_ASD_res_granger_graph %s% lncRTarget_priori_graph
Neumat_Normal_res_granger_priori_graph <- Neumat_Normal_res_granger_graph %s% lncRTarget_priori_graph

ASTPP_ASD_res_granger_priori_graph <- ASTPP_ASD_res_granger_graph %s% lncRTarget_priori_graph
ASTPP_Normal_res_granger_priori_graph <- ASTPP_Normal_res_granger_graph %s% lncRTarget_priori_graph

# seqICP
ASD_res_seqICP_priori_graph <- ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Normal_res_seqICP_priori_graph <- Normal_res_seqICP_graph %s% lncRTarget_priori_graph

ACC_ASD_res_seqICP_priori_graph <- ACC_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
ACC_Normal_res_seqICP_priori_graph <- ACC_Normal_res_seqICP_graph %s% lncRTarget_priori_graph
PFC_ASD_res_seqICP_priori_graph <- PFC_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
PFC_Normal_res_seqICP_priori_graph <- PFC_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

Lower18_ASD_res_seqICP_priori_graph <- Lower18_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Lower18_Normal_res_seqICP_priori_graph <- Lower18_Normal_res_seqICP_graph %s% lncRTarget_priori_graph
Larger18_ASD_res_seqICP_priori_graph <- Larger18_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Larger18_Normal_res_seqICP_priori_graph <- Larger18_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

Male_ASD_res_seqICP_priori_graph <- Male_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Male_Normal_res_seqICP_priori_graph <- Male_Normal_res_seqICP_graph %s% lncRTarget_priori_graph
Female_ASD_res_seqICP_priori_graph <- Female_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Female_Normal_res_seqICP_priori_graph <- Female_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

NeuNRGNII_ASD_res_seqICP_priori_graph <- NeuNRGNII_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
NeuNRGNII_Normal_res_seqICP_priori_graph <- NeuNRGNII_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

L56_ASD_res_seqICP_priori_graph <- L56_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
L56_Normal_res_seqICP_priori_graph <- L56_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

Oligodendrocytes_ASD_res_seqICP_priori_graph <- Oligodendrocytes_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Oligodendrocytes_Normal_res_seqICP_priori_graph <- Oligodendrocytes_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

OPC_ASD_res_seqICP_priori_graph <- OPC_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
OPC_Normal_res_seqICP_priori_graph <- OPC_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

ASTFB_ASD_res_seqICP_priori_graph <- ASTFB_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
ASTFB_Normal_res_seqICP_priori_graph <- ASTFB_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

Endothelial_ASD_res_seqICP_priori_graph <- Endothelial_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Endothelial_Normal_res_seqICP_priori_graph <- Endothelial_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

Microglia_ASD_res_seqICP_priori_graph <- Microglia_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Microglia_Normal_res_seqICP_priori_graph <- Microglia_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

NeuNRGNI_ASD_res_seqICP_priori_graph <- NeuNRGNI_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
NeuNRGNI_Normal_res_seqICP_priori_graph <- NeuNRGNI_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

INVIP_ASD_res_seqICP_priori_graph <- INVIP_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
INVIP_Normal_res_seqICP_priori_graph <- INVIP_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

L56CC_ASD_res_seqICP_priori_graph <- L56CC_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
L56CC_Normal_res_seqICP_priori_graph <- L56CC_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

INSV2C_ASD_res_seqICP_priori_graph <- INSV2C_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
INSV2C_Normal_res_seqICP_priori_graph <- INSV2C_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

L23_ASD_res_seqICP_priori_graph <- L23_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
L23_Normal_res_seqICP_priori_graph <- L23_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

INPV_ASD_res_seqICP_priori_graph <- INPV_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
INPV_Normal_res_seqICP_priori_graph <- INPV_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

L4_ASD_res_seqICP_priori_graph <- L4_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
L4_Normal_res_seqICP_priori_graph <- L4_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

INSST_ASD_res_seqICP_priori_graph <- INSST_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
INSST_Normal_res_seqICP_priori_graph <- INSST_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

Neumat_ASD_res_seqICP_priori_graph <- Neumat_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
Neumat_Normal_res_seqICP_priori_graph <- Neumat_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

ASTPP_ASD_res_seqICP_priori_graph <- ASTPP_ASD_res_seqICP_graph %s% lncRTarget_priori_graph
ASTPP_Normal_res_seqICP_priori_graph <- ASTPP_Normal_res_seqICP_graph %s% lncRTarget_priori_graph

########################################################################### Downstream analysis ################################################################################################
############################# 1.1. Validation with priori information #############################
lncRTarget_groundtruth <- as.matrix(read.csv("RegNetwork_high+medium.csv", header = TRUE, sep=","))
lncRTarget_groundtruth_graph <- make_graph(c(t(lncRTarget_groundtruth[, 1:2])), directed = FALSE)

# Clean
ASD_res_darkcausality_priori_validated <- as_data_frame(ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Normal_res_darkcausality_priori_validated <- as_data_frame(Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

ACC_ASD_res_darkcausality_priori_validated <- as_data_frame(ACC_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
ACC_Normal_res_darkcausality_priori_validated <- as_data_frame(ACC_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
PFC_ASD_res_darkcausality_priori_validated <- as_data_frame(PFC_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
PFC_Normal_res_darkcausality_priori_validated <- as_data_frame(PFC_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

Lower18_ASD_res_darkcausality_priori_validated <- as_data_frame(Lower18_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Lower18_Normal_res_darkcausality_priori_validated <- as_data_frame(Lower18_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Larger18_ASD_res_darkcausality_priori_validated <- as_data_frame(Larger18_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Larger18_Normal_res_darkcausality_priori_validated <- as_data_frame(Larger18_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

Male_ASD_res_darkcausality_priori_validated <- as_data_frame(Male_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Male_Normal_res_darkcausality_priori_validated <- as_data_frame(Male_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Female_ASD_res_darkcausality_priori_validated <- as_data_frame(Female_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Female_Normal_res_darkcausality_priori_validated <- as_data_frame(Female_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNII_ASD_res_darkcausality_priori_validated <- as_data_frame(NeuNRGNII_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNII_Normal_res_darkcausality_priori_validated <- as_data_frame(NeuNRGNII_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

L56_ASD_res_darkcausality_priori_validated <- as_data_frame(L56_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
L56_Normal_res_darkcausality_priori_validated <- as_data_frame(L56_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

Oligodendrocytes_ASD_res_darkcausality_priori_validated <- as_data_frame(Oligodendrocytes_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Oligodendrocytes_Normal_res_darkcausality_priori_validated <- as_data_frame(Oligodendrocytes_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

OPC_ASD_res_darkcausality_priori_validated <- as_data_frame(OPC_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
OPC_Normal_res_darkcausality_priori_validated <- as_data_frame(OPC_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

ASTFB_ASD_res_darkcausality_priori_validated <- as_data_frame(ASTFB_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
ASTFB_Normal_res_darkcausality_priori_validated <- as_data_frame(ASTFB_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

Endothelial_ASD_res_darkcausality_priori_validated <- as_data_frame(Endothelial_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Endothelial_Normal_res_darkcausality_priori_validated <- as_data_frame(Endothelial_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

Microglia_ASD_res_darkcausality_priori_validated <- as_data_frame(Microglia_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Microglia_Normal_res_darkcausality_priori_validated <- as_data_frame(Microglia_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNI_ASD_res_darkcausality_priori_validated <- as_data_frame(NeuNRGNI_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNI_Normal_res_darkcausality_priori_validated <- as_data_frame(NeuNRGNI_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

INVIP_ASD_res_darkcausality_priori_validated <- as_data_frame(INVIP_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
INVIP_Normal_res_darkcausality_priori_validated <- as_data_frame(INVIP_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

L56CC_ASD_res_darkcausality_priori_validated <- as_data_frame(L56CC_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
L56CC_Normal_res_darkcausality_priori_validated <- as_data_frame(L56CC_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

INSV2C_ASD_res_darkcausality_priori_validated <- as_data_frame(INSV2C_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
INSV2C_Normal_res_darkcausality_priori_validated <- as_data_frame(INSV2C_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

L23_ASD_res_darkcausality_priori_validated <- as_data_frame(L23_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
L23_Normal_res_darkcausality_priori_validated <- as_data_frame(L23_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

INPV_ASD_res_darkcausality_priori_validated <- as_data_frame(INPV_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
INPV_Normal_res_darkcausality_priori_validated <- as_data_frame(INPV_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

L4_ASD_res_darkcausality_priori_validated <- as_data_frame(L4_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
L4_Normal_res_darkcausality_priori_validated <- as_data_frame(L4_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

INSST_ASD_res_darkcausality_priori_validated <- as_data_frame(INSST_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
INSST_Normal_res_darkcausality_priori_validated <- as_data_frame(INSST_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

Neumat_ASD_res_darkcausality_priori_validated <- as_data_frame(Neumat_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
Neumat_Normal_res_darkcausality_priori_validated <- as_data_frame(Neumat_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

ASTPP_ASD_res_darkcausality_priori_validated <- as_data_frame(ASTPP_ASD_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)
ASTPP_Normal_res_darkcausality_priori_validated <- as_data_frame(ASTPP_Normal_res_darkcausality_priori_graph %s% lncRTarget_groundtruth_graph)

# Granger
ASD_res_granger_priori_validated <- as_data_frame(ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Normal_res_granger_priori_validated <- as_data_frame(Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

ACC_ASD_res_granger_priori_validated <- as_data_frame(ACC_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
ACC_Normal_res_granger_priori_validated <- as_data_frame(ACC_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
PFC_ASD_res_granger_priori_validated <- as_data_frame(PFC_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
PFC_Normal_res_granger_priori_validated <- as_data_frame(PFC_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

Lower18_ASD_res_granger_priori_validated <- as_data_frame(Lower18_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Lower18_Normal_res_granger_priori_validated <- as_data_frame(Lower18_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Larger18_ASD_res_granger_priori_validated <- as_data_frame(Larger18_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Larger18_Normal_res_granger_priori_validated <- as_data_frame(Larger18_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

Male_ASD_res_granger_priori_validated <- as_data_frame(Male_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Male_Normal_res_granger_priori_validated <- as_data_frame(Male_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Female_ASD_res_granger_priori_validated <- as_data_frame(Female_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Female_Normal_res_granger_priori_validated <- as_data_frame(Female_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNII_ASD_res_granger_priori_validated <- as_data_frame(NeuNRGNII_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNII_Normal_res_granger_priori_validated <- as_data_frame(NeuNRGNII_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

L56_ASD_res_granger_priori_validated <- as_data_frame(L56_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
L56_Normal_res_granger_priori_validated <- as_data_frame(L56_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

Oligodendrocytes_ASD_res_granger_priori_validated <- as_data_frame(Oligodendrocytes_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Oligodendrocytes_Normal_res_granger_priori_validated <- as_data_frame(Oligodendrocytes_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

OPC_ASD_res_granger_priori_validated <- as_data_frame(OPC_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
OPC_Normal_res_granger_priori_validated <- as_data_frame(OPC_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

ASTFB_ASD_res_granger_priori_validated <- as_data_frame(ASTFB_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
ASTFB_Normal_res_granger_priori_validated <- as_data_frame(ASTFB_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

Endothelial_ASD_res_granger_priori_validated <- as_data_frame(Endothelial_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Endothelial_Normal_res_granger_priori_validated <- as_data_frame(Endothelial_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

Microglia_ASD_res_granger_priori_validated <- as_data_frame(Microglia_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Microglia_Normal_res_granger_priori_validated <- as_data_frame(Microglia_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNI_ASD_res_granger_priori_validated <- as_data_frame(NeuNRGNI_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNI_Normal_res_granger_priori_validated <- as_data_frame(NeuNRGNI_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

INVIP_ASD_res_granger_priori_validated <- as_data_frame(INVIP_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
INVIP_Normal_res_granger_priori_validated <- as_data_frame(INVIP_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

L56CC_ASD_res_granger_priori_validated <- as_data_frame(L56CC_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
L56CC_Normal_res_granger_priori_validated <- as_data_frame(L56CC_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

INSV2C_ASD_res_granger_priori_validated <- as_data_frame(INSV2C_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
INSV2C_Normal_res_granger_priori_validated <- as_data_frame(INSV2C_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

L23_ASD_res_granger_priori_validated <- as_data_frame(L23_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
L23_Normal_res_granger_priori_validated <- as_data_frame(L23_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

INPV_ASD_res_granger_priori_validated <- as_data_frame(INPV_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
INPV_Normal_res_granger_priori_validated <- as_data_frame(INPV_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

L4_ASD_res_granger_priori_validated <- as_data_frame(L4_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
L4_Normal_res_granger_priori_validated <- as_data_frame(L4_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

INSST_ASD_res_granger_priori_validated <- as_data_frame(INSST_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
INSST_Normal_res_granger_priori_validated <- as_data_frame(INSST_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

Neumat_ASD_res_granger_priori_validated <- as_data_frame(Neumat_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
Neumat_Normal_res_granger_priori_validated <- as_data_frame(Neumat_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

ASTPP_ASD_res_granger_priori_validated <- as_data_frame(ASTPP_ASD_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)
ASTPP_Normal_res_granger_priori_validated <- as_data_frame(ASTPP_Normal_res_granger_priori_graph %s% lncRTarget_groundtruth_graph)

# seqICP
ASD_res_seqICP_priori_validated <- as_data_frame(ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Normal_res_seqICP_priori_validated <- as_data_frame(Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

ACC_ASD_res_seqICP_priori_validated <- as_data_frame(ACC_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
ACC_Normal_res_seqICP_priori_validated <- as_data_frame(ACC_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
PFC_ASD_res_seqICP_priori_validated <- as_data_frame(PFC_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
PFC_Normal_res_seqICP_priori_validated <- as_data_frame(PFC_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

Lower18_ASD_res_seqICP_priori_validated <- as_data_frame(Lower18_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Lower18_Normal_res_seqICP_priori_validated <- as_data_frame(Lower18_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Larger18_ASD_res_seqICP_priori_validated <- as_data_frame(Larger18_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Larger18_Normal_res_seqICP_priori_validated <- as_data_frame(Larger18_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

Male_ASD_res_seqICP_priori_validated <- as_data_frame(Male_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Male_Normal_res_seqICP_priori_validated <- as_data_frame(Male_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Female_ASD_res_seqICP_priori_validated <- as_data_frame(Female_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Female_Normal_res_seqICP_priori_validated <- as_data_frame(Female_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNII_ASD_res_seqICP_priori_validated <- as_data_frame(NeuNRGNII_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNII_Normal_res_seqICP_priori_validated <- as_data_frame(NeuNRGNII_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

L56_ASD_res_seqICP_priori_validated <- as_data_frame(L56_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
L56_Normal_res_seqICP_priori_validated <- as_data_frame(L56_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

Oligodendrocytes_ASD_res_seqICP_priori_validated <- as_data_frame(Oligodendrocytes_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Oligodendrocytes_Normal_res_seqICP_priori_validated <- as_data_frame(Oligodendrocytes_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

OPC_ASD_res_seqICP_priori_validated <- as_data_frame(OPC_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
OPC_Normal_res_seqICP_priori_validated <- as_data_frame(OPC_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

ASTFB_ASD_res_seqICP_priori_validated <- as_data_frame(ASTFB_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
ASTFB_Normal_res_seqICP_priori_validated <- as_data_frame(ASTFB_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

Endothelial_ASD_res_seqICP_priori_validated <- as_data_frame(Endothelial_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Endothelial_Normal_res_seqICP_priori_validated <- as_data_frame(Endothelial_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

Microglia_ASD_res_seqICP_priori_validated <- as_data_frame(Microglia_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Microglia_Normal_res_seqICP_priori_validated <- as_data_frame(Microglia_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNI_ASD_res_seqICP_priori_validated <- as_data_frame(NeuNRGNI_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNI_Normal_res_seqICP_priori_validated <- as_data_frame(NeuNRGNI_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

INVIP_ASD_res_seqICP_priori_validated <- as_data_frame(INVIP_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
INVIP_Normal_res_seqICP_priori_validated <- as_data_frame(INVIP_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

L56CC_ASD_res_seqICP_priori_validated <- as_data_frame(L56CC_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
L56CC_Normal_res_seqICP_priori_validated <- as_data_frame(L56CC_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

INSV2C_ASD_res_seqICP_priori_validated <- as_data_frame(INSV2C_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
INSV2C_Normal_res_seqICP_priori_validated <- as_data_frame(INSV2C_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

L23_ASD_res_seqICP_priori_validated <- as_data_frame(L23_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
L23_Normal_res_seqICP_priori_validated <- as_data_frame(L23_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

INPV_ASD_res_seqICP_priori_validated <- as_data_frame(INPV_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
INPV_Normal_res_seqICP_priori_validated <- as_data_frame(INPV_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

L4_ASD_res_seqICP_priori_validated <- as_data_frame(L4_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
L4_Normal_res_seqICP_priori_validated <- as_data_frame(L4_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

INSST_ASD_res_seqICP_priori_validated <- as_data_frame(INSST_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
INSST_Normal_res_seqICP_priori_validated <- as_data_frame(INSST_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

Neumat_ASD_res_seqICP_priori_validated <- as_data_frame(Neumat_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
Neumat_Normal_res_seqICP_priori_validated <- as_data_frame(Neumat_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

ASTPP_ASD_res_seqICP_priori_validated <- as_data_frame(ASTPP_ASD_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)
ASTPP_Normal_res_seqICP_priori_validated <- as_data_frame(ASTPP_Normal_res_seqICP_priori_graph %s% lncRTarget_groundtruth_graph)

############################# 1.2. Validation without priori information #############################
# Clean
ASD_res_darkcausality_validated <- as_data_frame(ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Normal_res_darkcausality_validated <- as_data_frame(Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

ACC_ASD_res_darkcausality_validated <- as_data_frame(ACC_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
ACC_Normal_res_darkcausality_validated <- as_data_frame(ACC_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
PFC_ASD_res_darkcausality_validated <- as_data_frame(PFC_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
PFC_Normal_res_darkcausality_validated <- as_data_frame(PFC_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

Lower18_ASD_res_darkcausality_validated <- as_data_frame(Lower18_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Lower18_Normal_res_darkcausality_validated <- as_data_frame(Lower18_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Larger18_ASD_res_darkcausality_validated <- as_data_frame(Larger18_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Larger18_Normal_res_darkcausality_validated <- as_data_frame(Larger18_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

Male_ASD_res_darkcausality_validated <- as_data_frame(Male_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Male_Normal_res_darkcausality_validated <- as_data_frame(Male_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Female_ASD_res_darkcausality_validated <- as_data_frame(Female_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Female_Normal_res_darkcausality_validated <- as_data_frame(Female_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNII_ASD_res_darkcausality_validated <- as_data_frame(NeuNRGNII_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNII_Normal_res_darkcausality_validated <- as_data_frame(NeuNRGNII_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

L56_ASD_res_darkcausality_validated <- as_data_frame(L56_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
L56_Normal_res_darkcausality_validated <- as_data_frame(L56_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

Oligodendrocytes_ASD_res_darkcausality_validated <- as_data_frame(Oligodendrocytes_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Oligodendrocytes_Normal_res_darkcausality_validated <- as_data_frame(Oligodendrocytes_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

OPC_ASD_res_darkcausality_validated <- as_data_frame(OPC_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
OPC_Normal_res_darkcausality_validated <- as_data_frame(OPC_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

ASTFB_ASD_res_darkcausality_validated <- as_data_frame(ASTFB_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
ASTFB_Normal_res_darkcausality_validated <- as_data_frame(ASTFB_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

Endothelial_ASD_res_darkcausality_validated <- as_data_frame(Endothelial_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Endothelial_Normal_res_darkcausality_validated <- as_data_frame(Endothelial_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

Microglia_ASD_res_darkcausality_validated <- as_data_frame(Microglia_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Microglia_Normal_res_darkcausality_validated <- as_data_frame(Microglia_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNI_ASD_res_darkcausality_validated <- as_data_frame(NeuNRGNI_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNI_Normal_res_darkcausality_validated <- as_data_frame(NeuNRGNI_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

INVIP_ASD_res_darkcausality_validated <- as_data_frame(INVIP_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
INVIP_Normal_res_darkcausality_validated <- as_data_frame(INVIP_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

L56CC_ASD_res_darkcausality_validated <- as_data_frame(L56CC_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
L56CC_Normal_res_darkcausality_validated <- as_data_frame(L56CC_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

INSV2C_ASD_res_darkcausality_validated <- as_data_frame(INSV2C_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
INSV2C_Normal_res_darkcausality_validated <- as_data_frame(INSV2C_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

L23_ASD_res_darkcausality_validated <- as_data_frame(L23_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
L23_Normal_res_darkcausality_validated <- as_data_frame(L23_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

INPV_ASD_res_darkcausality_validated <- as_data_frame(INPV_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
INPV_Normal_res_darkcausality_validated <- as_data_frame(INPV_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

L4_ASD_res_darkcausality_validated <- as_data_frame(L4_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
L4_Normal_res_darkcausality_validated <- as_data_frame(L4_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

INSST_ASD_res_darkcausality_validated <- as_data_frame(INSST_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
INSST_Normal_res_darkcausality_validated <- as_data_frame(INSST_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

Neumat_ASD_res_darkcausality_validated <- as_data_frame(Neumat_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
Neumat_Normal_res_darkcausality_validated <- as_data_frame(Neumat_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

ASTPP_ASD_res_darkcausality_validated <- as_data_frame(ASTPP_ASD_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)
ASTPP_Normal_res_darkcausality_validated <- as_data_frame(ASTPP_Normal_res_darkcausality_graph %s% lncRTarget_groundtruth_graph)

# Granger
ASD_res_granger_validated <- as_data_frame(ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Normal_res_granger_validated <- as_data_frame(Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

ACC_ASD_res_granger_validated <- as_data_frame(ACC_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
ACC_Normal_res_granger_validated <- as_data_frame(ACC_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)
PFC_ASD_res_granger_validated <- as_data_frame(PFC_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
PFC_Normal_res_granger_validated <- as_data_frame(PFC_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

Lower18_ASD_res_granger_validated <- as_data_frame(Lower18_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Lower18_Normal_res_granger_validated <- as_data_frame(Lower18_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)
Larger18_ASD_res_granger_validated <- as_data_frame(Larger18_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Larger18_Normal_res_granger_validated <- as_data_frame(Larger18_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

Male_ASD_res_granger_validated <- as_data_frame(Male_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Male_Normal_res_granger_validated <- as_data_frame(Male_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)
Female_ASD_res_granger_validated <- as_data_frame(Female_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Female_Normal_res_granger_validated <- as_data_frame(Female_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNII_ASD_res_granger_validated <- as_data_frame(NeuNRGNII_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNII_Normal_res_granger_validated <- as_data_frame(NeuNRGNII_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

L56_ASD_res_granger_validated <- as_data_frame(L56_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
L56_Normal_res_granger_validated <- as_data_frame(L56_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

Oligodendrocytes_ASD_res_granger_validated <- as_data_frame(Oligodendrocytes_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Oligodendrocytes_Normal_res_granger_validated <- as_data_frame(Oligodendrocytes_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

OPC_ASD_res_granger_validated <- as_data_frame(OPC_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
OPC_Normal_res_granger_validated <- as_data_frame(OPC_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

ASTFB_ASD_res_granger_validated <- as_data_frame(ASTFB_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
ASTFB_Normal_res_granger_validated <- as_data_frame(ASTFB_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

Endothelial_ASD_res_granger_validated <- as_data_frame(Endothelial_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Endothelial_Normal_res_granger_validated <- as_data_frame(Endothelial_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

Microglia_ASD_res_granger_validated <- as_data_frame(Microglia_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Microglia_Normal_res_granger_validated <- as_data_frame(Microglia_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNI_ASD_res_granger_validated <- as_data_frame(NeuNRGNI_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNI_Normal_res_granger_validated <- as_data_frame(NeuNRGNI_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

INVIP_ASD_res_granger_validated <- as_data_frame(INVIP_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
INVIP_Normal_res_granger_validated <- as_data_frame(INVIP_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

L56CC_ASD_res_granger_validated <- as_data_frame(L56CC_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
L56CC_Normal_res_granger_validated <- as_data_frame(L56CC_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

INSV2C_ASD_res_granger_validated <- as_data_frame(INSV2C_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
INSV2C_Normal_res_granger_validated <- as_data_frame(INSV2C_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

L23_ASD_res_granger_validated <- as_data_frame(L23_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
L23_Normal_res_granger_validated <- as_data_frame(L23_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

INPV_ASD_res_granger_validated <- as_data_frame(INPV_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
INPV_Normal_res_granger_validated <- as_data_frame(INPV_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

L4_ASD_res_granger_validated <- as_data_frame(L4_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
L4_Normal_res_granger_validated <- as_data_frame(L4_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

INSST_ASD_res_granger_validated <- as_data_frame(INSST_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
INSST_Normal_res_granger_validated <- as_data_frame(INSST_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

Neumat_ASD_res_granger_validated <- as_data_frame(Neumat_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
Neumat_Normal_res_granger_validated <- as_data_frame(Neumat_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

ASTPP_ASD_res_granger_validated <- as_data_frame(ASTPP_ASD_res_granger_graph %s% lncRTarget_groundtruth_graph)
ASTPP_Normal_res_granger_validated <- as_data_frame(ASTPP_Normal_res_granger_graph %s% lncRTarget_groundtruth_graph)

# seqICP
ASD_res_seqICP_validated <- as_data_frame(ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Normal_res_seqICP_validated <- as_data_frame(Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

ACC_ASD_res_seqICP_validated <- as_data_frame(ACC_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
ACC_Normal_res_seqICP_validated <- as_data_frame(ACC_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
PFC_ASD_res_seqICP_validated <- as_data_frame(PFC_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
PFC_Normal_res_seqICP_validated <- as_data_frame(PFC_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

Lower18_ASD_res_seqICP_validated <- as_data_frame(Lower18_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Lower18_Normal_res_seqICP_validated <- as_data_frame(Lower18_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Larger18_ASD_res_seqICP_validated <- as_data_frame(Larger18_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Larger18_Normal_res_seqICP_validated <- as_data_frame(Larger18_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

Male_ASD_res_seqICP_validated <- as_data_frame(Male_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Male_Normal_res_seqICP_validated <- as_data_frame(Male_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Female_ASD_res_seqICP_validated <- as_data_frame(Female_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Female_Normal_res_seqICP_validated <- as_data_frame(Female_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNII_ASD_res_seqICP_validated <- as_data_frame(NeuNRGNII_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNII_Normal_res_seqICP_validated <- as_data_frame(NeuNRGNII_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

L56_ASD_res_seqICP_validated <- as_data_frame(L56_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
L56_Normal_res_seqICP_validated <- as_data_frame(L56_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

Oligodendrocytes_ASD_res_seqICP_validated <- as_data_frame(Oligodendrocytes_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Oligodendrocytes_Normal_res_seqICP_validated <- as_data_frame(Oligodendrocytes_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

OPC_ASD_res_seqICP_validated <- as_data_frame(OPC_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
OPC_Normal_res_seqICP_validated <- as_data_frame(OPC_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

ASTFB_ASD_res_seqICP_validated <- as_data_frame(ASTFB_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
ASTFB_Normal_res_seqICP_validated <- as_data_frame(ASTFB_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

Endothelial_ASD_res_seqICP_validated <- as_data_frame(Endothelial_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Endothelial_Normal_res_seqICP_validated <- as_data_frame(Endothelial_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

Microglia_ASD_res_seqICP_validated <- as_data_frame(Microglia_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Microglia_Normal_res_seqICP_validated <- as_data_frame(Microglia_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

NeuNRGNI_ASD_res_seqICP_validated <- as_data_frame(NeuNRGNI_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
NeuNRGNI_Normal_res_seqICP_validated <- as_data_frame(NeuNRGNI_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

INVIP_ASD_res_seqICP_validated <- as_data_frame(INVIP_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
INVIP_Normal_res_seqICP_validated <- as_data_frame(INVIP_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

L56CC_ASD_res_seqICP_validated <- as_data_frame(L56CC_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
L56CC_Normal_res_seqICP_validated <- as_data_frame(L56CC_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

INSV2C_ASD_res_seqICP_validated <- as_data_frame(INSV2C_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
INSV2C_Normal_res_seqICP_validated <- as_data_frame(INSV2C_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

L23_ASD_res_seqICP_validated <- as_data_frame(L23_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
L23_Normal_res_seqICP_validated <- as_data_frame(L23_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

INPV_ASD_res_seqICP_validated <- as_data_frame(INPV_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
INPV_Normal_res_seqICP_validated <- as_data_frame(INPV_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

L4_ASD_res_seqICP_validated <- as_data_frame(L4_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
L4_Normal_res_seqICP_validated <- as_data_frame(L4_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

INSST_ASD_res_seqICP_validated <- as_data_frame(INSST_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
INSST_Normal_res_seqICP_validated <- as_data_frame(INSST_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

Neumat_ASD_res_seqICP_validated <- as_data_frame(Neumat_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
Neumat_Normal_res_seqICP_validated <- as_data_frame(Neumat_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

ASTPP_ASD_res_seqICP_validated <- as_data_frame(ASTPP_ASD_res_seqICP_graph %s% lncRTarget_groundtruth_graph)
ASTPP_Normal_res_seqICP_validated <- as_data_frame(ASTPP_Normal_res_seqICP_graph %s% lncRTarget_groundtruth_graph)

############################# 2. Heterogeneity analysis #########################
lncR_diagnostic_net <- list(ASD_res_darkcausality_priori_graph, Normal_res_darkcausality_priori_graph)
lncR_region_net <- list(ACC_ASD_res_darkcausality_priori_graph,
                        PFC_ASD_res_darkcausality_priori_graph,
                        ACC_Normal_res_darkcausality_priori_graph,                        
                        PFC_Normal_res_darkcausality_priori_graph)
lncR_age_net <- list(Lower18_ASD_res_darkcausality_priori_graph,
                     Larger18_ASD_res_darkcausality_priori_graph,
                     Lower18_Normal_res_darkcausality_priori_graph,
		     Larger18_Normal_res_darkcausality_priori_graph)
lncR_sex_net <- list(Male_ASD_res_darkcausality_priori_graph,
                     Female_ASD_res_darkcausality_priori_graph,
		     Male_Normal_res_darkcausality_priori_graph,
		     Female_Normal_res_darkcausality_priori_graph)
lncR_celltype_net <- list(NeuNRGNII_ASD_res_darkcausality_priori_graph,
                     L56_ASD_res_darkcausality_priori_graph,
                     Oligodendrocytes_ASD_res_darkcausality_priori_graph,
		     OPC_ASD_res_darkcausality_priori_graph,
                     ASTFB_ASD_res_darkcausality_priori_graph,
                     Endothelial_ASD_res_darkcausality_priori_graph,
                     Microglia_ASD_res_darkcausality_priori_graph,
		     NeuNRGNI_ASD_res_darkcausality_priori_graph,
		     INVIP_ASD_res_darkcausality_priori_graph,
		     L56CC_ASD_res_darkcausality_priori_graph,
		     INSV2C_ASD_res_darkcausality_priori_graph,
		     L23_ASD_res_darkcausality_priori_graph,
		     INPV_ASD_res_darkcausality_priori_graph,
		     L4_ASD_res_darkcausality_priori_graph,
		     INSST_ASD_res_darkcausality_priori_graph,
		     Neumat_ASD_res_darkcausality_priori_graph,
		     ASTPP_ASD_res_darkcausality_priori_graph,
                     NeuNRGNII_Normal_res_darkcausality_priori_graph,
                     L56_Normal_res_darkcausality_priori_graph,
                     Oligodendrocytes_Normal_res_darkcausality_priori_graph,
		     OPC_Normal_res_darkcausality_priori_graph,
                     ASTFB_Normal_res_darkcausality_priori_graph,
                     Endothelial_Normal_res_darkcausality_priori_graph,
                     Microglia_Normal_res_darkcausality_priori_graph,
		     NeuNRGNI_Normal_res_darkcausality_priori_graph,
		     INVIP_Normal_res_darkcausality_priori_graph,
		     L56CC_Normal_res_darkcausality_priori_graph,
		     INSV2C_Normal_res_darkcausality_priori_graph,
		     L23_Normal_res_darkcausality_priori_graph,
		     INPV_Normal_res_darkcausality_priori_graph,
		     L4_Normal_res_darkcausality_priori_graph,
		     INSST_Normal_res_darkcausality_priori_graph,
		     Neumat_Normal_res_darkcausality_priori_graph,
		     ASTPP_Normal_res_darkcausality_priori_graph
                     )
diagnostic_Het <- Het.network(lncR_diagnostic_net, lncR_diagnostic_net)
region_Het <- Het.network(lncR_region_net, lncR_region_net)
age_Het <- Het.network(lncR_age_net, lncR_age_net)
sex_Het <- Het.network(lncR_sex_net, lncR_sex_net)
celltype_Het <- Het.network(lncR_celltype_net, lncR_celltype_net)

############################# 3. ASD lncRNAs identification and enrichemnt analysis #############################
Dif_res_darkcausality_priori_graph <- (ASD_res_darkcausality_priori_graph %m% Normal_res_darkcausality_priori_graph) %u% (Normal_res_darkcausality_priori_graph %m% ASD_res_darkcausality_priori_graph)
Dif_lncRNA_list <- unique(tail_of(Dif_res_darkcausality_priori_graph, E(Dif_res_darkcausality_priori_graph))$name)
Dif_mRNA_list <- lapply(seq(Dif_lncRNA_list), function(i) neighbors(Dif_res_darkcausality_priori_graph, Dif_lncRNA_list[i])$name)

g_asd <- as_directed(ASD_res_darkcausality_priori_graph, mode = "arbitrary")
g_normal <- as_directed(Normal_res_darkcausality_priori_graph, mode = "arbitrary")

# Calculate edge rewiring
rewiring_results <- calculate_lncrna_rewiring(g_asd, g_normal)

# Calculate risk scores
risk_scores <- calculate_risk_scores_zscore(rewiring_results$rewiring_df)

# Identify ASD risk lncRNAs and their targets
significant_lncRNAs <- identify_significant_risk_lncrnas(risk_scores)
ASD_risk_lncRNAs <- significant_lncRNAs[significant_lncRNAs$significant, 1]
ASD_risk_lncRNAs_targets <- lapply(seq(ASD_risk_lncRNAs), function(i) Dif_mRNA_list[[which(Dif_lncRNA_list %in% ASD_risk_lncRNAs[i])]])

# Enrichment analysis of ASD risk lncRNA targets
enrichment_FEA_BP <- module_FA(ASD_risk_lncRNAs_targets, GOont = "BP", Analysis.type = "FEA")
enrichment_FEA_MF <- module_FA(ASD_risk_lncRNAs_targets, GOont = "MF", Analysis.type = "FEA")
enrichment_FEA_CC <- module_FA(ASD_risk_lncRNAs_targets, GOont = "CC", Analysis.type = "FEA")
enrichment_DEA <- module_FA(ASD_risk_lncRNAs_targets, Analysis.type = "DEA")

enrichment_FEA_BP_dataframe <- lapply(seq(enrichment_FEA_BP[[1]]), function(i) as.data.frame(enrichment_FEA_BP[[1]][[i]]))
enrichment_FEA_BP_num <- lapply(seq(enrichment_FEA_BP_dataframe), function(i) nrow(enrichment_FEA_BP_dataframe[[i]]))
enrichment_FEA_BP_res <- do.call(rbind, enrichment_FEA_BP_dataframe)

enrichment_FEA_MF_dataframe <- lapply(seq(enrichment_FEA_MF[[1]]), function(i) as.data.frame(enrichment_FEA_MF[[1]][[i]]))
enrichment_FEA_MF_num <- lapply(seq(enrichment_FEA_MF_dataframe), function(i) nrow(enrichment_FEA_MF_dataframe[[i]]))
enrichment_FEA_MF_res <- do.call(rbind, enrichment_FEA_MF_dataframe)

enrichment_FEA_CC_dataframe <- lapply(seq(enrichment_FEA_CC[[1]]), function(i) as.data.frame(enrichment_FEA_CC[[1]][[i]]))
enrichment_FEA_CC_num <- lapply(seq(enrichment_FEA_CC_dataframe), function(i) nrow(enrichment_FEA_CC_dataframe[[i]]))
enrichment_FEA_CC_res <- do.call(rbind, enrichment_FEA_CC_dataframe)

enrichment_FEA_KEGG_dataframe <- lapply(seq(enrichment_FEA_BP[[2]]), function(i) as.data.frame(enrichment_FEA_BP[[2]][[i]]))
enrichment_FEA_KEGG_num <- lapply(seq(enrichment_FEA_KEGG_dataframe), function(i) nrow(enrichment_FEA_KEGG_dataframe[[i]]))
enrichment_FEA_KEGG_res <- do.call(rbind, enrichment_FEA_KEGG_dataframe)

enrichment_FEA_Reactome_dataframe <- lapply(seq(enrichment_FEA_BP[[3]]), function(i) as.data.frame(enrichment_FEA_BP[[3]][[i]]))
enrichment_FEA_Reactome_num <- lapply(seq(enrichment_FEA_Reactome_dataframe), function(i) nrow(enrichment_FEA_Reactome_dataframe[[i]]))
enrichment_FEA_Reactome_res <- do.call(rbind, enrichment_FEA_Reactome_dataframe)

############################# 4. lncRNA biomarkers prediction (context-specific and context-enhanced biomarkers) and classification performance #############################
# Diagnostic: ASD and normal
diagnostic_graph <- list(ASD_res_darkcausality_priori_graph, Normal_res_darkcausality_priori_graph)
diagnostic_degree <- lapply(seq(diagnostic_graph), function(i) degree(diagnostic_graph[[i]]))
n_diagnostic <- lapply(seq(diagnostic_graph), function(i) length(diagnostic_degree[[i]]))
m_diagnostic <- lapply(seq(diagnostic_graph), function(i) ecount(diagnostic_graph[[i]]))
p1_diagnostic <- lapply(seq(diagnostic_graph), function(i) m_diagnostic[[i]]/(n_diagnostic[[i]]^2))
lamda_diagnostic <- lapply(seq(diagnostic_graph), function(i) n_diagnostic[[i]]*p1_diagnostic[[i]])
p_diagnostic <- lapply(seq(diagnostic_graph), function(i) ppois(diagnostic_degree[[i]] - 1, lambda = lamda_diagnostic[[i]], lower=FALSE))
hub_diagnostic <- lapply(seq(diagnostic_graph), function(i) p_diagnostic[[i]][which(p_diagnostic[[i]] < 0.05)])
hub_lncRNA_diagnostic <- lapply(seq(diagnostic_graph), function(i) hub_diagnostic[[i]][which(names(hub_diagnostic[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_diagnostic_specific <- Overlap.hub(list(names(hub_lncRNA_diagnostic[[1]]), names(hub_lncRNA_diagnostic[[2]])), overlap.num = 1, type = "equal")
biomarker_diagnostic_2 <- Overlap.hub(list(names(hub_lncRNA_diagnostic[[1]]), names(hub_lncRNA_diagnostic[[2]])), overlap.num = 2, type = "equal")
biomarker_diagnostic_2_fold <- diagnostic_degree[[1]][which(names(diagnostic_degree[[1]]) %in% biomarker_diagnostic_2)]/diagnostic_degree[[2]][which(names(diagnostic_degree[[2]]) %in% biomarker_diagnostic_2)]
biomarker_diagnostic_enhanced <- names(which(biomarker_diagnostic_2_fold > 2 | biomarker_diagnostic_2_fold < 0.5))

# Region: ACC_ASD, ACC_Normal, PFC_ASD, PFC_Normal
region_ASD_graph <- list(ACC_ASD_res_darkcausality_priori_graph, PFC_ASD_res_darkcausality_priori_graph)
region_ASD_degree <- lapply(seq(region_ASD_graph), function(i) degree(region_ASD_graph[[i]]))
n_region_ASD <- lapply(seq(region_ASD_graph), function(i) length(region_ASD_degree[[i]]))
m_region_ASD <- lapply(seq(region_ASD_graph), function(i) ecount(region_ASD_graph[[i]]))
p1_region_ASD <- lapply(seq(region_ASD_graph), function(i) m_region_ASD[[i]]/(n_region_ASD[[i]]^2))
lamda_region_ASD <- lapply(seq(region_ASD_graph), function(i) n_region_ASD[[i]]*p1_region_ASD[[i]])
p_region_ASD <- lapply(seq(region_ASD_graph), function(i) ppois(region_ASD_degree[[i]] - 1, lambda = lamda_region_ASD[[i]], lower=FALSE))
hub_region_ASD <- lapply(seq(region_ASD_graph), function(i) p_region_ASD[[i]][which(p_region_ASD[[i]] < 0.05)])
hub_lncRNA_region_ASD <- lapply(seq(region_ASD_graph), function(i) hub_region_ASD[[i]][which(names(hub_region_ASD[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_region_ASD_specific <- Overlap.hub(list(names(hub_lncRNA_region_ASD[[1]]), names(hub_lncRNA_region_ASD[[2]])), overlap.num = 1, type = "equal")
biomarker_region_ASD_2 <- Overlap.hub(list(names(hub_lncRNA_region_ASD[[1]]), names(hub_lncRNA_region_ASD[[2]])), overlap.num = 2, type = "equal")
biomarker_region_ASD_2_fold <- region_ASD_degree[[1]][which(names(region_ASD_degree[[1]]) %in% biomarker_region_ASD_2)]/region_ASD_degree[[2]][which(names(region_ASD_degree[[2]]) %in% biomarker_region_ASD_2)]
biomarker_region_ASD_enhanced <- names(which(biomarker_region_ASD_2_fold > 2 | biomarker_region_ASD_2_fold < 0.5))

region_Normal_graph <- list(ACC_Normal_res_darkcausality_priori_graph, PFC_Normal_res_darkcausality_priori_graph)
region_Normal_degree <- lapply(seq(region_Normal_graph), function(i) degree(region_Normal_graph[[i]]))
n_region_Normal <- lapply(seq(region_Normal_graph), function(i) length(region_Normal_degree[[i]]))
m_region_Normal <- lapply(seq(region_Normal_graph), function(i) ecount(region_Normal_graph[[i]]))
p1_region_Normal <- lapply(seq(region_Normal_graph), function(i) m_region_Normal[[i]]/(n_region_Normal[[i]]^2))
lamda_region_Normal <- lapply(seq(region_Normal_graph), function(i) n_region_Normal[[i]]*p1_region_Normal[[i]])
p_region_Normal <- lapply(seq(region_Normal_graph), function(i) ppois(region_Normal_degree[[i]] - 1, lambda = lamda_region_Normal[[i]], lower=FALSE))
hub_region_Normal <- lapply(seq(region_Normal_graph), function(i) p_region_Normal[[i]][which(p_region_Normal[[i]] < 0.05)])
hub_lncRNA_region_Normal <- lapply(seq(region_Normal_graph), function(i) hub_region_Normal[[i]][which(names(hub_region_Normal[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_region_Normal_specific <- Overlap.hub(list(names(hub_lncRNA_region_Normal[[1]]), names(hub_lncRNA_region_Normal[[2]])), overlap.num = 1, type = "equal")
biomarker_region_Normal_2 <- Overlap.hub(list(names(hub_lncRNA_region_Normal[[1]]), names(hub_lncRNA_region_Normal[[2]])), overlap.num = 2, type = "equal")
biomarker_region_Normal_2_fold <- region_Normal_degree[[1]][which(names(region_Normal_degree[[1]]) %in% biomarker_region_Normal_2)]/region_Normal_degree[[2]][which(names(region_Normal_degree[[2]]) %in% biomarker_region_Normal_2)]
biomarker_region_Normal_enhanced <- names(which(biomarker_region_Normal_2_fold > 2 | biomarker_region_Normal_2_fold < 0.5))

# Age: Lower18_ASD, Lower18_Normal, Larger18_ASD, Larger18_Normal
age_ASD_graph <- list(Lower18_ASD_res_darkcausality_priori_graph, Larger18_ASD_res_darkcausality_priori_graph)
age_ASD_degree <- lapply(seq(age_ASD_graph), function(i) degree(age_ASD_graph[[i]]))
n_age_ASD <- lapply(seq(age_ASD_graph), function(i) length(age_ASD_degree[[i]]))
m_age_ASD <- lapply(seq(age_ASD_graph), function(i) ecount(age_ASD_graph[[i]]))
p1_age_ASD <- lapply(seq(age_ASD_graph), function(i) m_age_ASD[[i]]/(n_age_ASD[[i]]^2))
lamda_age_ASD <- lapply(seq(age_ASD_graph), function(i) n_age_ASD[[i]]*p1_age_ASD[[i]])
p_age_ASD <- lapply(seq(age_ASD_graph), function(i) ppois(age_ASD_degree[[i]] - 1, lambda = lamda_age_ASD[[i]], lower=FALSE))
hub_age_ASD <- lapply(seq(age_ASD_graph), function(i) p_age_ASD[[i]][which(p_age_ASD[[i]] < 0.05)])
hub_lncRNA_age_ASD <- lapply(seq(age_ASD_graph), function(i) hub_age_ASD[[i]][which(names(hub_age_ASD[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_age_ASD_specific <- Overlap.hub(list(names(hub_lncRNA_age_ASD[[1]]), names(hub_lncRNA_age_ASD[[2]])), overlap.num = 1, type = "equal")
biomarker_age_ASD_2 <- Overlap.hub(list(names(hub_lncRNA_age_ASD[[1]]), names(hub_lncRNA_age_ASD[[2]])), overlap.num = 2, type = "equal")
biomarker_age_ASD_2_fold <- age_ASD_degree[[1]][which(names(age_ASD_degree[[1]]) %in% biomarker_age_ASD_2)]/age_ASD_degree[[2]][which(names(age_ASD_degree[[2]]) %in% biomarker_age_ASD_2)]
biomarker_age_ASD_enhanced <- names(which(biomarker_age_ASD_2_fold > 2 | biomarker_age_ASD_2_fold < 0.5))

age_Normal_graph <- list(Lower18_Normal_res_darkcausality_priori_graph, Larger18_Normal_res_darkcausality_priori_graph)
age_Normal_degree <- lapply(seq(age_Normal_graph), function(i) degree(age_Normal_graph[[i]]))
n_age_Normal <- lapply(seq(age_Normal_graph), function(i) length(age_Normal_degree[[i]]))
m_age_Normal <- lapply(seq(age_Normal_graph), function(i) ecount(age_Normal_graph[[i]]))
p1_age_Normal <- lapply(seq(age_Normal_graph), function(i) m_age_Normal[[i]]/(n_age_Normal[[i]]^2))
lamda_age_Normal <- lapply(seq(age_Normal_graph), function(i) n_age_Normal[[i]]*p1_age_Normal[[i]])
p_age_Normal <- lapply(seq(age_Normal_graph), function(i) ppois(age_Normal_degree[[i]] - 1, lambda = lamda_age_Normal[[i]], lower=FALSE))
hub_age_Normal <- lapply(seq(age_Normal_graph), function(i) p_age_Normal[[i]][which(p_age_Normal[[i]] < 0.05)])
hub_lncRNA_age_Normal <- lapply(seq(age_Normal_graph), function(i) hub_age_Normal[[i]][which(names(hub_age_Normal[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_age_Normal_specific <- Overlap.hub(list(names(hub_lncRNA_age_Normal[[1]]), names(hub_lncRNA_age_Normal[[2]])), overlap.num = 1, type = "equal")
biomarker_age_Normal_2 <- Overlap.hub(list(names(hub_lncRNA_age_Normal[[1]]), names(hub_lncRNA_age_Normal[[2]])), overlap.num = 2, type = "equal")
biomarker_age_Normal_2_fold <- age_Normal_degree[[1]][which(names(age_Normal_degree[[1]]) %in% biomarker_age_Normal_2)]/age_Normal_degree[[2]][which(names(age_Normal_degree[[2]]) %in% biomarker_age_Normal_2)]
biomarker_age_Normal_enhanced <- names(which(biomarker_age_Normal_2_fold > 2 | biomarker_age_Normal_2_fold < 0.5))

# Sex: Male_ASD, Male_Normal, Female_ASD, Female_Normal
sex_ASD_graph <- list(Male_ASD_res_darkcausality_priori_graph, Female_ASD_res_darkcausality_priori_graph)
sex_ASD_degree <- lapply(seq(sex_ASD_graph), function(i) degree(sex_ASD_graph[[i]]))
n_sex_ASD <- lapply(seq(sex_ASD_graph), function(i) length(sex_ASD_degree[[i]]))
m_sex_ASD <- lapply(seq(sex_ASD_graph), function(i) ecount(sex_ASD_graph[[i]]))
p1_sex_ASD <- lapply(seq(sex_ASD_graph), function(i) m_sex_ASD[[i]]/(n_sex_ASD[[i]]^2))
lamda_sex_ASD <- lapply(seq(sex_ASD_graph), function(i) n_sex_ASD[[i]]*p1_sex_ASD[[i]])
p_sex_ASD <- lapply(seq(sex_ASD_graph), function(i) ppois(sex_ASD_degree[[i]] - 1, lambda = lamda_sex_ASD[[i]], lower=FALSE))
hub_sex_ASD <- lapply(seq(sex_ASD_graph), function(i) p_sex_ASD[[i]][which(p_sex_ASD[[i]] < 0.05)])
hub_lncRNA_sex_ASD <- lapply(seq(sex_ASD_graph), function(i) hub_sex_ASD[[i]][which(names(hub_sex_ASD[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_sex_ASD_specific <- Overlap.hub(list(names(hub_lncRNA_sex_ASD[[1]]), names(hub_lncRNA_sex_ASD[[2]])), overlap.num = 1, type = "equal")
biomarker_sex_ASD_2 <- Overlap.hub(list(names(hub_lncRNA_sex_ASD[[1]]), names(hub_lncRNA_sex_ASD[[2]])), overlap.num = 2, type = "equal")
biomarker_sex_ASD_2_fold <- sex_ASD_degree[[1]][which(names(sex_ASD_degree[[1]]) %in% biomarker_sex_ASD_2)]/sex_ASD_degree[[2]][which(names(sex_ASD_degree[[2]]) %in% biomarker_sex_ASD_2)]
biomarker_sex_ASD_enhanced <- names(which(biomarker_sex_ASD_2_fold > 2 | biomarker_sex_ASD_2_fold < 0.5))

sex_Normal_graph <- list(Male_Normal_res_darkcausality_priori_graph, Female_Normal_res_darkcausality_priori_graph)
sex_Normal_degree <- lapply(seq(sex_Normal_graph), function(i) degree(sex_Normal_graph[[i]]))
n_sex_Normal <- lapply(seq(sex_Normal_graph), function(i) length(sex_Normal_degree[[i]]))
m_sex_Normal <- lapply(seq(sex_Normal_graph), function(i) ecount(sex_Normal_graph[[i]]))
p1_sex_Normal <- lapply(seq(sex_Normal_graph), function(i) m_sex_Normal[[i]]/(n_sex_Normal[[i]]^2))
lamda_sex_Normal <- lapply(seq(sex_Normal_graph), function(i) n_sex_Normal[[i]]*p1_sex_Normal[[i]])
p_sex_Normal <- lapply(seq(sex_Normal_graph), function(i) ppois(sex_Normal_degree[[i]] - 1, lambda = lamda_sex_Normal[[i]], lower=FALSE))
hub_sex_Normal <- lapply(seq(sex_Normal_graph), function(i) p_sex_Normal[[i]][which(p_sex_Normal[[i]] < 0.05)])
hub_lncRNA_sex_Normal <- lapply(seq(sex_Normal_graph), function(i) hub_sex_Normal[[i]][which(names(hub_sex_Normal[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_sex_Normal_specific <- Overlap.hub(list(names(hub_lncRNA_sex_Normal[[1]]), names(hub_lncRNA_sex_Normal[[2]])), overlap.num = 1, type = "equal")
biomarker_sex_Normal_2 <- Overlap.hub(list(names(hub_lncRNA_sex_Normal[[1]]), names(hub_lncRNA_sex_Normal[[2]])), overlap.num = 2, type = "equal")
biomarker_sex_Normal_2_fold <- sex_Normal_degree[[1]][which(names(sex_Normal_degree[[1]]) %in% biomarker_sex_Normal_2)]/sex_Normal_degree[[2]][which(names(sex_Normal_degree[[2]]) %in% biomarker_sex_Normal_2)]
biomarker_sex_Normal_enhanced <- names(which(biomarker_sex_Normal_2_fold > 2 | biomarker_sex_Normal_2_fold < 0.5))

# Cell type: NeuNRGNII_ASD, NeuNRGNII_Normal, L56_ASD, L56_Normal, Oligodendrocytes_ASD, Oligodendrocytes_Normal, OPC_ASD, OPC_Normal,
# ASTFB_ASD, ASTFB_Normal, Endothelial_ASD, Endothelial_Normal, Microglia_ASD, Microglia_Normal, NeuNRGNI_ASD, NeuNRGNI_Normal,
# INVIP_ASD, INVIP_Normal, L56CC_ASD, L56CC_Normal, INSV2C_ASD, INSV2C_Normal, L23_ASD, L23_Normal, INPV_ASD, INPV_Normal, L4_ASD, L4_Normal,
# INSST_ASD, INSST_Normal, Neumat_ASD, Neumat_Normal, ASTPP_ASD, ASTPP_Normal
celltype_ASD_graph <- list(NeuNRGNII_ASD_res_darkcausality_priori_graph, 
                       L56_ASD_res_darkcausality_priori_graph,
                       Oligodendrocytes_ASD_res_darkcausality_priori_graph,
                       OPC_ASD_res_darkcausality_priori_graph,
                       ASTFB_ASD_res_darkcausality_priori_graph,
                       Endothelial_ASD_res_darkcausality_priori_graph,
                       Microglia_ASD_res_darkcausality_priori_graph,
                       NeuNRGNI_ASD_res_darkcausality_priori_graph,
                       INVIP_ASD_res_darkcausality_priori_graph,
                       L56CC_ASD_res_darkcausality_priori_graph,
                       INSV2C_ASD_res_darkcausality_priori_graph,
                       L23_ASD_res_darkcausality_priori_graph,
                       INPV_ASD_res_darkcausality_priori_graph,
                       L4_ASD_res_darkcausality_priori_graph,
                       INSST_ASD_res_darkcausality_priori_graph,
                       Neumat_ASD_res_darkcausality_priori_graph,
                       ASTPP_ASD_res_darkcausality_priori_graph)
celltype_ASD_degree <- lapply(seq(celltype_ASD_graph), function(i) degree(celltype_ASD_graph[[i]]))
n_celltype_ASD <- lapply(seq(celltype_ASD_graph), function(i) length(celltype_ASD_degree[[i]]))
m_celltype_ASD <- lapply(seq(celltype_ASD_graph), function(i) ecount(celltype_ASD_graph[[i]]))
p1_celltype_ASD <- lapply(seq(celltype_ASD_graph), function(i) m_celltype_ASD[[i]]/(n_celltype_ASD[[i]]^2))
lamda_celltype_ASD <- lapply(seq(celltype_ASD_graph), function(i) n_celltype_ASD[[i]]*p1_celltype_ASD[[i]])
p_celltype_ASD <- lapply(seq(celltype_ASD_graph), function(i) ppois(celltype_ASD_degree[[i]] - 1, lambda = lamda_celltype_ASD[[i]], lower=FALSE))
hub_celltype_ASD <- lapply(seq(celltype_ASD_graph), function(i) p_celltype_ASD[[i]][which(p_celltype_ASD[[i]] < 0.05)])
hub_lncRNA_celltype_ASD <- lapply(seq(celltype_ASD_graph), function(i) hub_celltype_ASD[[i]][which(names(hub_celltype_ASD[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_celltype_ASD_specific <- Overlap.hub(list(names(hub_lncRNA_celltype_ASD[[1]]), 
                                                    names(hub_lncRNA_celltype_ASD[[2]]),
						    names(hub_lncRNA_celltype_ASD[[3]]),
						    names(hub_lncRNA_celltype_ASD[[4]]),
						    names(hub_lncRNA_celltype_ASD[[5]]),
						    names(hub_lncRNA_celltype_ASD[[6]]),
						    names(hub_lncRNA_celltype_ASD[[7]]),
						    names(hub_lncRNA_celltype_ASD[[8]]),
						    names(hub_lncRNA_celltype_ASD[[9]]),
						    names(hub_lncRNA_celltype_ASD[[10]]),
						    names(hub_lncRNA_celltype_ASD[[11]]),
						    names(hub_lncRNA_celltype_ASD[[12]]),
						    names(hub_lncRNA_celltype_ASD[[13]]),
						    names(hub_lncRNA_celltype_ASD[[14]]),
						    names(hub_lncRNA_celltype_ASD[[15]]),
						    names(hub_lncRNA_celltype_ASD[[16]]),
						    names(hub_lncRNA_celltype_ASD[[17]])), 
						    overlap.num = 1, type = "equal")
biomarker_celltype_ASD_2 <- Overlap.hub(list(names(hub_lncRNA_celltype_ASD[[1]]), 
                                                    names(hub_lncRNA_celltype_ASD[[2]]),
						    names(hub_lncRNA_celltype_ASD[[3]]),
						    names(hub_lncRNA_celltype_ASD[[4]]),
						    names(hub_lncRNA_celltype_ASD[[5]]),
						    names(hub_lncRNA_celltype_ASD[[6]]),
						    names(hub_lncRNA_celltype_ASD[[7]]),
						    names(hub_lncRNA_celltype_ASD[[8]]),
						    names(hub_lncRNA_celltype_ASD[[9]]),
						    names(hub_lncRNA_celltype_ASD[[10]]),
						    names(hub_lncRNA_celltype_ASD[[11]]),
						    names(hub_lncRNA_celltype_ASD[[12]]),
						    names(hub_lncRNA_celltype_ASD[[13]]),
						    names(hub_lncRNA_celltype_ASD[[14]]),
						    names(hub_lncRNA_celltype_ASD[[15]]),
						    names(hub_lncRNA_celltype_ASD[[16]]),
						    names(hub_lncRNA_celltype_ASD[[17]])), 
						    overlap.num = 2, type = "least")
biomarker_celltype_ASD_2_degree <- lapply(seq(celltype_ASD_graph), function(i) celltype_ASD_degree[[i]][which(names(celltype_ASD_degree[[i]]) %in% biomarker_celltype_ASD_2)])
biomarker_celltype_ASD_2_degree_mat <- do.call(cbind, biomarker_celltype_ASD_2_degree)
biomarker_celltype_ASD_2_fold <- apply(biomarker_celltype_ASD_2_degree_mat, 1, function(x) {
  max_val <- max(x)  
  max_idx <- which.max(x)  
  other_vals <- x[-max_idx]  
  mean_other <- mean(other_vals)  
  max_val / mean_other  
  }
)
biomarker_celltype_ASD_enhanced <- names(which(biomarker_celltype_ASD_2_fold > 2))

celltype_Normal_graph <- list(NeuNRGNII_Normal_res_darkcausality_priori_graph,
                       L56_Normal_res_darkcausality_priori_graph,
                       Oligodendrocytes_Normal_res_darkcausality_priori_graph,
                       OPC_Normal_res_darkcausality_priori_graph,
                       ASTFB_Normal_res_darkcausality_priori_graph,
                       Endothelial_Normal_res_darkcausality_priori_graph,
                       Microglia_Normal_res_darkcausality_priori_graph,
                       NeuNRGNI_Normal_res_darkcausality_priori_graph,
                       INVIP_Normal_res_darkcausality_priori_graph,                      
                       L56CC_Normal_res_darkcausality_priori_graph,                       
                       INSV2C_Normal_res_darkcausality_priori_graph,                       
                       L23_Normal_res_darkcausality_priori_graph,                       
                       INPV_Normal_res_darkcausality_priori_graph,                       
                       L4_Normal_res_darkcausality_priori_graph,                       
                       INSST_Normal_res_darkcausality_priori_graph,                       
                       Neumat_Normal_res_darkcausality_priori_graph,                       
                       ASTPP_Normal_res_darkcausality_priori_graph)
celltype_Normal_degree <- lapply(seq(celltype_Normal_graph), function(i) degree(celltype_Normal_graph[[i]]))
n_celltype_Normal <- lapply(seq(celltype_Normal_graph), function(i) length(celltype_Normal_degree[[i]]))
m_celltype_Normal <- lapply(seq(celltype_Normal_graph), function(i) ecount(celltype_Normal_graph[[i]]))
p1_celltype_Normal <- lapply(seq(celltype_Normal_graph), function(i) m_celltype_Normal[[i]]/(n_celltype_Normal[[i]]^2))
lamda_celltype_Normal <- lapply(seq(celltype_Normal_graph), function(i) n_celltype_Normal[[i]]*p1_celltype_Normal[[i]])
p_celltype_Normal <- lapply(seq(celltype_Normal_graph), function(i) ppois(celltype_Normal_degree[[i]] - 1, lambda = lamda_celltype_Normal[[i]], lower=FALSE))
hub_celltype_Normal <- lapply(seq(celltype_Normal_graph), function(i) p_celltype_Normal[[i]][which(p_celltype_Normal[[i]] < 0.05)])
hub_lncRNA_celltype_Normal <- lapply(seq(celltype_Normal_graph), function(i) hub_celltype_Normal[[i]][which(names(hub_celltype_Normal[[i]]) %in% colnames(ASD_model_data)[cause])])
biomarker_celltype_Normal_specific <- Overlap.hub(list(names(hub_lncRNA_celltype_Normal[[1]]), 
                                                       names(hub_lncRNA_celltype_Normal[[2]]),
						       names(hub_lncRNA_celltype_Normal[[3]]),
						       names(hub_lncRNA_celltype_Normal[[4]]),
						       names(hub_lncRNA_celltype_Normal[[5]]),
						       names(hub_lncRNA_celltype_Normal[[6]]),
						       names(hub_lncRNA_celltype_Normal[[7]]),
						       names(hub_lncRNA_celltype_Normal[[8]]),
						       names(hub_lncRNA_celltype_Normal[[9]]),
						       names(hub_lncRNA_celltype_Normal[[10]]),
						       names(hub_lncRNA_celltype_Normal[[11]]),
						       names(hub_lncRNA_celltype_Normal[[12]]),
						       names(hub_lncRNA_celltype_Normal[[13]]),
						       names(hub_lncRNA_celltype_Normal[[14]]),
						       names(hub_lncRNA_celltype_Normal[[15]]),
						       names(hub_lncRNA_celltype_Normal[[16]]),
						       names(hub_lncRNA_celltype_Normal[[17]])), 
						       overlap.num = 1, type = "equal")
biomarker_celltype_Normal_2 <- Overlap.hub(list(names(hub_lncRNA_celltype_Normal[[1]]), 
                                                    names(hub_lncRNA_celltype_Normal[[2]]),
						    names(hub_lncRNA_celltype_Normal[[3]]),
						    names(hub_lncRNA_celltype_Normal[[4]]),
						    names(hub_lncRNA_celltype_Normal[[5]]),
						    names(hub_lncRNA_celltype_Normal[[6]]),
						    names(hub_lncRNA_celltype_Normal[[7]]),
						    names(hub_lncRNA_celltype_Normal[[8]]),
						    names(hub_lncRNA_celltype_Normal[[9]]),
						    names(hub_lncRNA_celltype_Normal[[10]]),
						    names(hub_lncRNA_celltype_Normal[[11]]),
						    names(hub_lncRNA_celltype_Normal[[12]]),
						    names(hub_lncRNA_celltype_Normal[[13]]),
						    names(hub_lncRNA_celltype_Normal[[14]]),
						    names(hub_lncRNA_celltype_Normal[[15]]),
						    names(hub_lncRNA_celltype_Normal[[16]]),
						    names(hub_lncRNA_celltype_Normal[[17]])), 
						    overlap.num = 2, type = "least")
biomarker_celltype_Normal_2_degree <- lapply(seq(celltype_Normal_graph), function(i) celltype_Normal_degree[[i]][which(names(celltype_Normal_degree[[i]]) %in% biomarker_celltype_Normal_2)])
biomarker_celltype_Normal_2_degree_mat <- do.call(cbind, biomarker_celltype_Normal_2_degree)
biomarker_celltype_Normal_2_fold <- apply(biomarker_celltype_Normal_2_degree_mat, 1, function(x) {
  max_val <- max(x)  
  max_idx <- which.max(x)  
  other_vals <- x[-max_idx]  
  mean_other <- mean(other_vals)  
  max_val / mean_other  
  }
)
biomarker_celltype_Normal_enhanced <- names(which(biomarker_celltype_Normal_2_fold > 2))

# Set random seed for reproducibility
set.seed(123)

# Diagnostic
n_ASD <- ncol(ASD_lncRNAs_data)
n_Normal <- ncol(Control_lncRNAs_data)
labels_diagnostic <- factor(c(rep("ASD", n_ASD), rep("Normal", n_Normal)))
index_diagnostic_specific <- which(rownames(ASD_lncRNAs_data) %in% biomarker_diagnostic_specific)
index_diagnostic_enhanced <- which(rownames(ASD_lncRNAs_data) %in% biomarker_diagnostic_enhanced)
data_biomarker_diagnostic_specific <- rbind(t(ASD_lncRNAs_data[index_diagnostic_specific, ]), t(Control_lncRNAs_data[index_diagnostic_specific, ]))
data_biomarker_diagnostic_enhanced <- rbind(t(ASD_lncRNAs_data[index_diagnostic_enhanced, ]), t(Control_lncRNAs_data[index_diagnostic_enhanced, ]))
data_biomarker_diagnostic_specific <- data.frame(data_biomarker_diagnostic_specific, Class = labels_diagnostic)
data_biomarker_diagnostic_enhanced <- data.frame(data_biomarker_diagnostic_enhanced, Class = labels_diagnostic)
biomarker_diagnostic_specific_cv <- svm_cv(data_biomarker_diagnostic_specific, n_folds = 10, parallel = TRUE)
biomarker_diagnostic_specific_performance <- analyze_cv_results(biomarker_diagnostic_specific_cv, "binary")
biomarker_diagnostic_enhanced_cv <- svm_cv(data_biomarker_diagnostic_enhanced, n_folds = 10, parallel = TRUE)
biomarker_diagnostic_enhanced_performance <- analyze_cv_results(biomarker_diagnostic_enhanced_cv, "binary")
biomarker_diagnostic_cv <- svm_cv(cbind(data_biomarker_diagnostic_specific, data_biomarker_diagnostic_enhanced), n_folds = 10, parallel = TRUE)
biomarker_diagnostic_performance <- analyze_cv_results(biomarker_diagnostic_cv, "binary")

# Region
n_ACC_ASD <- ncol(ACC_ASD_lncRNAs_data)
n_PFC_ASD <- ncol(PFC_ASD_lncRNAs_data)
labels_region_ASD <- factor(c(rep("ACC", n_ACC_ASD), rep("PFC", n_PFC_ASD)))
index_region_ASD_specific <- which(rownames(ACC_ASD_lncRNAs_data) %in% biomarker_region_ASD_specific)
index_region_ASD_enhanced <- which(rownames(ACC_ASD_lncRNAs_data) %in% biomarker_region_ASD_enhanced)
data_biomarker_region_ASD_specific <- rbind(t(ACC_ASD_lncRNAs_data[index_region_ASD_specific, ]), t(PFC_ASD_lncRNAs_data[index_region_ASD_specific, ]))
data_biomarker_region_ASD_enhanced <- rbind(t(ACC_ASD_lncRNAs_data[index_region_ASD_enhanced, ]), t(PFC_ASD_lncRNAs_data[index_region_ASD_enhanced, ]))
data_biomarker_region_ASD_specific <- data.frame(data_biomarker_region_ASD_specific, Class = labels_region_ASD)
data_biomarker_region_ASD_enhanced <- data.frame(data_biomarker_region_ASD_enhanced, Class = labels_region_ASD)
biomarker_region_ASD_specific_cv <- svm_cv(data_biomarker_region_ASD_specific, n_folds = 10, parallel = TRUE)
biomarker_region_ASD_specific_performance <- analyze_cv_results(biomarker_region_ASD_specific_cv, "binary")
biomarker_region_ASD_enhanced_cv <- svm_cv(data_biomarker_region_ASD_enhanced, n_folds = 10, parallel = TRUE)
biomarker_region_ASD_enhanced_performance <- analyze_cv_results(biomarker_region_ASD_enhanced_cv, "binary")
biomarker_region_ASD_cv <- svm_cv(cbind(data_biomarker_region_ASD_specific, data_biomarker_region_ASD_enhanced), n_folds = 10, parallel = TRUE)
biomarker_region_ASD_performance <- analyze_cv_results(biomarker_region_ASD_cv, "binary")

n_ACC_Normal <- ncol(ACC_control_lncRNAs_data)
n_PFC_Normal <- ncol(PFC_control_lncRNAs_data)
labels_region_Normal <- factor(c(rep("ACC", n_ACC_Normal), rep("PFC", n_PFC_Normal)))
index_region_Normal_specific <- which(rownames(ACC_control_lncRNAs_data) %in% biomarker_region_Normal_specific)
index_region_Normal_enhanced <- which(rownames(ACC_control_lncRNAs_data) %in% biomarker_region_Normal_enhanced)
data_biomarker_region_Normal_specific <- rbind(t(ACC_control_lncRNAs_data[index_region_Normal_specific, ]), t(PFC_control_lncRNAs_data[index_region_Normal_specific, ]))
data_biomarker_region_Normal_enhanced <- rbind(t(ACC_control_lncRNAs_data[index_region_Normal_enhanced, ]), t(PFC_control_lncRNAs_data[index_region_Normal_enhanced, ]))
data_biomarker_region_Normal_specific <- data.frame(data_biomarker_region_Normal_specific, Class = labels_region_Normal)
data_biomarker_region_Normal_enhanced <- data.frame(data_biomarker_region_Normal_enhanced, Class = labels_region_Normal)
biomarker_region_Normal_specific_cv <- svm_cv(data_biomarker_region_Normal_specific, n_folds = 10, parallel = TRUE)
biomarker_region_Normal_specific_performance <- analyze_cv_results(biomarker_region_Normal_specific_cv, "binary")
biomarker_region_Normal_enhanced_cv <- svm_cv(data_biomarker_region_Normal_enhanced, n_folds = 10, parallel = TRUE)
biomarker_region_Normal_enhanced_performance <- analyze_cv_results(biomarker_region_Normal_enhanced_cv, "binary")
biomarker_region_Normal_cv <- svm_cv(cbind(data_biomarker_region_Normal_specific, data_biomarker_region_Normal_enhanced), n_folds = 10, parallel = TRUE)
biomarker_region_Normal_performance <- analyze_cv_results(biomarker_region_Normal_cv, "binary")

# Age
n_Lower18_ASD <- ncol(Lower18_ASD_lncRNAs_data)
n_Larger18_ASD <- ncol(Large18_ASD_lncRNAs_data)
labels_age_ASD <- factor(c(rep("Lower18", n_Lower18_ASD), rep("Larger18", n_Larger18_ASD)))
index_age_ASD_specific <- which(rownames(Lower18_ASD_lncRNAs_data) %in% biomarker_age_ASD_specific)
index_age_ASD_enhanced <- which(rownames(Lower18_ASD_lncRNAs_data) %in% biomarker_age_ASD_enhanced)
data_biomarker_age_ASD_specific <- rbind(t(Lower18_ASD_lncRNAs_data[index_age_ASD_specific, ]), t(Large18_ASD_lncRNAs_data[index_age_ASD_specific, ]))
data_biomarker_age_ASD_enhanced <- rbind(t(Lower18_ASD_lncRNAs_data[index_age_ASD_enhanced, ]), t(Large18_ASD_lncRNAs_data[index_age_ASD_enhanced, ]))
data_biomarker_age_ASD_specific <- data.frame(data_biomarker_age_ASD_specific, Class = labels_age_ASD)
data_biomarker_age_ASD_enhanced <- data.frame(data_biomarker_age_ASD_enhanced, Class = labels_age_ASD)
biomarker_age_ASD_specific_cv <- svm_cv(data_biomarker_age_ASD_specific, n_folds = 10, parallel = TRUE)
biomarker_age_ASD_specific_performance <- analyze_cv_results(biomarker_age_ASD_specific_cv, "binary")
biomarker_age_ASD_enhanced_cv <- svm_cv(data_biomarker_age_ASD_enhanced, n_folds = 10, parallel = TRUE)
biomarker_age_ASD_enhanced_performance <- analyze_cv_results(biomarker_age_ASD_enhanced_cv, "binary")
biomarker_age_ASD_cv <- svm_cv(cbind(data_biomarker_age_ASD_specific, data_biomarker_age_ASD_enhanced), n_folds = 10, parallel = TRUE)
biomarker_age_ASD_performance <- analyze_cv_results(biomarker_age_ASD_cv, "binary")

n_Lower18_Normal <- ncol(Lower18_control_lncRNAs_data)
n_Larger18_Normal <- ncol(Large18_control_lncRNAs_data)
labels_age_Normal <- factor(c(rep("Lower18", n_Lower18_Normal), rep("Larger18", n_Larger18_Normal)))
index_age_Normal_specific <- which(rownames(Lower18_control_lncRNAs_data) %in% biomarker_age_Normal_specific)
index_age_Normal_enhanced <- which(rownames(Lower18_control_lncRNAs_data) %in% biomarker_age_Normal_enhanced)
data_biomarker_age_Normal_specific <- rbind(t(Lower18_control_lncRNAs_data[index_age_Normal_specific, ]), t(Large18_control_lncRNAs_data[index_age_Normal_specific, ]))
data_biomarker_age_Normal_enhanced <- rbind(t(Lower18_control_lncRNAs_data[index_age_Normal_enhanced, ]), t(Large18_control_lncRNAs_data[index_age_Normal_enhanced, ]))
data_biomarker_age_Normal_specific <- data.frame(data_biomarker_age_Normal_specific, Class = labels_age_Normal)
data_biomarker_age_Normal_enhanced <- data.frame(data_biomarker_age_Normal_enhanced, Class = labels_age_Normal)
biomarker_age_Normal_specific_cv <- svm_cv(data_biomarker_age_Normal_specific, n_folds = 10, parallel = TRUE)
biomarker_age_Normal_specific_performance <- analyze_cv_results(biomarker_age_Normal_specific_cv, "binary")
biomarker_age_Normal_enhanced_cv <- svm_cv(data_biomarker_age_Normal_enhanced, n_folds = 10, parallel = TRUE)
biomarker_age_Normal_enhanced_performance <- analyze_cv_results(biomarker_age_Normal_enhanced_cv, "binary")
biomarker_age_Normal_cv <- svm_cv(cbind(data_biomarker_age_Normal_specific, data_biomarker_age_Normal_enhanced), n_folds = 10, parallel = TRUE)
biomarker_age_Normal_performance <- analyze_cv_results(biomarker_age_Normal_cv, "binary")

# Sex
n_Male_ASD <- ncol(Male_ASD_lncRNAs_data)
n_Female_ASD <- ncol(Female_ASD_lncRNAs_data)
labels_sex_ASD <- factor(c(rep("Male", n_Male_ASD), rep("Female", n_Female_ASD)))
index_sex_ASD_specific <- which(rownames(Male_ASD_lncRNAs_data) %in% biomarker_sex_ASD_specific)
index_sex_ASD_enhanced <- which(rownames(Male_ASD_lncRNAs_data) %in% biomarker_sex_ASD_enhanced)
data_biomarker_sex_ASD_specific <- rbind(t(Male_ASD_lncRNAs_data[index_sex_ASD_specific, ]), t(Female_ASD_lncRNAs_data[index_sex_ASD_specific, ]))
data_biomarker_sex_ASD_enhanced <- rbind(t(Male_ASD_lncRNAs_data[index_sex_ASD_enhanced, ]), t(Female_ASD_lncRNAs_data[index_sex_ASD_enhanced, ]))
data_biomarker_sex_ASD_specific <- data.frame(data_biomarker_sex_ASD_specific, Class = labels_sex_ASD)
data_biomarker_sex_ASD_enhanced <- data.frame(data_biomarker_sex_ASD_enhanced, Class = labels_sex_ASD)
biomarker_sex_ASD_specific_cv <- svm_cv(data_biomarker_sex_ASD_specific, n_folds = 10, parallel = TRUE)
biomarker_sex_ASD_specific_performance <- analyze_cv_results(biomarker_sex_ASD_specific_cv, "binary")
biomarker_sex_ASD_enhanced_cv <- svm_cv(data_biomarker_sex_ASD_enhanced, n_folds = 10, parallel = TRUE)
biomarker_sex_ASD_enhanced_performance <- analyze_cv_results(biomarker_sex_ASD_enhanced_cv, "binary")
biomarker_sex_ASD_cv <- svm_cv(cbind(data_biomarker_sex_ASD_specific, data_biomarker_sex_ASD_enhanced), n_folds = 10, parallel = TRUE)
biomarker_sex_ASD_performance <- analyze_cv_results(biomarker_sex_ASD_cv, "binary")

n_Male_Normal <- ncol(Male_control_lncRNAs_data)
n_Female_Normal <- ncol(Female_control_lncRNAs_data)
labels_sex_Normal <- factor(c(rep("Male", n_Male_Normal), rep("Female", n_Female_Normal)))
index_sex_Normal_specific <- which(rownames(Male_control_lncRNAs_data) %in% biomarker_sex_Normal_specific)
index_sex_Normal_enhanced <- which(rownames(Male_control_lncRNAs_data) %in% biomarker_sex_Normal_enhanced)
data_biomarker_sex_Normal_specific <- rbind(t(Male_control_lncRNAs_data[index_sex_Normal_specific, ]), t(Female_control_lncRNAs_data[index_sex_Normal_specific, ]))
data_biomarker_sex_Normal_enhanced <- rbind(t(Male_control_lncRNAs_data[index_sex_Normal_enhanced, ]), t(Female_control_lncRNAs_data[index_sex_Normal_enhanced, ]))
data_biomarker_sex_Normal_specific <- data.frame(data_biomarker_sex_Normal_specific, Class = labels_sex_Normal)
data_biomarker_sex_Normal_enhanced <- data.frame(data_biomarker_sex_Normal_enhanced, Class = labels_sex_Normal)
biomarker_sex_Normal_specific_cv <- svm_cv(data_biomarker_sex_Normal_specific, n_folds = 10, parallel = TRUE)
biomarker_sex_Normal_specific_performance <- analyze_cv_results(biomarker_sex_Normal_specific_cv, "binary")
biomarker_sex_Normal_enhanced_cv <- svm_cv(data_biomarker_sex_Normal_enhanced, n_folds = 10, parallel = TRUE)
biomarker_sex_Normal_enhanced_performance <- analyze_cv_results(biomarker_sex_Normal_enhanced_cv, "binary")
biomarker_sex_Normal_cv <- svm_cv(cbind(data_biomarker_sex_Normal_specific, data_biomarker_sex_Normal_enhanced), n_folds = 10, parallel = TRUE)
biomarker_sex_Normal_performance <- analyze_cv_results(biomarker_sex_Normal_cv, "binary")

# Cell type
n_celltype_ASD <- lapply(seq(ASD_celltype_lncRNAs_list), function(i) ncol(ASD_celltype_lncRNAs_list[[i]]))
labels_celltype_ASD <- unlist(lapply(seq(ASD_celltype_lncRNAs_list), function(i) factor(rep(names(ASD_celltype_lncRNAs_list)[i], n_celltype_ASD[[i]]))))
index_celltype_ASD_specific <- which(rownames(ASD_celltype_lncRNAs_list[[1]]) %in% biomarker_celltype_ASD_specific)
index_celltype_ASD_enhanced <- which(rownames(ASD_celltype_lncRNAs_list[[1]]) %in% biomarker_celltype_ASD_enhanced)
data_biomarker_celltype_ASD_specific <- do.call(rbind, lapply(seq(ASD_celltype_lncRNAs_list), function(i) t(ASD_celltype_lncRNAs_list[[i]][index_celltype_ASD_specific, ])))
data_biomarker_celltype_ASD_enhanced <- do.call(rbind, lapply(seq(ASD_celltype_lncRNAs_list), function(i) t(ASD_celltype_lncRNAs_list[[i]][index_celltype_ASD_enhanced, ])))
data_biomarker_celltype_ASD_specific <- data.frame(data_biomarker_celltype_ASD_specific, Class = labels_celltype_ASD)
data_biomarker_celltype_ASD_enhanced <- data.frame(data_biomarker_celltype_ASD_enhanced, Class = labels_celltype_ASD)
biomarker_celltype_ASD_specific_cv <- svm_cv(data_biomarker_celltype_ASD_specific, n_folds = 10, parallel = TRUE)
biomarker_celltype_ASD_specific_performance <- analyze_cv_results(biomarker_celltype_ASD_specific_cv, "multiclass")
biomarker_celltype_ASD_enhanced_cv <- svm_cv(data_biomarker_celltype_ASD_enhanced, n_folds = 10, parallel = TRUE)
biomarker_celltype_ASD_enhanced_performance <- analyze_cv_results(biomarker_celltype_ASD_enhanced_cv, "multiclass")
biomarker_celltype_ASD_cv <- svm_cv(cbind(data_biomarker_celltype_ASD_specific, data_biomarker_celltype_ASD_enhanced), n_folds = 10, parallel = TRUE)
biomarker_celltype_ASD_performance <- analyze_cv_results(biomarker_celltype_ASD_cv, "multiclass")

n_celltype_Normal <- lapply(seq(Control_celltype_lncRNAs_list), function(i) ncol(Control_celltype_lncRNAs_list[[i]]))
labels_celltype_Normal <- unlist(lapply(seq(Control_celltype_lncRNAs_list), function(i) factor(rep(names(Control_celltype_lncRNAs_list)[i], n_celltype_Normal[[i]]))))
index_celltype_Normal_specific <- which(rownames(Control_celltype_lncRNAs_list[[1]]) %in% biomarker_celltype_Normal_specific)
index_celltype_Normal_enhanced <- which(rownames(Control_celltype_lncRNAs_list[[1]]) %in% biomarker_celltype_Normal_enhanced)
data_biomarker_celltype_Normal_specific <- do.call(rbind, lapply(seq(Control_celltype_lncRNAs_list), function(i) t(Control_celltype_lncRNAs_list[[i]][index_celltype_Normal_specific, ])))
data_biomarker_celltype_Normal_enhanced <- do.call(rbind, lapply(seq(Control_celltype_lncRNAs_list), function(i) t(Control_celltype_lncRNAs_list[[i]][index_celltype_Normal_enhanced, ])))
data_biomarker_celltype_Normal_specific <- data.frame(data_biomarker_celltype_Normal_specific, Class = labels_celltype_Normal)
data_biomarker_celltype_Normal_enhanced <- data.frame(data_biomarker_celltype_Normal_enhanced, Class = labels_celltype_Normal)
biomarker_celltype_Normal_specific_cv <- svm_cv(data_biomarker_celltype_Normal_specific, n_folds = 10, parallel = TRUE)
biomarker_celltype_Normal_specific_performance <- analyze_cv_results(biomarker_celltype_Normal_specific_cv, "multiclass")
biomarker_celltype_Normal_enhanced_cv <- svm_cv(data_biomarker_celltype_Normal_enhanced, n_folds = 10, parallel = TRUE)
biomarker_celltype_Normal_enhanced_performance <- analyze_cv_results(biomarker_celltype_Normal_enhanced_cv, "multiclass")
biomarker_celltype_Normal_cv <- svm_cv(cbind(data_biomarker_celltype_Normal_specific, data_biomarker_celltype_Normal_enhanced), n_folds = 10, parallel = TRUE)
biomarker_celltype_Normal_performance <- analyze_cv_results(biomarker_celltype_Normal_cv, "multiclass")

############################# 5. Immune-related lncRNAs detection (Immune-specific and Immune-enhanced lncRNAs) and enrichment analysis #############################
# Microglia immune cells
Microglia_graph <- list(Microglia_ASD_res_darkcausality_priori_graph, Microglia_Normal_res_darkcausality_priori_graph)
Microglia_degree <- lapply(seq(Microglia_graph), function(i) degree(Microglia_graph[[i]]))
n_Microglia <- lapply(seq(Microglia_graph), function(i) length(Microglia_degree[[i]]))
m_Microglia <- lapply(seq(Microglia_graph), function(i) ecount(Microglia_graph[[i]]))
p1_Microglia <- lapply(seq(Microglia_graph), function(i) m_Microglia[[i]]/(n_Microglia[[i]]^2))
lamda_Microglia <- lapply(seq(Microglia_graph), function(i) n_Microglia[[i]]*p1_Microglia[[i]])
p_Microglia <- lapply(seq(Microglia_graph), function(i) ppois(Microglia_degree[[i]] - 1, lambda = lamda_Microglia[[i]], lower=FALSE))
hub_Microglia <- lapply(seq(Microglia_graph), function(i) p_Microglia[[i]][which(p_Microglia[[i]] < 0.05)])
hub_lncRNA_Microglia <- lapply(seq(Microglia_graph), function(i) hub_Microglia[[i]][which(names(hub_Microglia[[i]]) %in% colnames(ASD_model_data)[cause])])
immune_Microglia_specific <- Overlap.hub(list(names(hub_lncRNA_Microglia[[1]]), names(hub_lncRNA_Microglia[[2]])), overlap.num = 1, type = "equal")
immune_Microglia_2 <- Overlap.hub(list(names(hub_lncRNA_Microglia[[1]]), names(hub_lncRNA_Microglia[[2]])), overlap.num = 2, type = "equal")
immune_Microglia_2_fold <- Microglia_degree[[1]][which(names(Microglia_degree[[1]]) %in% immune_Microglia_2)]/Microglia_degree[[2]][which(names(Microglia_degree[[2]]) %in% immune_Microglia_2)]
immune_Microglia_enhanced <- names(which(immune_Microglia_2_fold > 2 | immune_Microglia_2_fold < 0.5))

Dif_Microglia_res_darkcausality_priori_graph <- (Microglia_ASD_res_darkcausality_priori_graph %m% Microglia_Normal_res_darkcausality_priori_graph) %u% (Microglia_Normal_res_darkcausality_priori_graph %m% Microglia_ASD_res_darkcausality_priori_graph)
Dif_Microglia_lncRNA_list <- unique(tail_of(Dif_Microglia_res_darkcausality_priori_graph, E(Dif_Microglia_res_darkcausality_priori_graph))$name)
Dif_Microglia_mRNA_list <- lapply(seq(Dif_Microglia_lncRNA_list), function(i) neighbors(Dif_Microglia_res_darkcausality_priori_graph, Dif_Microglia_lncRNA_list[i])$name)
immune_Microglia_specific_lncRNAs_targets <- lapply(seq(immune_Microglia_specific), function(i) Dif_Microglia_mRNA_list[[which(Dif_Microglia_lncRNA_list %in% immune_Microglia_specific[i])]])

Comm_Microglia_res_darkcausality_priori_graph <- Microglia_ASD_res_darkcausality_priori_graph %s% Microglia_Normal_res_darkcausality_priori_graph
Comm_Microglia_lncRNA_list <- unique(tail_of(Comm_Microglia_res_darkcausality_priori_graph, E(Comm_Microglia_res_darkcausality_priori_graph))$name)
Comm_Microglia_mRNA_list <- lapply(seq(Comm_Microglia_lncRNA_list), function(i) neighbors(Comm_Microglia_res_darkcausality_priori_graph, Comm_Microglia_lncRNA_list[i])$name)
immune_Microglia_enhanced_lncRNAs_targets <- lapply(seq(immune_Microglia_enhanced), function(i) Comm_Microglia_mRNA_list[[which(Comm_Microglia_lncRNA_list %in% immune_Microglia_enhanced[i])]])

immune_Microglia_specific_enrichment_FEA_BP <- module_FA(immune_Microglia_specific_lncRNAs_targets, GOont = "BP", Analysis.type = "FEA")
immune_Microglia_specific_enrichment_FEA_MF <- module_FA(immune_Microglia_specific_lncRNAs_targets, GOont = "MF", Analysis.type = "FEA")
immune_Microglia_specific_enrichment_FEA_CC <- module_FA(immune_Microglia_specific_lncRNAs_targets, GOont = "CC", Analysis.type = "FEA")
immune_Microglia_enhanced_enrichment_FEA_BP <- module_FA(immune_Microglia_enhanced_lncRNAs_targets, GOont = "BP", Analysis.type = "FEA")
immune_Microglia_enhanced_enrichment_FEA_MF <- module_FA(immune_Microglia_enhanced_lncRNAs_targets, GOont = "MF", Analysis.type = "FEA")
immune_Microglia_enhanced_enrichment_FEA_CC <- module_FA(immune_Microglia_enhanced_lncRNAs_targets, GOont = "CC", Analysis.type = "FEA")

immune_Microglia_specific_enrichment_FEA_BP_dataframe <- lapply(seq(immune_Microglia_specific_enrichment_FEA_BP[[1]]), function(i) as.data.frame(immune_Microglia_specific_enrichment_FEA_BP[[1]][[i]]))
immune_Microglia_specific_enrichment_FEA_BP_num <- lapply(seq(immune_Microglia_specific_enrichment_FEA_BP_dataframe), function(i) nrow(immune_Microglia_specific_enrichment_FEA_BP_dataframe[[i]]))
immune_Microglia_specific_enrichment_FEA_BP_res <- do.call(rbind, immune_Microglia_specific_enrichment_FEA_BP_dataframe)

immune_Microglia_specific_enrichment_FEA_MF_dataframe <- lapply(seq(immune_Microglia_specific_enrichment_FEA_MF[[1]]), function(i) as.data.frame(immune_Microglia_specific_enrichment_FEA_MF[[1]][[i]]))
immune_Microglia_specific_enrichment_FEA_MF_num <- lapply(seq(immune_Microglia_specific_enrichment_FEA_MF_dataframe), function(i) nrow(immune_Microglia_specific_enrichment_FEA_MF_dataframe[[i]]))
immune_Microglia_specific_enrichment_FEA_MF_res <- do.call(rbind, immune_Microglia_specific_enrichment_FEA_MF_dataframe)

immune_Microglia_specific_enrichment_FEA_CC_dataframe <- lapply(seq(immune_Microglia_specific_enrichment_FEA_CC[[1]]), function(i) as.data.frame(immune_Microglia_specific_enrichment_FEA_CC[[1]][[i]]))
immune_Microglia_specific_enrichment_FEA_CC_num <- lapply(seq(immune_Microglia_specific_enrichment_FEA_CC_dataframe), function(i) nrow(immune_Microglia_specific_enrichment_FEA_CC_dataframe[[i]]))
immune_Microglia_specific_enrichment_FEA_CC_res <- do.call(rbind, immune_Microglia_specific_enrichment_FEA_CC_dataframe)

immune_Microglia_specific_enrichment_FEA_KEGG_dataframe <- lapply(seq(immune_Microglia_specific_enrichment_FEA_BP[[2]]), function(i) as.data.frame(immune_Microglia_specific_enrichment_FEA_BP[[2]][[i]]))
immune_Microglia_specific_enrichment_FEA_KEGG_num <- lapply(seq(immune_Microglia_specific_enrichment_FEA_KEGG_dataframe), function(i) nrow(immune_Microglia_specific_enrichment_FEA_KEGG_dataframe[[i]]))
immune_Microglia_specific_enrichment_FEA_KEGG_res <- do.call(rbind, immune_Microglia_specific_enrichment_FEA_KEGG_dataframe)

immune_Microglia_specific_enrichment_FEA_Reactome_dataframe <- lapply(seq(immune_Microglia_specific_enrichment_FEA_BP[[3]]), function(i) as.data.frame(immune_Microglia_specific_enrichment_FEA_BP[[3]][[i]]))
immune_Microglia_specific_enrichment_FEA_Reactome_num <- lapply(seq(immune_Microglia_specific_enrichment_FEA_Reactome_dataframe), function(i) nrow(immune_Microglia_specific_enrichment_FEA_Reactome_dataframe[[i]]))
immune_Microglia_specific_enrichment_FEA_Reactome_res <- do.call(rbind, immune_Microglia_specific_enrichment_FEA_Reactome_dataframe)

immune_Microglia_enhanced_enrichment_FEA_BP_dataframe <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_BP[[1]]), function(i) as.data.frame(immune_Microglia_enhanced_enrichment_FEA_BP[[1]][[i]]))
immune_Microglia_enhanced_enrichment_FEA_BP_num <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_BP_dataframe), function(i) nrow(immune_Microglia_enhanced_enrichment_FEA_BP_dataframe[[i]]))
immune_Microglia_enhanced_enrichment_FEA_BP_res <- do.call(rbind, immune_Microglia_enhanced_enrichment_FEA_BP_dataframe)

immune_Microglia_enhanced_enrichment_FEA_MF_dataframe <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_MF[[1]]), function(i) as.data.frame(immune_Microglia_enhanced_enrichment_FEA_MF[[1]][[i]]))
immune_Microglia_enhanced_enrichment_FEA_MF_num <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_MF_dataframe), function(i) nrow(immune_Microglia_enhanced_enrichment_FEA_MF_dataframe[[i]]))
immune_Microglia_enhanced_enrichment_FEA_MF_res <- do.call(rbind, immune_Microglia_enhanced_enrichment_FEA_MF_dataframe)

immune_Microglia_enhanced_enrichment_FEA_CC_dataframe <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_CC[[1]]), function(i) as.data.frame(immune_Microglia_enhanced_enrichment_FEA_CC[[1]][[i]]))
immune_Microglia_enhanced_enrichment_FEA_CC_num <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_CC_dataframe), function(i) nrow(immune_Microglia_enhanced_enrichment_FEA_CC_dataframe[[i]]))
immune_Microglia_enhanced_enrichment_FEA_CC_res <- do.call(rbind, immune_Microglia_enhanced_enrichment_FEA_CC_dataframe)

immune_Microglia_enhanced_enrichment_FEA_KEGG_dataframe <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_BP[[2]]), function(i) as.data.frame(immune_Microglia_enhanced_enrichment_FEA_BP[[2]][[i]]))
immune_Microglia_enhanced_enrichment_FEA_KEGG_num <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_KEGG_dataframe), function(i) nrow(immune_Microglia_enhanced_enrichment_FEA_KEGG_dataframe[[i]]))
immune_Microglia_enhanced_enrichment_FEA_KEGG_res <- do.call(rbind, immune_Microglia_enhanced_enrichment_FEA_KEGG_dataframe)

immune_Microglia_enhanced_enrichment_FEA_Reactome_dataframe <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_BP[[3]]), function(i) as.data.frame(immune_Microglia_enhanced_enrichment_FEA_BP[[3]][[i]]))
immune_Microglia_enhanced_enrichment_FEA_Reactome_num <- lapply(seq(immune_Microglia_enhanced_enrichment_FEA_Reactome_dataframe), function(i) nrow(immune_Microglia_enhanced_enrichment_FEA_Reactome_dataframe[[i]]))
immune_Microglia_enhanced_enrichment_FEA_Reactome_res <- do.call(rbind, immune_Microglia_enhanced_enrichment_FEA_Reactome_dataframe)

############################# 6. Cell-cell interaction network construction #############################
## Calculating interaction strength based on similarity matrix
LR <- as.matrix(read.csv("LR_CellPhoneDB+CellChatDB+CellTalkDB+Cellinker+CellCall.csv", header = TRUE, sep=","))
LR_graph <- make_graph(c(t(LR[, 1:2])), directed = FALSE)

lncR_celltype_ASD_net <- list(NeuNRGNII_ASD_res_darkcausality_priori_graph,
                     L56_ASD_res_darkcausality_priori_graph,
                     Oligodendrocytes_ASD_res_darkcausality_priori_graph,
		     OPC_ASD_res_darkcausality_priori_graph,
                     ASTFB_ASD_res_darkcausality_priori_graph,
                     Endothelial_ASD_res_darkcausality_priori_graph,
                     Microglia_ASD_res_darkcausality_priori_graph,
		     NeuNRGNI_ASD_res_darkcausality_priori_graph,
		     INVIP_ASD_res_darkcausality_priori_graph,
		     L56CC_ASD_res_darkcausality_priori_graph,
		     INSV2C_ASD_res_darkcausality_priori_graph,
		     L23_ASD_res_darkcausality_priori_graph,
		     INPV_ASD_res_darkcausality_priori_graph,
		     L4_ASD_res_darkcausality_priori_graph,
		     INSST_ASD_res_darkcausality_priori_graph,
		     Neumat_ASD_res_darkcausality_priori_graph,
		     ASTPP_ASD_res_darkcausality_priori_graph
                     )

lncR_celltype_ASD_gene_list <- lapply(seq(lncR_celltype_ASD_net), function(i) unique(head_of(lncR_celltype_ASD_net[[i]], E(lncR_celltype_ASD_net[[i]]))$name))
lncR_celltype_ASD_existing_genes <- lapply(seq(lncR_celltype_ASD_net), function(i) lncR_celltype_ASD_gene_list[[i]][lncR_celltype_ASD_gene_list[[i]] %in% V(LR_graph)$name])
lncR_celltype_ASD_LR_graph_sub <- lapply(seq(lncR_celltype_ASD_net), function(i) induced_subgraph(graph = LR_graph, vids = which(V(LR_graph)$name %in% lncR_celltype_ASD_existing_genes[[i]])))
lncR_celltype_ASD_net_union <- lapply(seq(lncR_celltype_ASD_net), function(i) lncR_celltype_ASD_net[[i]] %u% lncR_celltype_ASD_LR_graph_sub[[i]])

lncR_celltype_Normal_net <- list(NeuNRGNII_Normal_res_darkcausality_priori_graph,
                     L56_Normal_res_darkcausality_priori_graph,
                     Oligodendrocytes_Normal_res_darkcausality_priori_graph,
		     OPC_Normal_res_darkcausality_priori_graph,
                     ASTFB_Normal_res_darkcausality_priori_graph,
                     Endothelial_Normal_res_darkcausality_priori_graph,
                     Microglia_Normal_res_darkcausality_priori_graph,
		     NeuNRGNI_Normal_res_darkcausality_priori_graph,
		     INVIP_Normal_res_darkcausality_priori_graph,
		     L56CC_Normal_res_darkcausality_priori_graph,
		     INSV2C_Normal_res_darkcausality_priori_graph,
		     L23_Normal_res_darkcausality_priori_graph,
		     INPV_Normal_res_darkcausality_priori_graph,
		     L4_Normal_res_darkcausality_priori_graph,
		     INSST_Normal_res_darkcausality_priori_graph,
		     Neumat_Normal_res_darkcausality_priori_graph,
		     ASTPP_Normal_res_darkcausality_priori_graph
                     )

lncR_celltype_Normal_gene_list <- lapply(seq(lncR_celltype_Normal_net), function(i) unique(head_of(lncR_celltype_Normal_net[[i]], E(lncR_celltype_Normal_net[[i]]))$name))
lncR_celltype_Normal_existing_genes <- lapply(seq(lncR_celltype_Normal_net), function(i) lncR_celltype_Normal_gene_list[[i]][lncR_celltype_Normal_gene_list[[i]] %in% V(LR_graph)$name])
lncR_celltype_Normal_LR_graph_sub <- lapply(seq(lncR_celltype_Normal_net), function(i) induced_subgraph(graph = LR_graph, vids = which(V(LR_graph)$name %in% lncR_celltype_Normal_existing_genes[[i]])))
lncR_celltype_Normal_net_union <- lapply(seq(lncR_celltype_Normal_net), function(i) lncR_celltype_Normal_net[[i]] %u% lncR_celltype_Normal_LR_graph_sub[[i]])

Str_comm_ASD <- Sim.network(lncR_celltype_ASD_net_union, lncR_celltype_ASD_net_union)
Str_comm_Normal <- Sim.network(lncR_celltype_Normal_net_union, lncR_celltype_Normal_net_union)


