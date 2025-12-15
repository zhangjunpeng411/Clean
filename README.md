# :hammer: Clean
**Hidden causal inference delineates dynamic lncRNA regulation in autism spectrum disorder**

## :boom: Background
Long non-coding RNAs (lncRNAs) are increasingly implicated in autism spectrum disorder (ASD), yet their causal regulatory mechanisms remain poorly characterized. Conventional static causal inference methods fail to capture the dynamic interplay of lncRNAs across diverse brain biological contexts. Here, we develop a hidden causality-based method **Clean** to infer dynamic lncRNA causal regulation in ASD.

A schematic illustration of **Clean** is shown in the folowing.

<p align="center">
  <img src="https://github.com/zhangjunpeng411/Clean/blob/main/Clean_schematic_illustration.png" alt="Schematic illustration of Clean" border="0.1">
</p>

According to brain biological contexts (i.e. diagnosis, region, age, sex and cell type), the snRNA-seq data of autism and normal brain cells is splitted into 48 data slices. For each data slice, metacells are built for removing technical noise and preserving the biological information of it. For the identified metacells, **Clean** uses trajectory inference to order them along a linear trajectory. Given the processed 48 data slices, **Clean** applies hidden causal inference and incorporates priori information of lncRNA targets to infer dynamic lncRNA causal regulation. Based on the identified dynamic lncRNA causal regulation, **Clean** further conducts heterogeneity analysis to capture conditional shifts in lncRNA regulation.

## :book: Description of each file in R folders
- **Clean.R**: Utility functions for delineating dynamic lncRNA regulation in ASD.

- **Case_study.R**: Case study for identifying dynamic lncRNA regulation in ASD.

## :gear: The usage of Clean
Paste all files into a single folder (set the folder as the directory of R environment). The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("R/Case_study.R")
```

## :zap: Quick example to use Clean
For inferring dynamic lncRNA regulation, users should prepare matched lncRNA and mRNA expression data and putative lncRNA targets in diffferent contexts. Paste the datasets and our source file (**Clean.R**) into a single folder (set the folder as the directory of R environment), users can use the following scripts to infer lncRNA regulation in two contexts (ASD and Normal). For convenience, our ASD single-nucleus RNA sequencing data prepared for users are from [here](https://drive.google.com/file/d/1sdLk7AcZszXYdV47_EMEVVhf0p1E_vDo/view?usp=drive_link).

```{r echo=FALSE, results='hide', message=FALSE}
## Load ASD dataset and source script
load("data/ASD_lncR_mR.RData")
source("R/Clean.R")

## Load packages
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

## Index of cause and effect for causal inference
cause <- 5134:5946
effect <- 1:5133

## Priori information of lncRNA targets
lncRTarget_priori <- as.matrix(read.csv("LncTar.csv", header = TRUE, sep=","))
lncRTarget_priori_graph <- make_graph(c(t(lncRTarget_priori[, 1:2])), directed = FALSE)

## Identifying metacells
# ASD
ASD_SC30 <- SCimplify(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 30, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
ASD_SC30_GE <- supercell_GE(as.matrix(rbind(ASD_mRNAs_data, ASD_lncRNAs_data)), ASD_SC30$membership)
colnames(ASD_SC30_GE) <- paste("MC", seq(ncol(ASD_SC30_GE)), sep = "")

# Normal
Normal_SC30 <- SCimplify(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)),  # gene expression matrix 
                      k.knn = 5, # number of nearest neighbors to build kNN network
                      gamma = 30, # graining level
                      n.var.genes = 2000 # number of the top variable genes to use for dimentionality reduction
		      )
Normal_SC30_GE <- supercell_GE(as.matrix(rbind(Control_mRNAs_data, Control_lncRNAs_data)), Normal_SC30$membership)
colnames(Normal_SC30_GE) <- paste("MC", seq(ncol(Normal_SC30_GE)), sep = "")

## Cell trajectory inference
# ASD and Normal
ASD_space <- reduce_dimensionality(t(as.matrix(ASD_SC30_GE)), ndim =2)
Normal_space <- reduce_dimensionality(t(as.matrix(Normal_SC30_GE)), ndim = 2)

ASD_model <- infer_trajectory(ASD_space)
Normal_model <- infer_trajectory(Normal_space)

ASD_model_data <- t(as.matrix(ASD_SC30_GE))[match(names(sort(ASD_model$time)), colnames(ASD_SC30_GE)), ]
Normal_model_data <- t(as.matrix(Normal_SC30_GE))[match(names(sort(Normal_model$time)), colnames(Normal_SC30_GE)), ]

## Causal inference using Clean
# ASD and Normal
ASD_res_darkcausality <- darkcausality_parallel(ASD_model_data, cause, effect, num.cores = 48)
Normal_res_darkcausality <- darkcausality_parallel(Normal_model_data, cause, effect, num.cores = 48)
```    
