library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(reticulate)
library(enrichR)
library(R.utils)
library(remotes)
library(ALRA)
library(readr)
library(rliger)
library(qs)
library(MetaNeighbor)
library(SummarizedExperiment)
library(networkD3)
library(ggraph)
library(igraph)



setwd("~/PhD/scSeq/Zach_Cam paper")

Mouse_DAneuron <- read_rds(file = "Zach_DAneuron")

#Correcting for dropouts using ALRA on Zach data
counts.1 <- as.matrix(Mouse_DAneuron[["RNA"]]$counts.1)
A_norm <- t(counts.1) #Transpose the expression matrix, because the cells need to be rows and genes need to be columns
A_norm <- normalize_data(A_norm) #Normalize the cells on the total amount of UMI and, multiply by 10000 and log transform it
result.completed <- alra(A_norm) #Perform ALRA on extracted expression matrix
A_norm_completed <- result.completed[[3]]
A_norm_completed <- t(A_norm_completed) #Transpose the ALRA output
cell_names <- rownames(A_norm)
colnames(A_norm_completed) <- cell_names #Assign the column names as cells ID
A_norm.1 <- A_norm_completed
saveRDS(A_norm.1, file = "Zachdata_ALRA_A_norm1.rds")

counts.2 <- as.matrix(Mouse_DAneuron[["RNA"]]$counts.2)
A_norm <- t(counts.2) #Transpose the expression matrix, because the cells need to be rows and genes need to be columns
A_norm <- normalize_data(A_norm) #Normalize the cells on the total amount of UMI and, multiply by 10000 and log transform it
result.completed <- alra(A_norm) #Perform ALRA on extracted expression matrix
A_norm_completed <- result.completed[[3]]
A_norm_completed <- t(A_norm_completed) #Transpose the ALRA output
cell_names <- rownames(A_norm)
colnames(A_norm_completed) <- cell_names #Assign the column names as cells ID
A_norm.2 <- A_norm_completed
saveRDS(A_norm.2, file = "Zachdata_ALRA_A_norm2.rds")

counts.3 <- as.matrix(Mouse_DAneuron[["RNA"]]$counts.3)
A_norm <- t(counts.3) #Transpose the expression matrix, because the cells need to be rows and genes need to be columns
A_norm <- normalize_data(A_norm) #Normalize the cells on the total amount of UMI and, multiply by 10000 and log transform it
result.completed <- alra(A_norm) #Perform ALRA on extracted expression matrix
A_norm_completed <- result.completed[[3]]
A_norm_completed <- t(A_norm_completed) #Transpose the ALRA output
cell_names <- rownames(A_norm)
colnames(A_norm_completed) <- cell_names #Assign the column names as cells ID
A_norm.3 <- A_norm_completed
saveRDS(A_norm.3, file = "Zachdata_ALRA_A_norm3.rds")

counts.4 <- as.matrix(Mouse_DAneuron[["RNA"]]$counts.4)
A_norm <- t(counts.4) #Transpose the expression matrix, because the cells need to be rows and genes need to be columns
A_norm <- normalize_data(A_norm) #Normalize the cells on the total amount of UMI and, multiply by 10000 and log transform it
result.completed <- alra(A_norm) #Perform ALRA on extracted expression matrix
A_norm_completed <- result.completed[[3]]
A_norm_completed <- t(A_norm_completed) #Transpose the ALRA output
cell_names <- rownames(A_norm)
colnames(A_norm_completed) <- cell_names #Assign the column names as cells ID
A_norm.4 <- A_norm_completed
saveRDS(A_norm.4, file = "Zachdata_ALRA_A_norm4.rds")

A_norm_completed.1 <- cbind(A_norm.1, A_norm.2) 
A_norm_completed.2 <- cbind(A_norm.3, A_norm.4)
A_norm_completed <- cbind(A_norm_completed.1, A_norm_completed.2)
saveRDS(A_norm_completed, file = "Zachdata_ALRA_A_norm_complet.rds")

A_norm_completed <- readRDS(file = "Zachdata_ALRA_A_norm_complet.rds")

Mouse_DAneuron$ALRA <- CreateAssayObject(counts = A_norm_completed) #Add the ALRA expression matrix to a new assay in the Seurat object

saveRDS(Mouse_DAneuron, file = "Zach_DAneuron_ALRA.rds")
