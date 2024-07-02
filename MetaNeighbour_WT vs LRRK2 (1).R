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
library(ggpubr)




setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/scSeq/Zach_Cam paper")

Mouse_DAneuron <- read_rds(file = "Zach_DAneuron_ALRA.rds")

#Exclude the clusters 16-21
Mouse_DAneuron <- subset(Mouse_DAneuron, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20"))

#Subset the control samples coming from control and lrrk2 samples
control_DAneuron <- subset(Mouse_DAneuron, subset = group == "control")
lrrk2_DAneuron <- subset(Mouse_DAneuron, subset = group == "lrrk2")

#------------------------------------------------------------------------------------------------------------------
#Generate the SummarizedExperiment object recquired for MetaNeighbour
#Generate the expression matrix, where rows are the gene names and the columns are the Cell_Type-cellID
control <- GetAssayData(control_DAneuron, assay = "ALRA", layer = "data")
control <- as.matrix(control)

lrrk2 <- GetAssayData(lrrk2_DAneuron, assay = "ALRA", layer = "data")
lrrk2 <- as.matrix(lrrk2)

expression_matrix <- cbind(control, lrrk2)
row.names(expression_matrix) #Inspect if the rows are genes
colnames(expression_matrix) #Inspect if the columns are barcodes

#Generate the gene sets
Genes <- row.names(control)

#Provide the additional information
#Sample_ID
control_Sample_ID <- as.character(control_DAneuron@meta.data[["sample"]])
lrrk2_Sample_ID <- as.character(lrrk2_DAneuron@meta.data[["sample"]])
Sample_ID <- c(control_Sample_ID, lrrk2_Sample_ID)

#Study_ID
control_Study_ID <- rep("control",ncol(control_DAneuron))
lrrk2_Study_ID <- rep("lrrk2",ncol(lrrk2_DAneuron)) 
Study_ID <- c(control_Study_ID, lrrk2_Study_ID)

#Cell_Type_ID
control_Cell_Type_ID <- control_DAneuron@meta.data[["seurat_clusters"]]
lrrk2_Cell_Type_ID <- lrrk2_DAneuron@meta.data[["seurat_clusters"]]
Cell_Type_ID <- c(control_Cell_Type_ID, lrrk2_Cell_Type_ID)

#Generate the colData object, where rows are cells, and columns are Sample_ID, Study_ID and Cell_Type_ID
barcodes <- c(colnames(control_DAneuron), colnames(lrrk2_DAneuron))
df <- data.frame(Sample_ID, Study_ID, Cell_Type_ID)
row.names(df) <- colnames(expression_matrix)

# Generate the cell_labels matrix
cell_labels <- data.frame(barcodes, Cell_Type_ID)
cell_labels <- xtabs(~barcodes+Cell_Type_ID, data=cell_labels) #Generate the binary matrix


#Generate the summarizedexperiment object
SE <- SummarizedExperiment(assays = SimpleList(expression_matrix),
                     rowData = Genes, colData = DataFrame(df), metadata = cell_labels
)

var_genes = variableGenes(dat = SE, exp_labels = SE$Study_ID)
celltype_NV <- MetaNeighborUS(dat = SE,
                             var_genes = var_genes,
                             study_id = Study_ID,
                             cell_type = Cell_Type_ID,
                             fast_version = TRUE)

celltype_NV <- celltype_NV[1:20,21:40]

#Organize the clusters by similarities
celltype_NV <- celltype_NV[c("control|18", "control|14", "control|5", "control|11", "control|6", "control|12", "control|13", "control|9", "control|19", "control|7", "control|10", "control|8", "control|2", "control|4", "control|3", "control|1", "control|20", "control|0", "control|15", "control|17"),
            c("lrrk2|18", "lrrk2|14", "lrrk2|5", "lrrk2|11", "lrrk2|6", "lrrk2|12", "lrrk2|13", "lrrk2|9", "lrrk2|19", "lrrk2|7", "lrrk2|10", "lrrk2|8", "lrrk2|2", "lrrk2|4", "lrrk2|3", "lrrk2|1", "lrrk2|20", "lrrk2|0", "lrrk2|15", "lrrk2|17")]
#------------------------------------------------------------------------------------------------------------------
#Data visualization
#Heatmap
png(file='MetaNeighbour_control vs lrrk2 Heatmap without 16 and 21.png', width=1000, height=1000, pointsize=25)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(
  celltype_NV,
  keysize = 1,
  key.xlab = "AUROC",
  key.title = NULL,
  trace = "none",
  density.info = "none",
  col = cols,
  breaks = breaks,
  offsetRow = 0.1,
  offsetCol = 0.1,
  cexRow = 0.7,
  cexCol = 0.7,
  margins = c(13, 13),  # Set your desired margin values here
  Colv=NA, Rowv=NA #Remove the dendrogram
)
dev.off()


#------------------------------------------------------------------------------------------------------------------
#Boxplot of cluster similarities
reciprocal_hits <- data.frame(row.names(celltype_NV),
                              colnames(celltype_NV)[apply(celltype_NV,1,which.max)], #Return the colname for each row where the AUROC score is the highest
                              apply(celltype_NV, 1, max, na.rm=TRUE)
)
colnames(reciprocal_hits) <- c("control Cluster", "lrrk2 Cluster", "AUROC Score")
View(reciprocal_hits)
View(celltype_NV)

same_category <- diag(celltype_NV)
cross_category <- celltype_NV[!celltype_NV %in% same_category]

# Finding maximum length 
data <- as.data.frame(
c(same_category, cross_category),
c(rep("Same Category", length(same_category)), rep("Cross-Category", length(cross_category))))
colnames(data) <- c("AUROC Score")
data$Category <- c(rep("Same Category", length(same_category)), rep("Cross-Category", length(cross_category)))


png(file='Boxplot_control vs lrrk2 AUROC Score without 16 and 21.png', width=300, height=500)

ggboxplot(data, x = "Category", y = "AUROC Score",
          color = "Category", palette = "jco",
          add = "jitter", xlab = "", ylab = "AUROC Score", ylim = c(0,1.05)) + stat_compare_means(hjust = 0, vjust = -3.4, method = "t.test")

dev.off()

