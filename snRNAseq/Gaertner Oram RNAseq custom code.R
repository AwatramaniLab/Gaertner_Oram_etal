#Custom code for RNAseq analysis in Gaertner & Oram et al 2024
#Code by Zack Gaertner, 2024
#Email: zachary.gaertner@northwestern.edu


#Code for generation of Seurat object from original data is below
#follow these steps for generating a seurat object necersary for all subsequent analyses
#
#load libraries
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(patchwork)
library(dplyr)
library(ggplot2)
library(glmGamPoi)
library(SeuratWrappers)
library(ggplot2)
library(future)
#
#Read in files
control1 <- Read10X(data.dir = "path/to/your/datasets")
lrrk1 <- Read10X(data.dir = "path/to/your/datasets")
control2 <- Read10X(data.dir = "path/to/your/datasets")
lrrk2 <- Read10X(data.dir = "path/to/your/datasets")
#
#Make objects
control1 <- CreateSeuratObject(control1)
lrrk1 <- CreateSeuratObject(lrrk1)
control2 <- CreateSeuratObject(control2)
lrrk2 <- CreateSeuratObject(lrrk2)
#
#Add metadata slots
control1[['group']] <- 'control'
lrrk1[['group']] <- 'lrrk2'
control2[['group']] <- 'control'
lrrk2[['group']] <- 'lrrk2'
control1[['sample']] <- 'control1'
lrrk1[['sample']] <- 'lrrk1'
control2[['sample']] <- 'control2'
lrrk2[['sample']] <- 'lrrk2'
control1 <- PercentageFeatureSet(control1, pattern = "^mt-", col.name = "percent.mt")
lrrk1 <- PercentageFeatureSet(lrrk1, pattern = "^mt-", col.name = "percent.mt")
control2 <- PercentageFeatureSet(control2, pattern = "^mt-", col.name = "percent.mt")
lrrk2 <- PercentageFeatureSet(lrrk2, pattern = "^mt-", col.name = "percent.mt")
control1 <- PercentageFeatureSet(control1, pattern = "^Rp[sl]", col.name = "percent.ribo")
lrrk1 <- PercentageFeatureSet(lrrk1, pattern = "^Rp[sl]", col.name = "percent.ribo")
control2 <- PercentageFeatureSet(control2, pattern = "^Rp[sl]", col.name = "percent.ribo")
lrrk2 <- PercentageFeatureSet(lrrk2, pattern = "^Rp[sl]", col.name = "percent.ribo")
#
#QC filtering
control1 <- subset(control1, subset = nFeature_RNA > 1200 & nFeature_RNA < 7800 & percent.mt < 0.5 & nCount_RNA < 29000 & percent.ribo < 0.5)
lrrk1 <- subset(lrrk1, subset = nFeature_RNA > 1200 & nFeature_RNA < 6500 & percent.mt < 0.5 & nCount_RNA < 26000 & percent.ribo < 0.5)
control2 <- subset(control2, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent.mt < 0.5 & nCount_RNA < 15000 & percent.ribo < 0.5)
lrrk2 <- subset(lrrk2, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 0.5 & nCount_RNA < 20000 & percent.ribo < 0.5)
#
#V5 integration pipeline
lrrkdata.v5 <- merge(control1, list(control2, lrrk1, lrrk2))
DefaultAssay(lrrkdata.v5) <- "RNA"
lrrkdata.v5 <- SCTransform(lrrkdata.v5, vst.flavor = "v2", vars.to.regress = c('percent.mt', 'percent.ribo'))
lrrkdata.v5 <- RunPCA(lrrkdata.v5)
lrrkdata.v5 <- IntegrateLayers(object = lrrkdata.v5, method = CCAIntegration, orig.reduction = "pca", new.reduction = 'cca', assay = "SCT", normalization.method = "SCT")
lrrkdata.v5 <- IntegrateLayers(object = lrrkdata.v5, method = RPCAIntegration, orig.reduction = "pca", new.reduction = 'rpca', assay = "SCT", normalization.method = "SCT")
lrrkdata.v5 <- FindNeighbors(lrrkdata.v5, dims = 1:35, reduction = "cca")
lrrkdata.v5 <- FindClusters(lrrkdata.v5, resolution = 0.8)
lrrkdata.v5 <- RunUMAP(lrrkdata.v5, reduction = "cca", dims = 1:35, reduction.name = "umap.cca")
#Isolate EW nucleus cells
ewcells <- WhichCells(lrrkdata.v5, idents = '23')
#Reprocess dataset without EW nuc cells, custom settings to try and maximize separation of closely related cell types
lrrkdata.v5 <- merge(control1, list(control2, lrrk1, lrrk2))
lrrkdata.v5 <- subset(lrrkdata.v5, cells = ewcells, invert = T)
lrrkdata.v5 <- SCTransform(lrrkdata.v5, vst.flavor = "v2")
lrrkdata.v5 <- RunPCA(lrrkdata.v5)
lrrkdata.v5 <- IntegrateLayers(object = lrrkdata.v5, method = CCAIntegration, orig.reduction = "pca", new.reduction = 'cca', assay = "SCT", normalization.method = "SCT")
lrrkdata.v5 <- FindNeighbors(lrrkdata.v5, dims = 1:32, reduction = "cca", n.trees = 500, k.param = 40)
lrrkdata.v5 <- FindClusters(lrrkdata.v5, resolution = 0.8, algorithm = 1, group.singletons = TRUE, graph.name = "SCT_snn")
lrrkdata.v5 <- RunUMAP(lrrkdata.v5, reduction = "cca", dims = 1:32, reduction.name = "umap.cca", n.epochs = 500, min.dist = .2, n.neighbors = 1000)


#
#Code below this comment is for specific custom pipelines used in analyses in this paper ONLY
#Below code does NOT include an exhaustive reproduction of standard analysis pipelines from third party R packages, as these do not require any custom code to reproduce
#Thus, any analyses not shown below utilize standard commands in Seurat or other R packages as discussed in the source paper and methods
#For any custom arguments or parameters used in standard pipelines, please see the methods section of Gaertner & Oram et al 2024
#
#
#
#Calculation of sex proportions for each sample
df <- data.frame(lrrkdata.v5@assays$SCT$counts["Uty",])
df[1] <- lrrkdata.v5@assays$SCT$counts["Uty",]
df[2] <- lrrkdata.v5@assays$SCT$counts["Eif2s3y",]
df[3] <- lrrkdata.v5@assays$SCT$counts["Tsix",]
df[4] <- lrrkdata.v5@assays$SCT$counts["Xist",]
df[5] <- df[1] + df[2]
df[6] <- df[3] + df[4]
df[7] <- (df[5]+0.0000001)/(df[6]+0.0000001)
colnames(df) <- c('uty', 'eif','tsix', 'xist', 'male', 'female', 'ratio')
rownames(df) <- colnames(lrrkdata.v5)
hist(log2(df$'ratio'), xlab = "Normalized Sex Gene Ratio", main = "Distribution of Cells by Inferred Sex")

#Prepping for downstream DEG testing
lrrkdata.v5$cluster.group <- paste(lrrkdata.v5$seurat_clusters, lrrkdata.v5$group, sep = "_")

#Main figure dendrogram generation
Idents(lrrkdata.v5) <- "seurat_clusters"
tree.usingPCs <- BuildClusterTree(lrrkdata.v5, dims = 1:32, reduction = "cca")
PlotClusterTree(tree.usingPCs, direction = 'rightwards', show.tip.label = TRUE, show.node.label = FALSE)

#Calculating markers for dendrogram division points
for (i in 21:39) {
  # Run FindMarkers for each cluster
  markers <- FindMarkers(tree.usingPCs, ident.1 = "clustertree", ident.2 = as.character(i), 
                         min.diff.pct = 0.25, logfc.threshold = 1, min.pct = 0.2, recorrect_umi = FALSE)
  # Select top 25 with highest positive avg_log2FC
  top_positive <- markers %>%
    filter(avg_log2FC > 0) %>%
    top_n(n = 25, wt = avg_log2FC)
  # Select top 25 with lowest negative avg_log2FC
  top_negative <- markers %>%
    filter(avg_log2FC < 0) %>%
    top_n(n = 25, wt = -avg_log2FC)
  # Combine the two sets
  final_markers <- bind_rows(top_positive, top_negative)
  # Assign the result to a new variable with dynamic name
  assign(paste0("nodemarkers.", i), final_markers, envir = .GlobalEnv)
}

#Heatmap of top marker genes from each cluster
Idents(lrrkdata.v5) <- "seurat_clusters"
# Initialize an empty dataframe to store the results
results_df <- data.frame(Gene = character(),
                         Cluster = integer(),
                         Rank = integer(),
                         stringsAsFactors = FALSE)
# Loop through clusters 0 to 20, excluding 16
for (cluster in setdiff(0:20, 16)) {
  # Run the FindConservedMarkers function for the current cluster
  result <- FindConservedMarkers(lrrkdata.v5, assay = 'SCT', ident.1 = cluster, 
                                 grouping.var = 'sample', logfc.threshold = 0.15, 
                                 min.pct = 0.1, recorrect_umi = FALSE, only.pos = TRUE)
  # Extract the top 3 gene names (rownames) from the result
  top_genes <- rownames(result)[1:10]
  # Create a temporary dataframe for the current cluster
  temp_df <- data.frame(Gene = top_genes,
                        Cluster = rep(cluster, length(top_genes)),
                        Rank = 1:length(top_genes),
                        stringsAsFactors = FALSE)
  # Append the temporary dataframe to the results dataframe
  results_df <- rbind(results_df, temp_df)
}
# Print the final dataframe
print(results_df)
# Load the dplyr package
library(dplyr)
final_df <- results_df %>%
  group_by(Gene) %>%               # Group by Gene
  arrange(Rank) %>%                # Arrange by Rank within each group
  slice(1) %>%                     # Select the first row (lowest Rank) in each group
  ungroup() %>%                    # Ungroup to return to a normal dataframe
  arrange(Cluster, Rank) %>%       # Reorder the dataframe by Cluster, then by Rank
  group_by(Cluster) %>%            # Group by Cluster
  slice_head(n = 3) %>%            # Keep only the top 3 entries for each cluster
  ungroup()                        # Ungroup to return to a normal dataframe
# Print the final dataframe
print(final_df)
#Utilize normalized RNA assay for heatmap visualization
DefaultAssay(lrrkdata.v5) <- "RNA"
lrrkdata.v5.viz <- NormalizeData(lrrkdata.v5) %>% ScaleData()
DoHeatmap(subset(lrrkdata.v5.viz, downsample = 40), features = final_df$Gene, assay = 'RNA') + theme(text = element_text(size = 25)) + NoLegend()

#Calculating distinct markers for clusters within each family
#these assume you have created separate objects comprised of only the cluster family in question
#Explanation of steps: find markers using the arguments as per below, then filter out genes with high expression% in the 'negative' group, then calculate ratio of expression difference to 'negative'. Sort by that, picking the first one from each cluster but skip over un-named genes (gmxxxx, etc)
markers.sox6 <- FindAllMarkers(lrrkdata.sox6, min.pct = .2, logfc.threshold = 1, only.pos = T, min.diff.pct = .3, recorrect_umi = FALSE)
markers.sox6$diff <- (markers.sox6$pct.1) - (markers.sox6$pct.2)
markers.sox6 <- filter(markers.sox6, pct.2 < 0.3)
markers.sox6$diffratio <- (markers.sox6$diff) / (markers.sox6$pct.2)
markers.calb1 <- FindAllMarkers(lrrkdata.calb1, min.pct = .2, logfc.threshold = 1, only.pos = T, min.diff.pct = .3, recorrect_umi = FALSE)
markers.calb1$diff <- (markers.calb1$pct.1) - (markers.calb1$pct.2)
markers.calb1 <- filter(markers.calb1, pct.2 < 0.3)
markers.calb1$diffratio <- (markers.calb1$diff) / (markers.calb1$pct.2)
markers.gad2 <- FindAllMarkers(lrrkdata.gad2, min.pct = .2, logfc.threshold = 1, only.pos = T, min.diff.pct = .3, recorrect_umi = FALSE)
markers.gad2$diff <- (markers.gad2$pct.1) - (markers.gad2$pct.2)
markers.gad2 <- filter(markers.gad2, pct.2 < 0.3)
markers.gad2$diffratio <- (markers.gad2$diff) / (markers.gad2$pct.2)

#Supp fig cluster stability calculations
Idents(lrrkdata.v5) <- "seurat_clusters"
#set downsample size
sample.size <- as.integer(0.8 * ncol(lrrkdata.v5))
set.seed(42)
#helper function for a single run of downsampling and calculating max jaccard for each old cluster
ClusterStability.singlerun <- function(object){
  object.downsample <- subset(object, cells = Clusterstability.cellnames(object), seed = NULL)
  message(head(cellnames))
  object.downsample <- RunPCA(object.downsample, verbose = F, seed.use = NULL)
  object.downsample <- FindNeighbors(object.downsample, dims = 1:32, reduction = "cca", n.trees = 500, k.param = 40, verbose = F)
  object.downsample <- FindClusters(object.downsample, resolution = 0.8, algorithm = 1, group.singletons = TRUE, graph.name = "SCT_snn", verbose = F)
  nclusters.old <- length(levels(object))
  nclusters.new <- length(levels(object.downsample))
  jac.max <- NULL
  jac.max <- vector(mode = "numeric", length = nclusters.old)
  for (i in 1:nclusters.old) {
    for (z in 1:nclusters.new) {
      cellNames.old <- WhichCells(object, idents = (i-1), seed = NULL)
      cellNames.New <- WhichCells(object.downsample, idents = (z-1), seed = NULL)
      jac.max[i] <- max(jac.max[i], calculateJac(cellNames.old, cellNames.New))
    }
  }
  jac.max
}
#help function for calculating jaccard similarity index
calculateJac <- function(set1, set2){
  in_length = length(intersect(set1, set2))
  un_length = length(union(set1, set2))
  jaccard = in_length/un_length
  return(jaccard)
}
#helper function for getting cell names that are used in downsampled dataset
Clusterstability.cellnames <- function(object){
  allCells <- Cells(object)
  cellnames <<- sample(allCells, size = sample.size, replace = F)
  cellnames
}
#wrapper function for calculating jaccard n times and saving it as a dataframe
Clusterstability.replicated <- function(object, reps){
  message("Running replication number 1")
  outputs <- data.frame(ClusterStability.singlerun(object))
  rownames(outputs) <- levels(object)
  for (i in 1:(reps - 1)) {
    message("Running replication number ", i+1)
    outputs <- cbind(outputs, ClusterStability.singlerun(object))
  }
  colnames(outputs) <- 1:reps
  outputs
}
#execute code for our dataset
#assumes a temporary copy of the seurat object exists for manipulation
Idents(templrrkdata) <- "seurat_clusters"
stability <- Clusterstability.replicated(templrrkdata, 100)
stability.norm <- (stability / 0.8)
medians <- rowMedians(as.matrix(stability))
medians.norm <- (medians / 0.8)

#Plot scDRS scores for cells
# Step 1: Extract UMAP coordinates and feature values
umap_coords <- Embeddings(scdrs, "umap.cca")
feature_values <- FetchData(scdrs, vars = "parkinsons_score")
# Combine into a data frame
plot_data <- data.frame(umap_coords, feature = feature_values$parkinsons_score)
# Step 2: Calculate opacity based on the absolute value of the feature, normalized to a 0-1 scale
plot_data$opacity <- abs(plot_data$feature) / max(abs(plot_data$feature))
# Step 3: Plot the data using ggplot2 with dynamic opacity
ggplot(plot_data, aes(x = umapcca_1, y = umapcca_2, color = feature)) +
  geom_point(aes(alpha = opacity), size = .9) +  # Use calculated opacity
  scale_color_gradient2(low = "red3", mid = "grey95", high = "blue3", midpoint = 0, limits = c(-4, 4)) +
  scale_alpha(range = c(0.01, .99), guide = 'none') +  # Adjust the range of opacity and remove the legend for alpha
  theme_void() +
  labs(title = "Parkinson's Disease GWAS Risk Score", color = "Score") +
  coord_fixed()

#Bootstrapping CIs for scDRS scores for each cluster
#the below pipeline can be further modified to address whatever identities/groups necessary, as was done in the paper to explore scores by cluster family etc
library(Seurat)
library(boot)
library(ggplot2)
# Extract the data from Seurat object
scores <- scdrs@meta.data$parkinsons_score
clusters <- scdrs@meta.data$seurat_clusters
# Create a data frame
data <- data.frame(score = scores, cluster = as.factor(clusters))
# Bootstrapping function to calculate the mean
boot_mean <- function(data, indices) {
  d <- data[indices, ]  # Sample with replacement
  return(mean(d$score))
}
# Running bootstrap for each cluster and storing results
set.seed(42)
results <- lapply(split(data, data$cluster), function(cluster_data) {
  boot(cluster_data, boot_mean, R=10000)  # 10000 bootstrap replicates
})
# Extracting the bootstrap results and calculating the 95% CI
ci_data <- do.call(rbind, lapply(names(results), function(name) {
  ci <- boot.ci(results[[name]], type="perc")
  data.frame(
    Cluster = name,
    Mean = results[[name]]$t0,
    Lower = ci$percent[4],
    Upper = ci$percent[5]
  )
}))
