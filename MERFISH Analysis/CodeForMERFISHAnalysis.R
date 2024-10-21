# Code by Cameron Oram for Gaertner & Oram et al 2024
# For questions please email: cameron.oram@mail.mcgill.ca

#Load in all the relevant libraries

library(Seurat)
library(SeuratData)
library(SeuratObject)
library(patchwork)
library(dplyr)
library(ggplot2)
library(rapport)
library(rapportools)
library(tidyverse)
library(progressr)
library(ape)
library(sctransform)
library(scales)
library(SeuratDisk)
#Necessary to make the RObjects V5
Seurat.object.assay.version = "v5"
set.seed(100)

## load libs ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(BiocParallel)    
})

mol.type <- "microns"
z.stack <- 3L # z stack to use

# Change the directories to match the location of the cell_by_gene.csv and cell_metadata.csv files
folder.A<-"E:/Cameron/MERSCOPEData/OutputFiles/202303131238_CO-D852SampVer-13Mar2023-500genePanel_VMSC07702_cellpose/"
folder.B<-"E:/Cameron/MERSCOPEData/OutputFiles/202303201243_D852CMidbrainPeptide20Mar23_VMSC07702_cellpose/"
folder.C<-"E:/Cameron/MERSCOPEData/OutputFiles/202303221620_D852-DMidbrainPeptidesRun2_VMSC07702_cellpose/"
folder.D<-"E:/Cameron/MERSCOPEData/OutputFiles/202304041053_CO-WTA-A-Midbrain-04Apr2023_VMSC07702_cellpose/"
folder.E<-"E:/Cameron/MERSCOPEData/3rdOutput/OutputFilesWTD-C_29Jul2023/"
folder.F<-"E:/Cameron/MERSCOPEData/2ndBatchOutput/WTD-D/JustFiles/"
folder.G<-"E:/Cameron/MERSCOPEData/3rdOutput/OutputFilesWTD-E_31Jul2023/"

# Load the csvs containing the different combination of gene names to use:
# Commissural Genes are the genes that show similar expression levels between sn-Seq and spatial.
subGeneNames<-read.csv(file="E:/Cameron/MERSCOPEData/GeneNames/CommissuralGenes_21Aug2023.csv")
subGeneNames<-subGeneNames$Genes
#No glial genes are the full panel of 500 genes minus the glial cytoskeleton/myelin associated genes
noGlialNames<-read.csv(file="E:/Cameron/MERSCOPEData/GeneNames/geneNamesNoGlia.csv")
noGlialNames<-noGlialNames$Genes

#RemovedGenes are the gene names selected by calculating the difference between transcripts detected in snSeq and spatSeq datasets
removedGeneNames<-read.csv(file="E:/Cameron/MERSCOPEData/GeneNames/GeneNamesRemove.csv")
removedGeneNames<-removedGeneNames$Gene.Name

# load Each slice individually
SliceA <-
  LoadVizgen_opt(data.dir = folder.A,
                 fov = "SliceA", 
                 assay = "RNA",
                 metadata = c("volume", "fov"), # add cell volume info
                 mol.type = mol.type, # molecule coords in µm spa
                 type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                 z = z.stack,
                 add.zIndex = TRUE, # add z slice section to the object
                 update.object = TRUE,
                 use.BiocParallel = TRUE, # using `BiocParallel` instead of `future`
                 DTthreads.pct = NULL, # percentage of total threads to use for `data.table::fread`
                 verbose=TRUE
  )
SliceB <-
  LoadVizgen_opt(data.dir = folder.B,
                 fov = "SliceB", 
                 assay = "RNA",
                 metadata = c("volume", "fov"), # add cell volume info
                 mol.type = mol.type, # molecule coords in µm spa
                 type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                 z = z.stack,
                 add.zIndex = TRUE, # add z slice section to the object
                 update.object = TRUE,
                 use.BiocParallel = TRUE, # using `BiocParallel` instead of `future`
                 DTthreads.pct = NULL, # percentage of total threads to use for `data.table::fread`
                 verbose=TRUE
  )
SliceC <-
  LoadVizgen_opt(data.dir = folder.C,
                 fov = "SliceC", 
                 assay = "RNA",
                 metadata = c("volume", "fov"), # add cell volume info
                 mol.type = mol.type, # molecule coords in µm spa
                 type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                 z = z.stack,
                 add.zIndex = TRUE, # add z slice section to the object
                 update.object = TRUE,
                 use.BiocParallel = TRUE, # using `BiocParallel` instead of `future`
                 DTthreads.pct = NULL, # percentage of total threads to use for `data.table::fread`
                 verbose=TRUE
  )
SliceD <-
  LoadVizgen_opt(data.dir = folder.D,
                 fov = "SliceD", 
                 assay = "RNA",
                 metadata = c("volume", "fov"), # add cell volume info
                 mol.type = mol.type, # molecule coords in µm spa
                 type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                 z = z.stack,
                 add.zIndex = TRUE, # add z slice section to the object
                 update.object = TRUE,
                 use.BiocParallel = TRUE, # using `BiocParallel` instead of `future`
                 DTthreads.pct = NULL, # percentage of total threads to use for `data.table::fread`
                 verbose=TRUE
  )
SliceE <-
  LoadVizgen_opt(data.dir = folder.E,
                 fov = "SliceE", 
                 assay = "RNA",
                 metadata = c("volume", "fov"), # add cell volume info
                 mol.type = mol.type, # molecule coords in µm spa
                 type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                 z = z.stack,
                 add.zIndex = TRUE, # add z slice section to the object
                 update.object = TRUE,
                 use.BiocParallel = TRUE, # using `BiocParallel` instead of `future`
                 DTthreads.pct = NULL, # percentage of total threads to use for `data.table::fread`
                 verbose=TRUE
  )
SliceF <-
  LoadVizgen_opt(data.dir = folder.F,
                 fov = "SliceF", 
                 assay = "RNA",
                 metadata = c("volume", "fov"), # add cell volume info
                 mol.type = mol.type, # molecule coords in µm spa
                 type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                 z = z.stack,
                 add.zIndex = TRUE, # add z slice section to the object
                 update.object = TRUE,
                 use.BiocParallel = TRUE, # using `BiocParallel` instead of `future`
                 DTthreads.pct = NULL, # percentage of total threads to use for `data.table::fread`
                 verbose=TRUE
  )

SliceG <-
  LoadVizgen_opt(data.dir = folder.G,
                 fov = "SliceG", 
                 assay = "RNA",
                 metadata = c("volume", "fov"), # add cell volume info
                 mol.type = mol.type, # molecule coords in µm spa
                 type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                 z = z.stack,
                 add.zIndex = TRUE, # add z slice section to the object
                 update.object = TRUE,
                 use.BiocParallel = TRUE, # using `BiocParallel` instead of `future`
                 DTthreads.pct = NULL, # percentage of total threads to use for `data.table::fread`
                 verbose=TRUE
  )

DefaultAssay(SliceA) <- "RNA"
DefaultAssay(SliceB) <- "RNA"
DefaultAssay(SliceC) <- "RNA"
DefaultAssay(SliceD) <- "RNA"
DefaultAssay(SliceE) <- "RNA"
DefaultAssay(SliceF) <- "RNA"
DefaultAssay(SliceG) <- "RNA"
SliceA[['sample']] <- 'sliceA'
SliceB[['sample']] <- 'sliceB'
SliceC[['sample']] <- 'sliceC'
SliceD[['sample']] <- 'sliceD'
SliceE[['sample']] <- 'sliceE'
SliceF[['sample']] <- 'sliceF'
SliceG[['sample']] <- 'sliceG'

SliceA[['RNA']] <- as(object = SliceA[['RNA']], Class = "Assay5")
SliceB[['RNA']] <- as(object = SliceB[['RNA']], Class = "Assay5")
SliceC[['RNA']] <- as(object = SliceC[['RNA']], Class = "Assay5")
SliceD[['RNA']] <- as(object = SliceD[['RNA']], Class = "Assay5")
SliceE[['RNA']] <- as(object = SliceE[['RNA']], Class = "Assay5")
SliceF[['RNA']] <- as(object = SliceF[['RNA']], Class = "Assay5")
SliceG[['RNA']] <- as(object = SliceG[['RNA']], Class = "Assay5")

# Filter out cells based on 1) # of transcripts per cell, 2) Size of Cell
# Final is Volume <4000. This can be changed as necessary
SliceA<-subset(SliceA,subset=nFeature_RNA>40 & volume>500 & volume<4000)
SliceB<-subset(SliceB,subset=nFeature_RNA>40 & volume>500 & volume<4000)
SliceC<-subset(SliceC,subset=nFeature_RNA>40 & volume>500 & volume<4000)
SliceD<-subset(SliceD,subset=nFeature_RNA>40 & volume>500 & volume<4000)
SliceE<-subset(SliceE,subset=nFeature_RNA>40 & volume>500 & volume<4000)
SliceF<-subset(SliceF,subset=nFeature_RNA>40 & volume>500 & volume<4000)
SliceG<-subset(SliceG,subset=nFeature_RNA>40 & volume>500 & volume<4000)

merscope.v5 <-merge(SliceA, list(SliceB, SliceC, SliceD, SliceE, SliceF, SliceG))

# Now remove the old variables to save space
rm(SliceA)
rm(SliceB)
rm(SliceC)
rm(SliceD)
rm(SliceE)
rm(SliceF)
rm(SliceG)

# Run normalization and integration using RPCAIntegration. Other options include Harmony and CCAIntegration
merscope.v5 <- SCTransform(merscope.v5, vst.flavor = "v2")
merscope.v5 <- RunPCA(merscope.v5)

merscope.v5 <- IntegrateLayers(object = merscope.v5, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", assay = "SCT", normalization.method = "SCT")
merscope.v5 <- RunUMAP(merscope.v5, reduction = "integrated.rpca", dims = 1:37, reduction.name = "umap.rpca")
merscope.v5 <- FindNeighbors(merscope.v5, dims = 1:37, reduction = "integrated.rpca")
merscope.v5 <- FindClusters(merscope.v5, resolution = 0.3, cluster.name = "rpca.clusters")

# Get a list of markers for each cluster to use for manual annotation
merscope.v5<-PrepSCTFindMarkers(merscope.v5,assay="SCT",verbose=TRUE)
markers.all.WholeBrain<-FindAllMarkers(merscope.v5,test="roc",assay="SCT")

#At this point I compiled an excel spreadsheet of top gene markers for each cluster

# Load simplified Annotations for the whole brain:
Idents(merscope.v5)<-"seurat_clusters"
ClusterAnnotations<-read.csv(file="E:/Cameron/MERSCOPEData/GeneMarkers/4000Vol/WholeBrainAnnotations.csv")

#Full annotation breaks down all 40 clusters
FullAnnotations<-ClusterAnnotations$Specific.Annotations
LowAnnotations<-ClusterAnnotations$Low.Annotation
SpecificAnnotations<-ClusterAnnotations$Specific.Annotations

#Simple annotations breaks them down into more generalized clusters
SimpleAnnotations<-ClusterAnnotations$Mid.Annotation
HighLevelAnnotations<-ClusterAnnotations$High.Level.Annotation
merscope.annotated<-merscope.v5
merscope.Specific<-merscope.v5

merscope.SimpAnnotated<-merscope.v5
merscope.HighAnnotation<-merscope.v5
names(FullAnnotations)<-levels(merscope.annotated)
merscope.annotated<-RenameIdents(merscope.annotated,FullAnnotations)

names(LowAnnotations)<-levels(merscope.v5)
merscope.v5<-RenameIdents(merscope.v5,LowAnnotations)
merscope.v5<-AddMetaData(merscope.v5,merscope.v5@active.ident,col.name="LowAnnot")

names(SimpleAnnotations)<-levels(merscope.v5)
merscope.v5<-RenameIdents(merscope.SimpAnnotated,SimpleAnnotations)
merscope.v5<-AddMetaData(merscope.v5,merscope.v5@active.ident,col.name="SimpAnnote")

names(HighLevelAnnotations)<-levels(merscope.HighAnnotation)
merscope.v5<-RenameIdents(merscope.v5,HighLevelAnnotations)
merscope.v5<-AddMetaData(merscope.v5,merscope.v5@active.ident,col.name="HighAnnote")

names(SpecificAnnotations)<-levels(merscope.v5)
merscope.v5<-RenameIdents(merscope.v5,SpecificAnnotations)
merscope.v5<-AddMetaData(merscope.v5,merscope.v5@active.ident,col.name="SpecificAnnote")

# Subset out just the Dopamine neurons from cluster 25:
merscope.DAsubset.v5 <- subset(merscope.v5, rpca.clusters == "25",features=noGlialNames)

# Now run normalization for all of these dopamine subsets
DefaultAssay(merscope.DAsubset.v5)<-"RNA"
merscope.DAsubset.v5<-SCTransform(merscope.DAsubset.v5,vst.flavor="v2",return.only.var.genes=FALSE)

#Based on the ElbowPlot ran the asymptote appears at Dim 19
merscope.DAsubset.v5<-RunPCA(merscope.DAsubset.v5,reduction.name="SCT.post")
merscope.DAsubset.v5<-RunUMAP(merscope.DAsubset.v5,reduction="SCT.post",dims=1:19,reduction.name="SCT.post.umap")
merscope.DAsubset.v5<-FindNeighbors(merscope.DAsubset.v5,dims=1:19,reduction="SCT.post",graph.name="SCT.post.graph")
merscope.DAsubset.v5<-FindClusters(merscope.DAsubset.v5,resolution=0.6,cluster.name="SCT.subclusters",graph.name="SCT.post.graph")

#In this dataset Cluster 10 contains 0 levels of Th and Slc6a3 and clusters very far from the rest of the cells. We removed these to isolate only DA
merscope.DAsubset.v5<-subset(merscope.DAsubset.v5,SCT.subclusters=="10",invert=TRUE)
DefaultAssay(merscope.DAsubset.v5)<-"RNA"
merscope.DAsubset.v5<-SCTransform(merscope.DAsubset.v5,vst.flavor="v2",return.only.var.genes=FALSE)
merscope.DAsubset.v5<-RunPCA(merscope.DAsubset.v5,reduction.name="SCT.post")
merscope.DAsubset.v5<-RunUMAP(merscope.DAsubset.v5,reduction="SCT.post",dims=1:19,reduction.name="SCT.post.umap")
merscope.DAsubset.v5<-FindNeighbors(merscope.DAsubset.v5,dims=1:19,reduction="SCT.post",graph.name="SCT.post.graph",k.param=20)
merscope.DAsubset.v5<-FindClusters(merscope.DAsubset.v5,resolution=0.55,cluster.name="SCT.subclusters.0.8",graph.name="SCT.post.graph")

merscope.DAsubset.v5<-PrepSCTFindMarkers(merscope.DAsubset.v5,assay="SCT",verbose=TRUE)
markers.DAsubset<-FindAllMarkers(merscope.DAsubset.v5,test="roc",assay="SCT")


# Add Spatial Localization information to each cell:

# Identify the markers associated with each DA cluster and add them to the subset dataset
merscope.DAsubset.v5<-PrepSCTFindMarkers(merscope.DAsubset.v5,assay="SCT",verbose=TRUE)
markers.DAsubset<-FindAllMarkers(merscope.DAsubset.v5,test="roc",assay="SCT")
# Add manual annotations to the DASubset 
ClusterAnnotations<-read.csv(file="E:/Cameron/MERSCOPEData/GeneMarkers/4000Vol/DASubtypeAnnotations.csv")
ClusterAnnotations<-ClusterAnnotations$Annotations.
Clustermarkers<-read.csv(file="E:/Cameron/MERSCOPEData/GeneMarkers/4000Vol/GenesForHeatmap.csv")
Clustermarkers<-Clustermarkers$Genes
merscope.DAsubset.v5.annotated<-merscope.DAsubset.v5
names(ClusterAnnotations)<-levels(merscope.DAsubset.v5.annotated)
merscope.DAsubset.v5.annotated<-RenameIdents(merscope.DAsubset.v5.annotated,ClusterAnnotations)

# Now load and integrate the sn RNA-seq with the MERFISh dataset:
ZachData<-readRDS(file="E:/Cameron/MERSCOPEData/ZachsData/Lrrk2FinalDataset")

DefaultAssay(ZachData)<-"RNA"
ZachData<-SCTransform(ZachData,vst.flavor = "v2")
ZachData <- RunPCA(ZachData)
ZachData <- RunUMAP(ZachData, dims = 1:30, reduction.name = "umap.rpca",return.model=TRUE)
ZachData <- FindNeighbors(ZachData, dims = 1:30,return.neighbor=TRUE)
ZachData <- FindClusters(ZachData, resolution = .8, cluster.name = "rpca.clusters")

subZach<-subset(x=ZachData,features=removedGeneNames)
DefaultAssay(subZach)<-"RNA"
subZach<-SCTransform(subZach,vst.flavor = "v2")
subZach <- RunPCA(subZach)
subZach <- RunUMAP(subZach, dims = 1:30, reduction.name = "umap.rpca",return.model=TRUE)
subZach <- FindNeighbors(subZach, dims = 1:30,return.neighbor=TRUE)
subZach <- FindClusters(subZach, resolution = 0.6, cluster.name = "rpca_clusters")

# Prep the DA subset data for integration:
merscope.DAsubset.v5<-JoinLayers(merscope.DAsubset.v5)
merscope.DAsubset.v5<-AddMetaData(merscope.DAsubset.v5,"spatial",col.name="group")

reference <- ZachData
reference <- RunUMAP(reference, reduction = "cca", dims = 1:32, reduction.name = "umap.cca", n.epochs = 500, min.dist = .2, n.neighbors = 1000,return.model=TRUE)
subreference<-subset(reference,features=subGeneNames)

reference<-JoinLayers(reference)
object3<-merscope.DAsubset.v5
annotations<-ZachData@meta.data$seurat_clusters

anchor2 <- FindTransferAnchors(
  reference = subreference,
  query = object3,
  reduction="cca",
  normalization.method = "SCT",
  recompute.residuals = FALSE
)

object3 <- MapQuery(
  anchorset = anchor2,
  query = object3,
  refdata=annotations,
  reference = subreference,
  reference.reduction="pca",
  reduction.model = "umap.cca",
  reference.dims=1:50
)

# Impute the single-cell data from the snSEQ to merfish:

object3<-TransferData(
  anchorset=anchor2,
  query=object3,
  refdata="RNA",
  reference=reference,
  weight.reduction="cca"
)

# Only keep cells with high transfer prediction score: 
confidentCells<-subset(object3,predicted.id.score>0.5)

## Figure Generation:

# Figure 3:
Dimplot(merscope.v5,group.by="SimpAnnote",reduction="rcpa.umap")
FeaturePlot(merscope.v5,features="Th")
FeaturePlot(merscope.v5,features="Slc6a3")
FeaturePlot(merscope.v5,features="Slc18a2")

ImageDimPlot(merscope.v5,fov="SliceA",group.by="SimpAnnote",border.size=NA,size=1.5)
ImageDimPlot(merscope.v5,fov="SliceE",group.by="SimpAnnote",border.size=NA,size=1.5)
ImageDimPlot(merscope.v5,fov="SliceB",group.by="SimpAnnote",border.size=NA,size=1.5)

ImageFeaturePlot(merscope.v5,fov="SliceA",features="Th")
ImageFeaturePlot(merscope.v5,fov="SliceE",features="Th")
ImageFeaturePlot(merscope.v5,fov="SliceB",features="Th")

# Figure 4:
# Figure 4b:
ImageDimPlot(confidentCells, group.by="predicted.id", fov="SliceA", border.size=NA, size=1.5)
ImageDimPlot(confidentCells, group.by="predicted.id", fov="SliceB", border.size=NA, size=1.5)
ImageDimPlot(confidentCells, group.by="predicted.id", fov="SliceE", border.size=NA, size=1.5)

# Figure 4C-E:
DefaultAssay(confidentCells)<-"predicted.id"
ImageFeaturePlot(confidentCells,features="Sox6",fov="SliceA")
ImageFeaturePlot(confidentCells,features="Sox6",fov="SliceB")
ImageFeaturePlot(confidentCells,features="Sox6",fov="SliceE")

ImageFeaturePlot(confidentCells,features="Calb1",fov="SliceA")
ImageFeaturePlot(confidentCells,features="Calb1",fov="SliceB")
ImageFeaturePlot(confidentCells,features="Calb1",fov="SliceE")

ImageFeaturePlot(confidentCells,features="Gad2",fov="SliceA")
ImageFeaturePlot(confidentCells,features="Gad2",fov="SliceB")
ImageFeaturePlot(confidentCells,features="Gad2",fov="SliceE")

# Figure 5:
# FeaturePlots of Single-nuclei dataset and spatial Feature Plots for cluster identifiers:
FeaturePlot(ZachData,features="Tacr3")
FeaturePlot(ZachData,features="March3")
FeaturePlot(ZachData,features="Hs6st3")
FeaturePlot(ZachData,features="Asic2")
FeaturePlot(ZachData,features="Kcnmb2")
FeaturePlot(ZachData,features="Tenm3")

ImageFeaturePlot(confidentCells,fov="SliceA",features="Tacr3",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="March3",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Hs6st3",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Asic2",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Kcnmb2",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Tenm3",border.size=NA,size=1.5)

# Figure 6:
FeaturePlot(ZachData,features="Ptprt")
FeaturePlot(ZachData,features="Sema5b")
FeaturePlot(ZachData,features="Mob3b")
FeaturePlot(ZachData,features="Col23a1")
FeaturePlot(ZachData,features="Kctd8")
FeaturePlot(ZachData,features="Pcdh7")
FeaturePlot(ZachData,features="Pde11a")
FeaturePlot(ZachData,features="Postn")

ImageFeaturePlot(confidentCells,fov="SliceA",features="Ptprt",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Sema5b",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Mob3b",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Col23a1",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Kctd8",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Pcdh7",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Pde11a",border.size=NA,size=1.5)
ImageFeaturePlot(confidentCells,fov="SliceA",features="Postn",border.size=NA,size=1.5)


# Supp Fig 4:

VlnPlot(merfish.v5,features="nCount_RNA")
VlnPlot(merfish.v5,features="nFeature_RNA")
VlnPlot(merfish.v5,features="volume")


# Supp Fig 5:
DimPlot(merfish.DASubset.v5,reduction="SCT.post")
FeaturePlot(merfish.DASubset.v5,features="Th")
FeaturePlot(merfish.DASubset.v5,features="Slc6a3")
FeaturePlot(merfish.DASubset.v5,features="Slc18a2")
DimPlot(merfish.DASubset.v5,group.by="MidbrainLocalization")

# Supp Fig 6:
DimPlot(ZachData,reduction="umap.rpca")
DimPlot(confidentCells,reduction="ref.umap",group.by="predicted.id")
FeaturePlot(confidentCells,features="predicted.id.score",reduction="ref.umap")

FeaturePlot(ZachData,reduction="umap.rpca",features="Gad2")
FeaturePlot(ZachData,reduction="umap.rpca",features="Calb1")
FeaturePlot(ZachData,reduction="umap.rpca",features="Aldh1a1")

DefaultAssay(confidentCells)<-"RNA"
FeaturePlot(confidentCells,reduction="ref.umap",features="Gad2")
FeaturePlot(confidentCells,reduction="ref.umap",features="Calb1")
FeaturePlot(confidentCells,reduction="ref.umap",features="Aldh1a1")
ImageFeaturePlot(confidentCells,fov="SliceE",features="Gad2")
ImageFeaturePlot(confidentCells,fov="SliceE",features="Calb1")
ImageFeaturePlot(confidentCells,fov="SliceE",features="Aldh1a1")

DefaultAssay(confidentCells)<-"predicted.id"
FeaturePlot(confidentCells,reduction="ref.umap",features="Gad2")
FeaturePlot(confidentCells,reduction="ref.umap",features="Calb1")
FeaturePlot(confidentCells,reduction="ref.umap",features="Aldh1a1")
ImageFeaturePlot(confidentCells,fov="SliceE",features="Gad2")
ImageFeaturePlot(confidentCells,fov="SliceE",features="Calb1")
ImageFeaturePlot(confidentCells,fov="SliceE",features="Aldh1a1")
