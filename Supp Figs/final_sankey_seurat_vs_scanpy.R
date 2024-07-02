library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsankey)
library(devtools)
library(stringr)
library(networkD3)

lrrk2_object<-readRDS('Lrrk2FinalDataset')

cluster1<-subset(lrrk2_object, idents=1)
cluster1_R<-data.frame(cluster1$seurat_clusters)

cluster_R<-data.frame(lrrk2_object$seurat_clusters)


cluster_R_edited<-sub("1_1", "control1", rownames(cluster_R))
cluster_R_edited<-sub("1_2", "control2", cluster_R_edited)
cluster_R_edited<-sub("1_3", "lrrk2_1", cluster_R_edited)
cluster_R_edited<-sub("1_4", "lrrk2_2", cluster_R_edited)

df_seurat<-data.frame(cluster_R_edited,cluster_R) #df of cluster names and renamed cellnames
df_seurat<-df_seurat %>% rename(cell_name=cluster_R_edited)

#loading in scanpy cellnames w according sample 
samplesleiden<- data.frame(read.csv('samples_cellname_leiden'))
#1-0-0-0 control1 , 1-1-0-0 lrrk2_1, 1-1-0 lrrk2_2, 1-1 control2


#loading in scanpy cluster numbers
cluster_scanpy_leiden<- data.frame(read.csv('clusters_scanpy_leiden.csv'))
rownames(cluster_scanpy_leiden)<- cluster_scanpy_leiden[,1]



cluster_scanpy_leiden_edited<-sub("1-0-0-0", "control1", rownames(cluster_scanpy_leiden))
cluster_scanpy_leiden_edited<-sub("1-1-0-0", "control2", cluster_scanpy_leiden_edited)
cluster_scanpy_leiden_edited<-sub("1-1-0", "lrrk2_1", cluster_scanpy_leiden_edited)
cluster_scanpy_leiden_edited<-sub("1-1", "lrrk2_2", cluster_scanpy_leiden_edited)

df_leiden<-data.frame(cluster_scanpy_leiden_edited,cluster_scanpy_leiden) 
rownames(df_leiden)<-cluster_scanpy_leiden_edited
df_leiden<-df_leiden %>% rename(cell_name=cluster_scanpy_leiden_edited)


all_data = merge(df_leiden, df_seurat, by ="cell_name")
all_data<-all_data%>% select(-'X.x')

#create sankey df for leiden

sankey_df_leiden= data.frame(matrix(ncol = 3, nrow = 0))
columns_sankey_leiden <- c('seurat', 'leiden','cell_counts')
colnames(sankey_df_leiden)<-columns_sankey_leiden

for (x in 0:21) {
  for (y in 0:19) {
    output<-length(rownames(subset(all_data, lrrk2_object.seurat_clusters == x & leiden == y)))
    sankey_df_leiden<-rbind(sankey_df_leiden,list(x,(y+22),output))
  }
}

columns_sankey_leiden <- c('seurat', 'leiden','cell_counts')
colnames(sankey_df_leiden)<-columns_sankey_leiden

#remove cell counts that are less than the .05(total)
  

for (x in 0:21) {
  sumcells<-0
  for (y in 22:41) {
   sumcells_index<-which(sankey_df_leiden$seurat==x &sankey_df_leiden$leiden==y)
   sumcells<-sankey_df_leiden$cell_counts[sumcells_index]+sumcells
  }
  for (y in 22:41) {
    sumcells_index<-which(sankey_df_leiden$seurat==x &sankey_df_leiden$leiden==y)
    if (sankey_df_leiden$cell_counts[sumcells_index] > (.05*sumcells)){
      sankey_df_leiden$cell_counts[sumcells_index]<-sankey_df_leiden$cell_counts[sumcells_index]
      }else{
        sankey_df_leiden$cell_counts[sumcells_index]<-NA
      }
    
  }
}

sankey_df_leiden <-na.omit(sankey_df_leiden)
#create sankey leiden vs seurat


nodes = data.frame("name" = c("Seurat Cluster 0 ", "Seurat Cluster 1","Seurat Cluster 2", "Seurat Cluster 3", "Seurat Cluster 4", 
"Seurat Cluster 5","Seurat Cluster 6", "Seurat Cluster 7", "8 Seurat Cluster 8", "9 Seurat Cluster 9",
"Seurat Cluster 10", "Seurat Cluster 11", "Seurat Cluster 12", "13 Seurat Cluster 13","14 Seurat Cluster 14",
"Seurat Cluster 15", "Seurat Cluster 16", "Seurat Cluster 17","Seurat Cluster 18", "Seurat Cluster 19",
"Seurat Cluster 20", "Seurat Cluster 21","Scanpy Cluster 0", "Scanpy Cluster 1","Scanpy Cluster 2",
"Scanpy Cluster 3", "Scanpy Cluster 4", "Scanpy Cluster 5","Scanpy Cluster 6", "Scanpy Cluster 7", "Scanpy Cluster 8",
"Scanpy Cluster 9","Scanpy Cluster 10", "Scanpy Cluster 11", "Scanpy Cluster 12", "Scanpy Cluster 13","Scanpy Cluster 14",
"Scanpy Cluster 15", "Scanpy Cluster 16", "Scanpy Cluster 17","Scanpy Cluster 18", "Scanpy Cluster 19"))

my_color <- 'd3.scaleOrdinal() .domain(["Seurat Cluster 0 ", "Seurat Cluster 1","Seurat Cluster 2", "Seurat Cluster 3", "Seurat Cluster 4",
"Seurat Cluster 5","Seurat Cluster 6", "Seurat Cluster 7", "8 Seurat Cluster 8", "9 Seurat Cluster 9",
"Seurat Cluster 10", "Seurat Cluster 11", "Seurat Cluster 12", "13 Seurat Cluster 13","14 Seurat Cluster 14",
"Seurat Cluster 15", "Seurat Cluster 16", "Seurat Cluster 17","Seurat Cluster 18", "Seurat Cluster 19",
"Seurat Cluster 20", "Seurat Cluster 21","Scanpy Cluster 0", "Scanpy Cluster 1","Scanpy Cluster 2",
"Scanpy Cluster 3", "Scanpy Cluster 4", "Scanpy Cluster 5","Scanpy Cluster 6", "Scanpy Cluster 7", "Scanpy Cluster 8",
"Scanpy Cluster 9","Scanpy Cluster 10", "Scanpy Cluster 11", "Scanpy Cluster 12", "Scanpy Cluster 13","Scanpy Cluster 14",
"Scanpy Cluster 15", "Scanpy Cluster 16", "Scanpy Cluster 17","Scanpy Cluster 18", "Scanpy Cluster 19"]) #.range(["blue"])'


#need to rename for asthetics for the paper
nodes = data.frame("name" = c("seu  0 ", "seu  1","seu  2", "seu  3", "seu  4", 
                              "seu  5","seu  6", "seu  7", "8 seu  8", "9 seu  9",
                              "seu  10", "seu  11", "-----", "seu  13","seu  14", 
                              "seu  15", "seu  16", "seu  17","-", "--", 
                              "----", "seu  21","sc  0", "sc  1","sc  2", 
                              "sc  3", "sc  4", "sc  5","sc  6", "sc  7", "sc  8", 
                              "sc  9","sc  10", "sc  11", "sc  12", "sc  13","sc  14", 
                              "sc  15", "....", "...","..", "."))

my_color <- 'd3.scaleOrdinal() .domain(["seu  0 ", "seu  1","seu  2", "seu  3", "seu  4", 
                              "seu  5","seu  6", "seu  7", "8 seu  8", "9 seu  9",
                              "seu  10", "seu  11", "-----", "13 seu  13","14 seu  14", 
                              "seu  15", "seu  16", "seu  17","-", "--", 
                              "----", "seu  21","sc  0", "sc  1","sc  2", 
                              "sc  3", "sc  4", "sc  5","sc  6", "sc  7", "sc  8", 
                              "sc  9","sc  10", "sc  11", "sc  12", "sc  13","sc  14", 
                              "sc  15", "....", "...","..", "."]) .range(["black"])'
sankey_df_leiden$group = c("type_1", 
                "type_2",
                "type_3", 
                "type_4",
                "type_5",
                "type_6",
                "type_7",
                "type_8",
                "type_9", 
                "type_10",
                "type_11", 
                "type_12",
                "type_13",
                "type_14",
                "type_15",
                "type_16",
                "type_17", 
                "type_18",
                "type_19", 
                "type_20",
                "type_21",
                "type_22",
                "type_23",
                "type_24",
                "type_25", 
                "type_26",
                "type_27", 
                "type_28",
                "type_29",
                "type_30",
                "type_31",
                "type_32","type_33","type_34")

sankeyNetwork(Links = sankey_df_leiden, Nodes = nodes,
              Source = 'seurat', Target = 'leiden',
              Value = 'cell_counts', NodeID = "name", LinkGroup="group")


