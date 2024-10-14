library(Seurat)
library(ShinyCell)

#Read in your Seurat object
lrrk2_object <- readRDS('Lrrk2FinalDataset')

nonDA<-WhichCells(lrrk2_object, idents = c('16','21'))
lrrk2_object <- subset(lrrk2_object, cells = nonDA, invert = T)

scConf = createConfig(lrrk2_object)
scConf = modMetaName(scConf, 
                     meta.to.mod = c("nCount_SCT", "nFeature_SCT", "percent.mt", "seurat_clusters"), 
                     new.name = c("No. UMIs", "No. detected genes",
                                  "% MT genes", "clusters"))
scConf = delMeta(scConf, 
                 meta.to.del = 'orig.ident')
makeShinyApp(lrrk2_object, scConf, gene.mapping = TRUE, gex.assay = 'SCT',
             shiny.title = "Dopabase") 

library(rsconnect)
rsconnect::setAccountInfo(name='insertnamehere',
                          token='inserttokenhere',
                          secret='insertsecrethere')
rsconnect::deployApp('shinyappname')