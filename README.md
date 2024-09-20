# Gaertner_Oram_etal
Code for the Gaenter_Oram et al. paper 
Includes code for all main and supplemental figures

In this project, we perform single-nucleus RNAseq on midbrain dopaminergic neurons in order to classify molecular subtypes of these cells. We then perform several analyses to explore the differential gene expression patterns among clusters, as well as comparing between wildtype and Lrrk2 G2019S mice. 

To get started with our code, follow the stepwise instructions (Gaertner Oram RNA seq custom code) to integrate samples (available from GEO) into a Seurat dataset. From there, custom code is procided and annotated for all figures which require custom analyses in our source paper. Figures that utilize base Seurat functions (for example, performing differential expression on marker genes using the Seurat FindMarkers() command) are not included in this file, unless specific modifications have been made to the default parameters of these commands or if the output is required for subsequent custom analyses.


Supplemental fig code for Scanpy / Sankey supplemental figures can be found in Supp Figs folder.
Begin with the .ipynb Jupyter file and follow the steps chronologically. 
Next, use the seurat object from Gaertner Oram RNA seq custom code and the dataframe from the .ipynb file to create the sankey diagram. 

