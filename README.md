# R_Shiny_Visualization
Visualization of RNA sequence count data

This R script produces an R shiny web application that visualizes RNA-seq count data. 

This application was created for the visualization of RNA-seq count data from O'Meara et. al. However, the application can be utilized for any count data. 

Any original count data can be uploaded and the metadata is extracted and shown in the first tab. The data must be formatted as experimental parameters as columns and the 
genes as rows. 

In the second tab, various visualization methods for the count data are available. The side panel provides options to adjust for variation and non-zero samples in the 
visualization. The tab presents statistics, including the count of genes passing and failing the filter. Additionally, it generates two scatter plots depicting median count
vs. variance and median count vs. the number of non-zero samples. A heatmap and a principal component analysis plot are also featured, allowing users to determine the 
necessary principal components for visualization. 

Moving to the third tab, it conducts a differential expression analysis using DESeq2. The output includes a table with a search feature, enabling users to find any gene and 
access its corresponding differential expression data. A customizable volcano plot is also part of this tab, allowing users to adjust the X and Y axes and magnitude.

The final tab offers diverse options for individual gene visualization. Users can search for a specific gene, choose the plot type, and select a category based on the 
extracted metadata.

Citation of the original paper:
O'Meara CC, Wamstad JA, Gladstone RA, Fomovsky GM, Butty VL, Shrikumar A, Gannon JB, Boyer LA, Lee RT. Transcriptional reversion of cardiac myocyte fate during mammalian 
cardiac regeneration. Circ Res. 2015 Feb 27;116(5):804-15. doi: 10.1161/CIRCRESAHA.116.304269. Epub 2014 Dec 4. PMID: 25477501; PMCID: PMC4344930.
