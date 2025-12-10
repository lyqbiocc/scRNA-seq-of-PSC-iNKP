Title: Dual preventions of tumor occurrence and relapse via pluripotent stem cell-derived NK progenitor cell therapy  
https://www.biorxiv.org/content/10.1101/2025.01.07.631650v1  
Summary:  
Raw data of scRNA-seq can be downloaded and processed by CellRanger software package (version 7.0.1, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation), then aggregated by the ‘aggr’ function of the CellRanger package and subjected to Seurat (version 4.2.0) for further analysis.  
Raw data of scRNA-seq:  
The scRNA-seq data of UCB_resting NK cells, PSC-derived cells from day 16-20 organoids, and CD19CAR-R4iNK cells from PB, BM, Liver, Lung, and Spleen have been deposited in the Genome Sequence Archive public database (HRA008662). The scRNA-seq data of natural NKP and NK are from the Gene Expression Omnibus database (GSE149938). The accession number scRNA-seq data of UCB_activated NK is HRA007978. The accession number scRNA-seq data of PSC-iNK, PB_resting NK cells, and PB_activated cells is HRA001609.  

# scRNA seq analysis of hPSC derived iNKP and iNK cells

## Analysis of all the GFP negative cells of the hPSC-derived organoids from day 16, 17, 18, 19 and 20
    Figure S1A-B
    01_Identify_the_cellpopulation_D16_D20.R
    plot_the_cellpopulation_D16_D20.R

## Locate the iNKP-like cells and identify iNK lineage cells 
    Figure 1B-C
    02_Identify_the_iNKP_D16_D20.R
    03_pseudotime_the_iNKP_D16_D20.R
    plot_the_iNKP_D16_D20.R
    plot_pseudotime_the_iNKP_D16_D20.R

## Analysis of hPSC-derived iNK isolated from different tissue of iNKP recipients
    Figure 2L-N
    04_analysis_iNKP_derived_iNK_invivo.R
    05_correlation_analysis.R
