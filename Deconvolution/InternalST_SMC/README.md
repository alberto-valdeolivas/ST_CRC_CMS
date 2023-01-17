# Deconvolution on our set of CRC samples using the SMC cohort

This folder contains the deconvolution results on  our set of CRC ST samples using as reference the annotations from the scRNA-seq generated from a Korean cohort (SMC dataset). It is organized as follows: 

* [07a_CEll2Loc_ExpressionSignature.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/InternalST_SMC/07a_CEll2Loc_ExpressionSignature.ipynb): extraccion of cell type specific signatures according to the annotations from the reference. 

* [07b_Cell2Loc_Deconvolution.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/InternalST_SMC/07b_Cell2Loc_Deconvolution.ipynb): the  deconvolution itself. Cell2Location assigns cell types to transcripts detected in every spot according to the signatures extracted in the previous script. 

* [07c_Cell2Loc_Deconvolution_Downstream_Analysis.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/InternalST_SMC/07c_Cell2Loc_Deconvolution_Downstream_Analysis.ipynb): downstream analysis on the deconvolution resutls containing a clustering based on cell abundance and a cell type co-localization analysis, among others. Among others, this script contains the co-localization analysis of the samples. The R script to plot these results in a format according with the remaining figures of the publication can be found [here](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/Deconvolution/InternalST_SMC/Pub_Reading_Factors.R). 

In addition, the notebooks containing the detailed results per sample for all the cell types are availabe at zenodo (Too heavy files for github). 
