# Deconvolution on our set of CRC samples using the KUL3 cohort

This folder contains the deconvolution results on  our set of CRC ST samples using as reference the annotations from the scRNA-seq generated from a Belgian cohort (KUL3 dataset). It is organized as follows: 

* [07a_CEll2Loc_ExpressionSignature-KUL3_cohort.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/InternalST_KUL3/07a_CEll2Loc_ExpressionSignature-KUL3_cohort.ipynb): extraccion of cell type specific signatures according to the annotations from the reference. 

* [07b_Cell2Loc_Deconvolution-KUL3_cohort.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/InternalST_KUL3/07b_Cell2Loc_Deconvolution-KUL3_cohort.ipynb): the  deconvolution itself. Cell2Location assigns cell types to transcripts detected in every spot according to the signatures extracted in the previous script. 

* [07c_Cell2Loc_Deconvolution_Downstream_Analysis-KUL3_cohort.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/InternalST_KUL3/07c_Cell2Loc_Deconvolution_Downstream_Analysis-KUL3_cohort.ipynb): downstream analysis on the deconvolution resutls containing a clustering based on cell abundance and a cell type co-localization analysis, among others. 


