# Deconvolution on the ST liver metastasis samples using the SMC cohort

This folder contains the deconvolution results on the external set of CRC ST samples presented in this [study](https://aacrjournals.org/cancerdiscovery/article/12/1/134/675646/Spatiotemporal-Immune-Landscape-of-Colorectal) using as reference the annotations from the scRNA-seq generated from a Korean cohort (SMC dataset). It is organized as follows: 

* [07a_CEll2Loc_ExpressionSignature.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/InternalST_SMC/07a_CEll2Loc_ExpressionSignature.ipynb): extraccion of cell type specific signatures according to the annotations from the reference. Of note, this script is common for for the internal and external datasets. 

* [07b_Cell2Loc_Deconvolution-CRC_LiverMetastasis.ipynb](https://github.com/alberto-valdeolivas/CRC_CMS_ST/blob/main/Deconvolution/ExternalST_SMC/07b_Cell2Loc_Deconvolution-CRC_LiverMetastasis.ipynb): the  deconvolution itself. Cell2Location assigns cell types to transcripts detected in every spot according to the signatures extracted in the previous script. 

