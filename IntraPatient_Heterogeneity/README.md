# Intra-patient heterogeneity in CMS2 tumors

This folder contains the scripts created to study the intra-patient heterogeneity in a couple of tumors classified as CMS2 in our set of samples. In particular, it contains the following scripts: 

* [34_GeneExpressionGradients_SampleA595688_Rep1](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/IntraPatient_Heterogeneity/34_GeneExpressionGradients_SampleA595688_Rep1.Rmd): This script extract the CMS2 tumor annotated spots from the S2_Col_R_Rep1 sample and assign them to different categories regarding their distance to non-tumor spots. Like this, we create a zonation model to distinguish different tumor regions (solid internal areas versus peripheral areas in contact with stromal regions). Among others, this script contains the code to generate the Figure 3g-h and the Supplementary Figures S26-S28 and Supplementary table 4. 
* [32_02_IntraPatient_BayesSpace_SampleA121573_Rep1](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/IntraPatient_Heterogeneity/32_02_IntraPatient_BayesSpace_SampleA121573_Rep1.Rmd): This script extract the CMS2 tumor annotated spots from the S5_Rec_Rep1 sample and performs sublustering on them using BayesSpace. In particular, it contains the code to generate the figure 3i and Supplementary Figures S29-S30 and the data from Supplementary table 5. 
* [32_03_IntraPatient_BayesSpace_SampleA121573_Rep1_PathwayActivity](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/IntraPatient_Heterogeneity/32_03_IntraPatient_BayesSpace_SampleA121573_Rep1_PathwayActivity.Rmd): This scripts uses BayesSpace to compute pathway activity at an enhanced resolution in the clusters resulting from the previous script. Among others, it contains the code to generate the Figure 3j and the Supplementary Figure S31.  







