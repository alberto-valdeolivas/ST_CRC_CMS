# Additional Analysis in the External ST CRC dataset

In this folder, we can find the scripts generated to analyse the external ST CRC dataset. In particular, these are the scripts you can find here: 

* [01_Clustering_Resolution_05](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/01_Clustering_Resolution_05.Rmd): Gene-expression based clustering of the spot of the different external ST CRC samples. 
* [02_PathwayActivity](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/02_PathwayActivity.Rmd): Per spot pathway activity and their associated spatial maps in the different external ST CRC samples. Among others, it contains the code to generate the figures 6g,i and some supplementary figures. 
* [03_Tf_Activity](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/03_Tf_Activity.Rmd): Per spot TF activity and their associated spatial maps in the different external ST CRC samples. Among others, it contains the code to generate the figures 6h,j and some supplementary figures. 
* [Pub_Proportion_Plots](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/Pub_Proportion_Plots.Rmd):This script contains the code to compute the per spot proportion of the different cell abundances. Among others, it contains the code to generate Figure 6a,b and the Supplementary Figures S49-S52.
* [Pub_Cell2LocPlots](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/Pub_Cell2LocPlots.Rmd): This script contains the code to generate the figures showing the overlay of the cell abundances (as estimated by the deconvolution) with the histological images. 
* [Pub_Pathway_CMS_correlation](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/Pub_Pathway_CMS_correlation.Rmd): This script contains the code to compute the per spot correlation between pathway activity and CMS cell type bundances. Among others, it cotains the code to generate Figure 6e. 
*[Pub_TF_CMS_correlation](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/Pub_TF_CMS_correlation.Rmd):This script contains the code to compute the per spot correlation between TF activity and CMS cell type bundances. Among others, it cotains the code to generate Figure 6f. 
* [07_Misty_TFs_DE_ligandsDE](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/07_Misty_TFs_DE_ligandsDE.Rmd): This script contains the code to generate misty predictions individually in each of the different external ST CRC samples. 
* [07_02_Misty_TFs_DE_Aggregated](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/07_02_Misty_TFs_DE_Aggregated.Rmd): This script collects and aggregates the indiviual results generated in the previous script in order to obtain the most consistent signals across samples. Among others, it contains the code to generate Figure 6k and Supplementary Figure S57. 
* [Pub_CMSCaller_PseudoBulkClassification](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/AdditionalAnalysis_External_ST_CRC_Dataset/Pub_CMSCaller_PseudoBulkClassification.Rmd): This script generates pseudo-bulk from the set of external ST CRC samples and tries to assing them to some CMS using CMScaller. It contains the code to generate the Supplementay Figure S53. 