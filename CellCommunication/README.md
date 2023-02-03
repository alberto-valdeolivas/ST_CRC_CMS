# Cell Communication events at the tumor-stroma interface of CMS2 carcinomas

This folder contains the scripts created to study cell communication events at the tumor-stroma interface in CMS2 carcinomas and their potential associations with tumor progression processes.


* [14_InterPatient_Hetereogeneity_Dorothea_profiles.Rmd](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/14_InterPatient_Hetereogeneity_Dorothea_profiles.Rmd): This script cluster the spots of the samples containing CMS2 tumors using their TF activity profiles. Among others, it contains the code to generate the Figure 4a-c and the Supplementary Figures S32-S33 and S35.
* [14_InterPatient_Hetereogeneity_Dorothea_profiles_PathoAnnotations](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/14_InterPatient_Hetereogeneity_Dorothea_profiles_PathoAnnotations.Rmd): This scripts explore the pathological annotations of the spots clustered under different groups when using their TF activity profiles. Among others, it contains the code to generate the Figure 4d and the Supplementary Figure S34. 
*[22_Misty_TFsDE_ligandsDE](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/22_Misty_TFsDE_ligandsDE.Rmd): This script individually runs Misty in every sample. It tries to model the activity of the TFs operating in the TME using the expression of ligands expressed in the tumor and the TME. 
*[22_Misty_TFsDE_ligandsDE_AggregatedResults](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/22_Misty_TFsDE_ligandsDE_AggregatedResults.Rmd): This scripts aggregates the individual Misty results in every sample to find the most consistent signals across samples. It contains the code for the Figure 4e,  Figure 5e-f, Figure 5h-i and several supplementary Figures. 






