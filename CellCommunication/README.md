# Cell Communication events at the tumor-stroma interface of CMS2 carcinomas

This folder contains the scripts created to study cell communication events at the tumor-stroma interface in CMS2 carcinomas and their potential associations with tumor progression processes.


* [14_InterPatient_Hetereogeneity_Dorothea_profiles.Rmd](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/14_InterPatient_Hetereogeneity_Dorothea_profiles.Rmd): This script cluster the spots of the samples containing CMS2 tumors using their TF activity profiles. Among others, it contains the code to generate the Figure 4a-c and the Supplementary Figures S32-S33 and S35. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/14_InterPatient_Hetereogeneity_Dorothea_profiles.html?download=1)
* [14_InterPatient_Hetereogeneity_Dorothea_profiles_PathoAnnotations](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/14_InterPatient_Hetereogeneity_Dorothea_profiles_PathoAnnotations.Rmd): This scripts explore the pathological annotations of the spots clustered under different groups when using their TF activity profiles. Among others, it contains the code to generate the Figure 4d and the Supplementary Figure S34. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/14_InterPatient_Hetereogeneity_Dorothea_profiles_PathoAnnotations.html?download=1)
* [22_Misty_TFsDE_ligandsDE](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/22_Misty_TFsDE_ligandsDE.Rmd): This script individually runs Misty in every sample. It tries to model the activity of the TFs operating in the TME using the expression of ligands expressed in the tumor and the TME. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/22_Misty_TFsDE_ligandsDE.html?download=1)
* [22_Misty_TFsDE_ligandsDE_AggregatedResults](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/22_Misty_TFsDE_ligandsDE_AggregatedResults.Rmd): This scripts aggregates the individual Misty results in every sample to find the most consistent signals across samples. It contains the code for the Figure 4e,  Figure 5e-f, Figure 5h-i and several supplementary Figures. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/22_Misty_TFsDE_ligandsDE_AggregatedResults.html?download=1)


In particular, to explore ligand-receptor interactions, these are the required scripts: 

* [14_LigandReceptor_DorotheaClusters](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/14_LigandReceptor_DorotheaClusters.Rmd): This script contains the code to run LIANA in our CMS2 ST samples when grouping the spots by the clustering based on the TF activities. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/14_LigandReceptor_DorotheaClusters.html?download=1)
* [14_LigandReceptor_DorotheaClusters_ExploringResults](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/14_LigandReceptor_DorotheaClusters_ExploringResults.Rmd): This script reads the results generated in the previous script and generates the figure 4f. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/14_LigandReceptor_DorotheaClusters_ExploringResults.html?download=1)
* [24_Liana_sc_CMS2](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/24_Liana_sc_CMS2.Rmd): This script contains the code to run LIANA in the patients labeled as CMS2 in the scRNA-seq dataset from Lee et al. It contains the code to generate Figure 5 b-d. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/24_Liana_sc_CMS2.html?download=1).

To explore TF activity in the different cell tyes from the scRNA-seq data from Lee et al: 

* [25_TFactivity_sc_CMS2](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/25_TFactivity_sc_CMS2.Rmd): This script computes the activity of the TFs operating in the TME in the scRNA-seq dataset generated by Lee et al. It contains the code to generate Figure 5a. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/25_TFactivity_sc_CMS2.html?download=1)

To infer the most likely signaling cascades between the predicted ligand-receptor interaction and their corresponding TFs:

* [23_graphGeneration_TopLigandReceptors](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/CellCommunication/23_graphGeneration_TopLigandReceptors.Rmd): This script connects the predicted ligands-TF associations by using a network-based appraoch. In particular, it explores the shortest paths between the potential receptors of the predicted ligands and the associated TFs. It contains the code to generate the Figure 4g. The HTML notebook is available at zenodo and can be downloaded [here](https://zenodo.org/record/7588156/files/23_graphGeneration_TopLigandReceptors.html?download=1)



