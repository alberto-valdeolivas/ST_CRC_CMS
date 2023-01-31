# General Integrated Analysis. 

This folder contains the scripts created to perform different analyses jointly considering all our ST samples, excepting the evaluation of the deconvolution results where only replicates where jointly analysed. In particular, these are the details about the scripts stored in this folder:

* [Batch Correction and clustering](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/General_Integrated_Analysis/Pub_03_Filtering_Normalization_JointClustering_CorrectSample.Rmd): In this script, we used Harmony to correct batch effect per sample. Among others, it contains the code to generate the UMAP shown in Figure 1b in the main manuscript. The HTML notebook is available at zenodo and can be downloaded here.
* [Deconvolution Results Versus Annotations](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/General_Integrated_Analysis/ProportionsCell2Loc_VS_Annotations.Rmd): In this script, we evaluated in which tissue regions (annotated by a pathologist) are located the different cell type proportions estimated by the deconvolution. Among others, it contains the code to generate the Supplementary figure S2. The HTML notebook is available at zenodo and can be downloaded here.
* [Deconvolution evaluation by comparison of results between replicates](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/General_Integrated_Analysis/Pub_XX_ClusteringReplicates_EvalDeconv_V2.Rmd): In this script, the results of the deconvolution were evaluated by comparing results between replicates. In particular, joint clustering between replicates was conducted with different resolution parameters in order to select small anatomical regions that can be considered equivalent between replicates. Then, we computed Pearson's correlation per cell type in those regions. Among other, this script contains the code to generate Figure 1d. The HTML notebook is available at zenodo and can be downloaded here.
* [Enrichment Depletion Plots](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/General_Integrated_Analysis/Enrichment_DepletionPlots.Rmd): This script contains the code to generate the Figures 1e and 2f. The HTML notebook is available at zenodo and can be downloaded here.
* [Correlation between Pathway activity and Cell abundance](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/General_Integrated_Analysis/Pub_XX_correlation_pathway_CMS.Rmd): This script computes the per spot correlation between the estimated pathway activities and the estimated cell abundances. Among others, it contains the code to generate Figure 2l. The HTML notebook is available at zenodo and can be downloaded here.
* [Correlation between TF activity and Cell abundance](https://github.com/alberto-valdeolivas/ST_CRC_CMS/blob/main/General_Integrated_Analysis/Pub_XX_correlation_tfs_CMS.Rmd): This script computes the per spot correlation between the estimated TF activities and the estimated cell abundances. Among others, it contains the code to generate Figure 2k. The HTML notebook is available at zenodo and can be downloaded here.









