# Charting the heterogeneity of colorectal cancer consensus molecular subtypes using Spatial transcriptomics

We processed fresh-frozen resection samples obtained from seven CRC patients for ST using 10x Genomics VISIUM aiming at exploring spatial molecular heterogeneity in CRC. We considered two serial sections per patient to generate technical replicates.The study outline is summarized in the following figure: 

![Study Outline](https://github.com/alberto-valdeolivas/CRC_CMS_ST/raw/main/Extras/StudyOutline.png)

This repository contains the scripts developed for the analysis of this set of CRC samples which results are presented in the following publication: 

**"Charting the Heterogeneity of Colorectal Cancer Consensus Molecular Subtypes using Spatial Transcriptomics"**   
_Alberto Valdeolivas, Bettina Amberg, Nicolas Giroud, Marion Richardson, Eric J.C. Gálvez, Solveig Badillo, Alice Julien-Laferrière, Demeter Turos, Lena Voith von Voithenberg, Isabelle Wells, Amy A. Lo, Emilio Yángüez, Meghna Das Thakur, Michael Bscheider, Marc Sultan, Nadine Kumpesa, Björn Jacobsen, Tobias Bergauer, Julio Saez-Rodriguez, Sven Rottenberg, Petra C. Schwalie, Kerstin Hahn  
bioRxiv, 2023; doi: https://doi.org/10.1101/2023.01.23.525135_

and it is organized as follows:

* [Quality Control](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/Quality_Control): Scripts to evaluate the global QC metrics compared across samples and batches. We also compare the number of reads in spots covered and non-covered by tissue to investigate potential contamination. 
* [Deconvolution](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/Deconvolution): Scripts to perform the deconvolution in the different ST datasets and using the two different available scRNA-seq references (Korean and Belgian cohort).
* [Sample Characterization](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/Sample_Characterization): Scripts to characterize the internal and external ST CRC datasets used in our manuscript based on different criteria. Among others, it contains the scripts required to generate the barplots describing the different cell type proportions per sample shown in the figures of the main manuscript. 
* [General Analysis per Sample](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/General_Analysis_perSample): Scripts to perform analysis, such as clustering, pathway and TF activity, individually on each sample individually. 
* [General Integrated analysis](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/General_Integrated_Analysis): Scripts to perform genenal analysis where the information from all the samples was merged in order to obain the results. This includes, Batch correction and UMAP embeddings, enrichment/depletion plots, Pathway and TF correlations with cell type abundances and some further validations of the deconvolution results. 
* [Inter-patient heterogeneity in CMS2 tumors](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/InterPatient_Heterogeneity): Scripts to explore the differences between CMS2 tumors originating from different patietns in our dataset. It contains the code for the analysis performed in the results section 2.3 of the manuscript. 
* [Intra-patient heterogeneity in CMS2 tumors](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/IntraPatient_Heterogeneity): Scripts to explore the intra-tumor heterogeneity of two different CMS2 carcinomas from two CRC patients. It contains the code for the analysis performed in the results section 2.3 of the manuscript.   
* [Cell-Communication events at the tumor-stroma interface of CMS2 carcinomas](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/CellCommunication): Scripts to explore cell communication events at the tumor-stroma interface of CMS2 carcinomas. In particular, it contains the scripts to perform the clustering of our CMS2 tumor samples based on TF activity. It also contains the scripts to infer which ligands can potentially have an influence in the activity of TFs active in tumor-surrounding regions. Then, it also contains the scipts to explore ligand-receptor interactions in out ST dataste and in an external scRNAseq dataset. Finally, it contains the script to infer the most likely signaling cascades connecting these inter-cellular communication events. These are basically the scripts needed for the results section 2.4 and figures 4 and 5. 
* [CancerSEA module Scores](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/CancerSEA_Scores): It contains the scripts to compute the module score of different cancaer related functions according to the annotations of [CancerSEA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6324047/). 
* [Additional analysis in the External ST CRC Dataset](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/AdditionalAnalysis_External_ST_CRC_Dataset):It contains the scripts with all the additional analysis performed on the external ST CRC dataset. 
* [Figure arrangement](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/FigureArrangement): Scripts containing the code to generate the figures displaying the spatial maps of differnet molecular features across our different set of samples (Gene expression, cell abundance maps as estimated by the deconvolution, pathway and TF activity...)
* [Wrapper Functions](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/WrapperFunction): It contains some functions and wrappers of Seurat functions used in the other scripts. 

The raw data (count matrices, H&E figures...) and the intermediates files generated during the analysis and required to re-run the analysis are available at zenodo: 

[https://doi.org/10.5281/zenodo.7551712](https://doi.org/10.5281/zenodo.7551712)

The notebooks associated with the scripts deposit in this repository are also available at zenodo: 

[https://doi.org/10.5281/zenodo.7440182](https://doi.org/10.5281/zenodo.7440182)







