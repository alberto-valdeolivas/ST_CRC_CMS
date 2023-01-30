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
* [Sample Characterization](https://github.com/alberto-valdeolivas/ST_CRC_CMS/tree/main/Sample_Characterization): Scripts to characterize the internal and external ST CRC datasets used in our manuscript based on different criteria.
* [General Analysis per Sample]: Scripts to perform analysis, such as clustering, pathway and TF activity, individually on each sample. 
* General Integrated analysis
* Figure arrangement 
* Additional analaysis of the external ST CRC samples
* Wrapper Functions

The raw data (count matrices, H&E figures...) and the intermediates files generated during the analysis and required to re-run the analysis are available at zenodo: 

[https://doi.org/10.5281/zenodo.7551712](https://doi.org/10.5281/zenodo.7551712)

The notebooks associated with the scripts deposit in this repository are also available at zenodo: 

[https://doi.org/10.5281/zenodo.7440182](https://doi.org/10.5281/zenodo.7440182)







