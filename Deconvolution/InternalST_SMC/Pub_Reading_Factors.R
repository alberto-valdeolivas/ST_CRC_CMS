#### Factors plus correlations Duplicates. 

library(readr)
library(pheatmap)
library(RColorBrewer)
setwd("/pstore/data/biomics/_pre_portfolio/_platform_evaluation/SpatialTranscriptomics/Experiment_CRC_AllSamples/analysis_alberto/")
alberto <- read_csv("Cell2Location/results/colocalizationFactors/R3.csv") %>% 
  dplyr::rename(CellType = "...1") %>% 
  dplyr::filter(CellType != "Unknown")

alberto <- read_csv("Cell2Location/results/colocalizationFactors/R3.csv") %>% 
  dplyr::rename(CellType = "...1") %>% 
  dplyr::filter(CellType != "Unknown")


celltypesTOchange <- data.frame(
  CellType = alberto$CellType,
  CellTypeNew = c("CD19+CD20+ B", "CD4+ T cells", "CD8+ T cells", "CMS1",
                  "CMS2", "CMS3", "CMS4", "Enteric glial",
                  "Goblet cells", "IgA+ Plasma", "IgG+ Plasma", "Intermediate",
                  "Lymphatic ECs", "Mast cells", "Enterocytes 1", "Enterocytes 2",
                  "Myofibroblasts", "NK cells", "Pericytes", "Pro-inflammatory",
                  "Proliferating", "Proliferative ECs", "Tregs", "SPP1+",
                  "Smooth muscle", "Stalk-like ECs", "Stem-like/TA", "Stromal 1",
                  "Stromal 2", "Stromal 3", "T foll helper", "T helper 17",
                  "Tip-like ECs", "cDC", "GammaDelta T cells"))


matrix_to_plot <- 
  inner_join(alberto, celltypesTOchange) %>% 
  dplyr::select(-CellType) %>% 
  dplyr::rename(CellType = CellTypeNew) %>%
  tidyr::pivot_longer(!CellType,  names_to = "Factor", values_to = "value") 
  
order <- matrix_to_plot %>%
  tidyr::pivot_wider(names_from = Factor, values_from = value) %>%
  tibble::column_to_rownames(var="CellType")

hc <- hclust(dist(order))
# plot(hc)
  

matrix_to_plot$CellType <- factor(matrix_to_plot$CellType, levels=hc$labels[hc$order])
  
# mypalette <- c("#FFFFFF",brewer.pal(n = 5, name ="Reds"))
  
# colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(42)

# par(mar=c(10,10,10,10))
    
# p <- pheatmap(t(matrix_to_plot), treeheight_row=0,treeheight_col=0, 
#         color =colorRampPalette(mypalette)(100),cluster_rows = FALSE,
#         fontsize = 18) 

p1 <- ggplot(matrix_to_plot, aes(CellType, Factor , fill= value)) + 
  geom_tile() + 
  theme(panel.background = element_blank(), 
  axis.text = element_text(size=12, face ="bold", family="Arial"), 
  axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2),
  axis.title = element_blank()) + 
  scale_fill_gradient(high = "#A50F15", low= "#FFFFFF")  + 
  theme(legend.position = "top", 
        legend.text = element_text(angle=270, size=12, face="bold", family="Arial", hjust=0.9,vjust=0),
        legend.title = element_blank())
  
  

df_correlations_all <- readRDS("/pstore/data/biomics/_pre_portfolio/_platform_evaluation/SpatialTranscriptomics/Experiment_CRC_AllSamples/analysis_alberto/IntermediaryFiles/Df_deconvolution_correlations.rds")

p2 <- df_correlations_all %>%
  ggplot( aes(x=patient_ID.y, y=correlations, fill=patient_ID.y)) +
  geom_boxplot(outlier.size = 0.75, alpha= 0.75) +
  # scale_fill_brewer(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.75, alpha=1) +
  theme_linedraw() +  
  theme(legend.position="none", axis.title.x = element_blank(), 
        axis.text = element_text(size=12, face ="bold", family="Arial"),
        axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2), 
        axis.title.y = element_text(size=16,family="Arial"),
        axis.ticks = element_blank(),  panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ylab("Correlation")

p2 + p1 +  plot_layout(widths = c(1, 2))
