get_sample_pseudo <- function(slide) {
  
  # Creates pseudobulk profile for each slide.
  bulk_p_data <- sumCountsAcrossCells(x = GetAssayData(slide, 
                                                       assay = "Spatial", slot = "counts"), ids = slide@meta.data$orig.ident)
  
  return(bulk_p_data)
}

#' Counts per million
#' @param expression_matrix: count matrix with samples in columns and genes in rows
#' @param scale_factor: used to generate the final counts
cpm_norm <- function(expression_matrix, scale_factor = 10000) {
  # First we transpose because of R conventions
  # Correct by coverage
  norm_mtx <- t(t(expression_matrix + 1)/colSums(expression_matrix))
  norm_mtx <- log1p(norm_mtx * scale_factor)
  return(norm_mtx)
}


#' Filter expression matrix
#' @param expression_matrix: count matrix with samples in columns and genes in rows
#' @param min.count: numeric. Minimum count required for at least some samples.
#' @param min.prop: numeric. Minimum proportion of samples in the smallest group that express the gene.
edgeR_filtering <- function(expression_matrix, 
                            min.count = 10,
                            min.prop = 0.1,
                            min.total.count = 15) {
  
  meta_data <- data.frame("sample_names" = colnames(expression_matrix),
                          "lib.size" = colSums(expression_matrix),
                          stringsAsFactors = F)
  
  #edgeR pipeline
  bulk_data <- edgeR::DGEList(expression_matrix, 
                              samples= meta_data)
  
  #Filtering lowly expressed genes in the context of the groups
  keep = edgeR::filterByExpr(bulk_data, 
                             min.count = min.count,
                             min.prop = min.prop,
                             min.total.count = min.total.count,
                             group = rep("same_group",ncol(bulk_data)))
  
  bulk_data = bulk_data[keep,]
  
  
  return(bulk_data$counts)
}