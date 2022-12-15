ComputeClustering <- function(seurat_object, resolution_parameter){
  
  
  seurat_object <- seurat_object %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(dims = 1:30, k.param = 10, reduction = "harmony") %>%
    FindClusters(verbose = FALSE, resolution = resolution_parameter) %>% 
    identity()
  
}

ComputeCorrelation <- function(seurat_object, results_C2L){
  all_clusters <- levels(Seurat::Idents(seurat_object))
  
  spots_cluster_info <- 
    dplyr::select(seurat_object@meta.data, orig.ident, seurat_clusters) %>%
    tibble::rownames_to_column(var = "spot_id") %>%
    dplyr::mutate(spot_id = str_remove(.$spot_id, pattern = "_[1-2]")) %>%
    dplyr::mutate(spot_id = paste0(.$orig.ident, "_", .$spot_id))
  
  matching_cluster_c2l <- 
    dplyr::inner_join(results_C2L, spots_cluster_info) %>%
    dplyr::relocate(orig.ident, seurat_clusters, .after=spot_id)
  
  length_df <- length(cell_types) * length(all_clusters)
  df_results <- data.frame(cluster_id = character(length = length_df), 
                           cell_type= character(length = length_df),
                           mean_rep1 = numeric(length = length_df),
                           mean_rep2 = numeric(length = length_df))
  
  a <- 1
  
  for (current_cluster in all_clusters){
    

    for (current_cellType in cell_types){
      
      mean_UMI_rep1 <- matching_cluster_c2l %>% 
        dplyr::filter(orig.ident == rep1_name, seurat_clusters == current_cluster) %>% 
        dplyr::pull(current_cellType) %>% mean()
      
      mean_UMI_rep2 <- matching_cluster_c2l %>% 
        dplyr::filter(orig.ident == rep2_name, seurat_clusters == current_cluster) %>% 
        dplyr::pull(current_cellType) %>% mean()
      
      df_results$cluster_id[a] <-current_cluster
      df_results$cell_type[a] <- current_cellType
      df_results$mean_rep1[a] <- mean_UMI_rep1 
      df_results$mean_rep2[a] <- mean_UMI_rep2
      
      a <- a + 1
    }
    
  }
  
  return(cor(df_results$mean_rep1, df_results$mean_rep2, use = "pairwise.complete.obs"))
  
}