# Copied from Alberto V.,Adapted some by AJL 


####################### QC 
#' get.seurat.objects
#' read multiples spaceranger output folder for subsequent analysis
#'
#' @param datasets : a vector of datasets, 
#' should be paths to the outs of spaceranger; each folder is expected to contain either  raw_feature_bc_matrix.h5 ; filtered_feature_bc_matrix.h5
#' @param sample_names : vector of sample names; will be the names of the returned list.
#' @param raw : boolean; whether the raw or the filtered matrix shold be read (filtered by spots under tissues)
#'
#' @return a list of seurat objects, one object per slice (per dataset)
#' @export
#'
#' @examples
get.seurat.objects <- 
  function(datasets, sample_names, raw = FALSE){
  ## This is because Load10X_Spatial with filter.matrix = TRUE do not filter fully the spots weirdly
  # In his casxe, spatial representation is good, but not the violin plot for the features.
  filename <- ifelse( raw, 'raw_feature_bc_matrix.h5', 'filtered_feature_bc_matrix.h5')
  seurat_objects <- datasets %>% 
    map(function(x){
      slice_name <- gsub(".*Count_","",x)
      slice_name <- gsub("/outs","",slice_name)
      Load10X_Spatial(x, filename = filename, filter.matrix = ! raw, 
                      slice = slice_name) 
    }) 
  names(seurat_objects) <- sample_names
  for (i in seq_along(sample_names)){
    seurat_objects[[i]][["orig.ident"]] <- sample_names[i]
    Project(seurat_objects[[i]]) <- sample_names[i]
  }
  return(seurat_objects)
}


#' add_tissue_position
#' read information of tissue from the space ranger output and return an updated metadata file
#' 
#' @param oneSeuratObject : seurat object with initial metadata
#' @param dataset : where spacerange output is stored (should contain spatial/tissue_positions_list.csv)
#'
#' @return an updated metadata table
#' @export
#'
#' @examples
add_tissue_position <- function( oneSeuratObject, dataset){
  # info_pos  <- read.csv( paste0( dataset, "/spatial/tissue_positions_list.csv"), header = FALSE,
  #                        col.names = c( "Barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")) %>%
  #  tibble::column_to_rownames("Barcode")
  
  ## Modified for the weird folder structure of 
  
  info_pos <- 
    read.csv( paste0(dataset, "/spatial/", "tissue_positions_list.csv"),
      header = FALSE, col.names = c( "Barcode", "in_tissue", "array_row", 
      "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")) %>%
    tibble::column_to_rownames("Barcode")
  
  meta.data <- merge( oneSeuratObject@meta.data, info_pos, by = 0)
  return(meta.data) 
}




#' get.image.spatial
#' returns a list of figures containing the raw image (H&E)
#' 
#'
#' @param list_seurat : list of seurat objects
#' @param crop : see Seurat::SpatialFeaturePlot
#'
#' @return  list of images (patchwork/ggplot)
#' @export
#'
#' @examples
get.image.spatial <- function(list_seurat, crop=FALSE){
  
  image_plots <- 
    lapply(list_seurat, function(x){
      Seurat::SpatialFeaturePlot(x, 
                                 features = c("nCount_Spatial"), 
                                 pt.size.factor = 0, crop = crop) + 
        labs(title = names(x)) 
    })
  
  return(image_plots)
}



#' get.qc.spatial
#'
#' @param list_seurat list of seurat objects
#' @param features default c("nCount_Spatial", "nFeature_Spatial"), see Seurat::SpatialFeaturePlot
#' @param pt.size.factor see  Seurat::SpatialFeaturePlot
#' @param alpha see  Seurat::SpatialFeaturePlot
#'
#' @return list of Seurat::SpatialFeaturePlot (one per Seurat Objects); 
#' @export
#'
#' @examples
get.qc.spatial <- function(list_seurat, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size.factor = 1.6, alpha = 1){
  
  spatial_qc_plots <- 
    lapply(list_seurat, function(x){
      Seurat::SpatialFeaturePlot(x, 
                                 features = features , 
                                 alpha = alpha, 
                                 pt.size.factor = pt.size.factor) + 
        labs(title = names(x)) 
    })
  
  return(spatial_qc_plots)
}


#' get.qc.violin
#'
#' @param list_seurat 
#' @param features default c("nCount_Spatial", "nFeature_Spatial"), see Seurat::VlnPlot
#' @param split_by see Seurat::VlnPlot
#'
#' @return list of Seurat::VlnPlot
#' @export
#'
#' @examples
get.qc.violin <- function(list_seurat, 
                          features = c("nCount_Spatial", "nFeature_Spatial"),
                          split_by = FALSE){
  
  violin_qc_plots <- 
    lapply(list_seurat, function(x){
      Seurat::VlnPlot(x, 
                      features = features) + 
        labs(title = names(x),
             split_by = split_by) + NoLegend()
    })
  
  return(violin_qc_plots)
}


#' get.undertissue.violon
#' return violin plot from Seurat object CONTAINING the image metadata 
#' violin plot will then be splited by spot condition (under tissue/not under tissue)
#' (see add_tissue_position function to obtain the right info.)
#'
#' @param list_metadata 
#'
#' @return
#' @export
#'
#' @examples
get.undertissue.violon <- function( list_metadata){
  raw.qc.plots <- lapply(list_metadata, function(x  ){
    sample <- x %>%
      mutate( in_tissue = factor(plyr::revalue(as.character(in_tissue), c( "0"="noTissue", "1"="underTissue")) ) ) %>%
      dplyr::select(nCount_Spatial,  nFeature_Spatial, in_tissue) %>% 
      pivot_longer(., cols = c("nCount_Spatial",  "nFeature_Spatial"))
    
    ggplot(sample, aes( x = in_tissue, y = value, fill = in_tissue) ) + 
      geom_violin()  + geom_jitter(size= 0.5, alpha = 0.5)+
      facet_wrap(~name, scales =  'free_y')
  } )
  
  
  return(raw.qc.plots)
}


#' get.qc.mt.percentage
#' This function computes the percentage of mitochondrial genes per sample 
#'
#' @param list_seurat list containing seurat objects
#' @param mypattern Pattern of the mitochondrial genes. It may change by species.
#'
#' @return list of Seurat::SpatialFeaturePlot 
#' @export
#'
#' @examples

get.qc.mt.percentage <- function(list_seurat, mypattern = "^MT-",...){
  
  
  seurat_objects <- lapply(list_seurat, function(x){
    percentage_mt <- PercentageFeatureSet(x, pattern = mypattern)
    AddMetaData(x, percentage_mt, col.name = "percent.mt")
  })
    

  # spatial_mt_plots <- 
  #  lapply(list_seurat, function(x){
  #    
  #    count_matrix <- GetAssayData(x, slot = "counts", assay = "Spatial")
  #    mito_genes <- 
  #      rownames(count_matrix)[startsWith(rownames(count_matrix), 'Mt-')]
  #    
  #    x@meta.data$percentage_mt <- 
  #      round(colSums(count_matrix[mito_genes,])/x@meta.data$nCount_Spatial * 100,2)
  #    
  #    Seurat::SpatialFeaturePlot(x, 
  #        features = "percentage_mt", 
  #        ...) + 
  #      labs(title = names(x)) 
  #  })
  
  return(seurat_objects)
}



#' get.FeatureScatter.plots
#' This function gets Feature Scatter plots from the Seurat package for a list
#' of Seurat objects and features. 
#' 
#'
#' @param list_seurat list containing seurat objects
#' @param ... Additional parameters for the Seurat::FeatureScatter
#'
#' @return list of Seurat::FeatureScatter 
#' @export
#'
#' @examples

get.FeatureScatter.plots <- function(list_seurat, ...){
  FeatureScatter_plots <- lapply(list_seurat, function(x){
    FeatureScatter(x, ...) + NoLegend()
  })
  return(FeatureScatter_plots)
}


#################### CLUSTERING



normalize.cluster <- function(list_seurat, k.param = 20, resolution = 0.8,
  min.dist = 0.3, spread = 1){
  
  seurat_objects <- list_seurat %>% 
    map(function(x) {
      SCTransform(x, assay = "Spatial", verbose = FALSE) %>% 
        RunPCA(assay = "SCT", verbose = FALSE) %>% 
        FindNeighbors(reduction = "pca", dims = 1:30, k.param= k.param) %>% 
        FindClusters(verbose = FALSE, resolution = resolution)   %>% 
        RunUMAP(reduction = "pca", dims = 1:30, 
                min.dist = min.dist, spread = spread)  
    })
  
  
  return(seurat_objects)
} 

get.umap.spatialClusters <- function(list_seurat, label = TRUE, ...){
  
  umap_spatial_plots <- lapply(list_seurat, function(x){
    umap_plot <- Seurat::DimPlot(x, reduction = "umap", label = TRUE)
    spatial_plot <- Seurat::SpatialDimPlot(x, label = TRUE, ...)
    umap_plot + spatial_plot
  })
  
  return(umap_spatial_plots)
}


get.all.markers <- function(list_seurat, mytreshold = 0.01,...){
  
  all_markers <- lapply(list_seurat, function(x){
    Seurat::FindAllMarkers(object = x, ...) %>% 
      dplyr::filter(p_val_adj  < mytreshold)
  }) 
  
  return(all_markers)
}


get.markers.table <- function(markers_list, nr_markers = 2){
  all_markers_table <- lapply(markers_list, function(x){
    x %>% 
      as_tibble() %>% 
      group_by(cluster) %>% 
      top_n(n = nr_markers, wt = avg_log2FC ) %>% 
      as.data.frame(col.names = colnames(.)) 
  })
}


get.enrichment.results <- function(markers_list, background_genes, ...){
  
  
  #all_enrichment_results <- 
  #  mapply(function(x,y, SIMPLIFY = FALSE) {
  #    clusters <- dplyr::pull(x, cluster) %>% unique()
  #    query_results_df <- data.frame()
  #    
  #    for (i in clusters){
  #      query <- dplyr::filter(x, cluster  == i) %>% pull(gene)
  #      query_results <- 
  #        gost(query, custom_bg = y, ...)
  #      query_results_df <- rbind.data.frame(query_results_df, 
  #        query_results$result %>% dplyr::mutate(cluster = i))  %>% 
  #        as_tibble()
  #    }
  #    
  #    query_results_df_final <- query_results_df  %>% 
  #      group_by(cluster) %>% nest() %>% as_tibble()
  #  }, x=markers_list, y=background_genes)
  
  counter <- 0
  
  all_enrichment_results <- lapply(markers_list, function(x){
    
    clusters <- dplyr::pull(x, cluster) %>% unique()
    query_results_df <- data.frame()
    
    counter <<- counter + 1
    
    for (i in clusters){
      query <- dplyr::filter(x, cluster  == i) %>% pull(gene)
      query_results <- 
        gost(query, custom_bg = unlist(background_genes[counter], use.names = FALSE)) # , ...)
      query_results_df <- rbind.data.frame(query_results_df, 
                                           query_results$result %>% dplyr::mutate(cluster = i))  
    }
    
    query_results_df_final <- query_results_df  %>% 
      group_by(cluster) %>% nest() %>% as_tibble()
    
  })
  
  return(all_enrichment_results)
}


get.spatial.features <- function(list_seurat, nr_variable_features = 1000, 
  method = "markvariogram", dorothea = FALSE, mito = FALSE){ 
  
  seurat_objects <- lapply(list_seurat, function(x){
    
    if (dorothea) {
      
      Tfs <- unique(rownames(GetAssayData(x, assay = "dorothea")))
      FindSpatiallyVariableFeatures(x, assay = "dorothea", 
                                    features = Tfs, slot = "data",
                                    selection.method = method)
    } else {
      
      if (mito == TRUE){
        features_selected <-  VariableFeatures(x)[1:nr_variable_features]
      } else {
        idx <- grep("^MT[-]",VariableFeatures(x))
        features_selected <-  VariableFeatures(x)[-idx]
        features_selected <- features_selected[1:nr_variable_features]
          
      }  
        
      FindSpatiallyVariableFeatures(x, assay = "SCT", 
          features = features_selected, selection.method = method)
    }
    
  })
  
  return(seurat_objects)
  
}

write.spatial.features <- 
  function(list_seurat, printname, assay = "SCT", nr_features = 1000, 
           method = "markvariogram"){
  
  lapply(list_seurat, function(x){
    top.features <- 
      SpatiallyVariableFeatures(x, selection.method = method, assay = assay)
    top.features.print <- top.features[1:nr_features]
    write_csv(as.data.frame(top.features.print), 
              file = paste0(printname,"_", names(x), ".csv"))
    
  })  
      
}

get.spatial.plots <- function(list_seurat, method = "markvariogram", 
  nr_features = 6, alpha = c(1, 1), ncol = 3, pt.size.factor = 1.6, 
  assay = "SCT", dorothea = FALSE, ...){
  
  variable_spatial_plots <- lapply(list_seurat, function(x){
    
    top.features <- SpatiallyVariableFeatures(x, selection.method = method,
      assay = assay)
    
    top.features.plot <-head(top.features, nr_features)
    
    if (dorothea){
      DefaultAssay(x) <- "dorothea"
      SpatialFeaturePlot(x, features = top.features.plot, ncol = ncol, 
                         alpha = alpha,
                         pt.size.factor = pt.size.factor, ...)  
    } else {
      SpatialFeaturePlot(x, features = top.features.plot, ncol = ncol, 
                         alpha = alpha,
                         pt.size.factor = pt.size.factor, ...)  
    }
    
  })
  
  return(variable_spatial_plots)
} 

get.progeny.scores <- 
  function(list_seurat,scale=FALSE, organism="Human", top=500, perm=1,
           return_assay = TRUE, assay_name = "SCT"){
    
    
    seurat_objects <- lapply(list_seurat, function(x) {
      seurat_objects <- 
        progeny(x, scale=scale, organism=organism, top=top, perm=perm,
                return_assay = TRUE, assay_name = assay_name)
    })
    
    seurat_objects <- lapply(seurat_objects, function(x) {
      seurat_objects <- Seurat::ScaleData(x, assay = "progeny")
      
      
    })
    
    return(seurat_objects)
  }  


get.pathway.plots <- function(list_seurat, ncol = NULL, pt.size.factor = 1.6, ...){
  
  pathway_plots <- lapply(list_seurat, function(x){
    pathways <- unique(rownames(GetAssayData(x, assay = "progeny")))
    DefaultAssay(x) <- "progeny"
    SpatialFeaturePlot(object = x, features = pathways, ncol = ncol, 
                       pt.size.factor = pt.size.factor, ...)
    
  })
  
}

get.dorothea.scores <- 
  function(list_seurat, confidence_levels=c("A","B","C"), organism="Human", 
           assay_name = "SCT", minsize = 4, eset.filter = FALSE, cores = 1){
    
    
    if (organism == "Human"){
      dorothea_regulon <- 
        get(data("dorothea_hs", package = "dorothea"))  
    } else {
      if (organism == "Mouse") {
        dorothea_regulon <- 
          get(data("dorothea_mm", package = "dorothea"))  
      } else {
        stop("Specie not supported")
      }
    }
    
    regulon_filtered <- dorothea_regulon %>%
      dplyr::filter(confidence %in% confidence_levels)
    
    seurat_objects <- lapply(list_seurat, function(x) {
      seurat_objects <- 
        run_viper(x, regulon_filtered, assay_key = assay_name, 
                  options = list(
                    method = "scale", minsize = minsize, 
                    eset.filter = eset.filter, cores = 1, verbose = FALSE))
    })
    
    return(seurat_objects)
  }  



get.enrichment.barplots <- 
  function(enrichment_df, nr_results_ontology = 10, threshold = 0.01) {
    
    clusters <- dplyr::pull(enrichment_df, cluster) %>% unique()
    
    plot_list <- list()
    
    for (i in clusters){
      Toplot_df <- enrichment_df %>% 
        dplyr::filter(cluster == i) %>% 
        unnest(cols = data) %>% 
        # dplyr::filter(p_value < threshold) %>% 
        group_by(source) %>% 
        dplyr::top_n(nr_results_ontology, desc(p_value)) %>% 
        ungroup() %>% 
        mutate(logPvalue = -log(p_value)) 
      
      plot_enrichment <- ggplot(Toplot_df, aes(logPvalue, 
                                               reorder(term_name, logPvalue))) +
        geom_col(aes(fill=logPvalue)) + facet_wrap(~ source) + 
        scale_fill_gradient(low = "#FFA07A", high = "#800000", 
                            na.value = "whitesmoke") +   
        theme_minimal() + 
        theme(legend.position = "bottom", 
              legend.text = element_text(size = 7, hjust = 1, angle = 90),
              axis.text.x = element_text(angle = 0, hjust = 0, size = 12, face="bold"), 
              axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) + 
        # ylab(paste0("Pathways and Ontologies")) +
        ylab(paste0("")) + 
        xlab(paste0("-Log (p-Value)")) + 
        labs(fill = "-Log (p-Value)") + 
        scale_x_continuous(position = "bottom") + 
        facet_grid(vars(source), scales = "free", space = "free") + 
        geom_vline(xintercept = -log(threshold), linetype = "dashed", 
                   colour = "steelblue", size = 1) + 
        ggtitle(paste0("cluster: ", i))
      
      
      plot_list[[i]] <- plot_enrichment
      
    }
    names(plot_list) <- clusters
    
    return(plot_list)
  }

get.global.QCmetrics <- function(dataset_vec, samples_name, metrics_name){
  
  read_files <- function(x){
    a <- read_csv(file = paste0(x, metrics_name))
    return(a)
  }
  return(dataset_vec %>% map_dfr(read_files) %>% add_column(samples_name))
  
}

get.barplot.qc <- function(qc_df, column_to_plot)  {
  
  df_to_plot <- qc_df %>% 
    dplyr::select(samples_name, all_of(column_to_plot)) %>% 
    dplyr::rename(nr_to_plot = column_to_plot) %>% 
    tibble::add_column(batch = str_remove(.$samples_name, pattern = "_.*")) %>% 
    dplyr::mutate(replicate = str_remove(.$samples_name, pattern = ".*_")) %>% 
    dplyr::mutate(sample_name_show = str_remove(.$samples_name, pattern = "SN[0-9]+_"))
  
  p <- ggplot(df_to_plot, aes(fill=batch, y=nr_to_plot , x=sample_name_show)) + 
    geom_bar(position="dodge", stat="identity") + theme_minimal() + 
    theme(legend.position = "top", legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12, hjust = 1, angle = 0),
          axis.text.x = element_text(angle = 90, hjust = 0, size = 12), 
          axis.text.y = element_text(angle = 0, hjust = 1, size = 12), 
          axis.title.y = element_text(size = 16, face = "bold")) + 
    xlab("") + ylab(column_to_plot) 
  
  return(p)
}

get.enrichment.results.cp <- function(markers_list, background_genes, ...){
  
  counter <- 0
  
  all_enrichment_results <- lapply(markers_list, function(x){
    
    clusters <- dplyr::pull(x, cluster) %>% unique()
    query_results_df <- data.frame()
    
    counter <<- counter + 1
    
    for (i in clusters){
      query <- dplyr::filter(x, cluster  == i) %>% pull(gene)
      # print(i)
      # print(counter)
      query_results <- 
        enricher(query, universe = unlist(background_genes[counter], 
                                          use.names = FALSE), ...) 
      
      if (!is.null(query_results)){
        query_results_df <- rbind.data.frame(query_results_df, 
                                             query_results@result %>% dplyr::mutate(cluster = i))  
      }
    }
    
    query_results_df_final <- query_results_df  %>% 
      group_by(cluster) %>% nest() %>% as_tibble()
    
  })
  
  return(all_enrichment_results)
}

get.enrichment.barplots.cp <- 
  function(enrichment_df, nr_results = 12, threshold = 0.01) {
    
    clusters <- dplyr::pull(enrichment_df, cluster) %>% unique()
    
    plot_list <- list()
    
    for (i in clusters){
      Toplot_df <- enrichment_df %>% 
        dplyr::filter(cluster == i) %>% 
        unnest(cols = data) %>% 
        # dplyr::filter(p_value < threshold) %>% 
        # group_by(source) %>% 
        dplyr::top_n(nr_results, desc(pvalue )) %>% 
        # ungroup() %>% 
        mutate(logPvalue = -log(pvalue)) %>%
        mutate(ID_toplot = substr(ID, 1,50))
      
      plot_enrichment <- ggplot(Toplot_df, aes(logPvalue, 
                                               reorder(ID_toplot, logPvalue))) +
        geom_col(aes(fill=logPvalue)) + # facet_wrap(~ source) + 
        scale_fill_gradient(low = "#FFA07A", high = "#800000", 
                            na.value = "whitesmoke") +   
        theme_minimal() + 
        theme(legend.position = "bottom", 
              legend.text = element_text(size = 8, hjust = 1, angle = 90),
              axis.text.x = element_text(angle = 0, hjust = 0, size = 8, face="bold"), 
              axis.text.y = element_text(angle = 0, hjust = 1, size = 8, face="bold")) + 
        # ylab(paste0("Pathways and Ontologies")) +
        ylab(paste0("")) + 
        xlab(paste0("-Log (p-Value)")) + 
        labs(fill = "-Log (p-Value)") + 
        scale_x_continuous(position = "bottom") + 
        # facet_grid(vars(source), scales = "free", space = "free") + 
        geom_vline(xintercept = -log(threshold), linetype = "dashed", 
                   colour = "steelblue", size = 1) + 
        ggtitle(paste0("cluster: ", i))
      
      
      plot_list[[i]] <- plot_enrichment
      
    }
    names(plot_list) <- clusters
    
    return(plot_list)
  }


#### Additional functions

write.all.markers <- function(markers_list, path){
  samples <- names(markers_list)
  for (i in samples){
    write_csv(markers_list[[i]], file = paste0(path, "markers_", i, ".csv"))
  }
  
}

write.loupe.file <- function(seurat_list, path){
  samples <- names(seurat_list)
  for (i in samples){
    df <- data.frame(barcode = rownames(seurat_list[[i]]@meta.data), 
      Cluster= seurat_list[[i]]@meta.data$seurat_clusters)
    write_csv(df, file = paste0(path, "clustersLoupe_", i, ".csv"))
  }
}  


######################### 

#+++++++++++++++++++++++++
# Computing of correlation matrix
#+++++++++++++++++++++++++
# Required package : corrplot
# x : matrix
# type: possible values are "lower" (default), "upper", "full" or "flatten";
#display lower or upper triangular of the matrix, full  or flatten matrix.
# graph : if TRUE, a correlogram or heatmap is plotted
# graphType : possible values are "correlogram" or "heatmap"
# col: colors to use for the correlogram
# ... : Further arguments to be passed to cor or cor.test function
# Result is a list including the following components :
# r : correlation matrix, p :  p-values
# sym : Symbolic number coding of the correlation matrix
rquery.cormat<-function(x,
                        type=c('lower', 'upper', 'full', 'flatten'),
                        graph=TRUE,
                        graphType=c("correlogram", "heatmap"),
                        col=NULL, ...)
{
  library(corrplot)
  # Helper functions
  #+++++++++++++++++
  # Compute the matrix of correlation p-values
  cor.pmat <- function(x, ...) {
    mat <- as.matrix(x)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # Get lower triangle of the matrix
  getLower.tri<-function(mat){
    upper<-mat
    upper[upper.tri(mat)]<-""
    mat<-as.data.frame(upper)
    mat
  }
  # Get upper triangle of the matrix
  getUpper.tri<-function(mat){
    lt<-mat
    lt[lower.tri(mat)]<-""
    mat<-as.data.frame(lt)
    mat
  }
  # Get flatten matrix
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  # Define color
  if (is.null(col)) {
    col <- colorRampPalette(
      c("#67001F", "#B2182B", "#D6604D", "#F4A582",
        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
        "#4393C3", "#2166AC", "#053061"))(200)
    col<-rev(col)
  }
  
  # Correlation matrix
  cormat<-signif(cor(x, use = "complete.obs", ...),2)
  pmat<-signif(cor.pmat(x, ...),2)
  # Reorder correlation matrix
  ord<-corrMatOrder(cormat, order="hclust")
  cormat<-cormat[ord, ord]
  pmat<-pmat[ord, ord]
  # Replace correlation coeff by symbols
  sym<-symnum(cormat, abbr.colnames=FALSE)
  # Correlogram
  if(graph & graphType[1]=="correlogram"){
    corrplot(cormat, type=ifelse(type[1]=="flatten", "lower", type[1]),
             tl.col="black", tl.srt=45,col=col,...)
  }
  else if(graphType[1]=="heatmap")
    heatmap(cormat, col=col, symm=TRUE)
  # Get lower/upper triangle
  if(type[1]=="lower"){
    cormat<-getLower.tri(cormat)
    pmat<-getLower.tri(pmat)
  }
  else if(type[1]=="upper"){
    cormat<-getUpper.tri(cormat)
    pmat<-getUpper.tri(pmat)
    sym=t(sym)
  }
  else if(type[1]=="flatten"){
    cormat<-flattenCorrMatrix(cormat, pmat)
    pmat=NULL
    sym=NULL
  }
  list(r=cormat, p=pmat, sym=sym)
}


merge_assay <- function(seurat_list, selected_features=NULL){
  
  merge_assay_data <- data.frame()
  if (is.null(selected_features)){
    for (i in names(seurat_list)){
      common_features <- 
        Reduce(intersect,lapply(seurat_list, function(x) rownames(GetAssayData(x))) )
      current_assay_data <- GetAssayData(seurat_list[[i]]) 
      current_assay_data <- 
        current_assay_data[rownames(current_assay_data) %in% common_features,] %>%
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "spot_id") %>%
        dplyr::mutate(spot_id = paste0(i, "_", .$spot_id)) 

      merge_assay_data <- rbind(merge_assay_data, current_assay_data)
    
    }
  } else {
    for (i in names(seurat_list)){
      current_assay_data <- GetAssayData(seurat_list[[i]]) 
      current_assay_data <- 
        current_assay_data[rownames(current_assay_data) %in% selected_features,] %>%
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "spot_id") %>%
        dplyr::mutate(spot_id = paste0(i, "_", .$spot_id)) 
      
      merge_assay_data <- rbind(merge_assay_data, current_assay_data)
    }  
    
  }
  
  return(merge_assay_data)
  
}

### Function to subset CMS2 sections.

subset_CMS2_sections <- function (seurat_list, resuts_decon, CMS2_cutoff = 0.5){
  lapply(seurat_list, function(x){
    current_sample <- x@project.name
    
    resuts_decon_current_sample <- resuts_decon %>%
      dplyr::filter(str_detect(spot_id, current_sample)) %>% 
      dplyr::mutate(spot_id = str_remove(spot_id,paste0(current_sample, "_")))
    
    x@meta.data <- x@meta.data %>% 
      tibble::rownames_to_column("spot_id") %>% 
      dplyr::left_join(resuts_decon_current_sample) %>% 
      tibble::column_to_rownames("spot_id")
    
    spots_toKeep <- x@meta.data %>% 
      dplyr::filter(CMS2 > CMS2_cutoff) %>%  rownames()
    
    DefaultAssay(x) <- "SCT"
    x <- subset(x, cells = spots_toKeep)
    names(x@images) <- current_sample
    x <- FindVariableFeatures(x)
    
    
  })
}