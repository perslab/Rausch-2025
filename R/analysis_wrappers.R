#' Create edgeR object
#'
#' @title
#' @import edgeR
#' @return
#' @author dylanmr
#' @export
edger_prep <- function(seur, trt_group = "geno", lib.size = 1e4,
                       filter_column = NULL, filter_value = NULL,
                       assay_name = "RNA", feature_pattern_exclude = "^Gm|Rik$", 
                       dataset = NULL, celltypes = celltypes, celltype_column = NULL,
                       min.count = 1, min.total.count = 5) {
  
  # Validate inputs
  stopifnot("Seurat" %in% class(seur), assay_name %in% names(seur@assays))
  DefaultAssay(seur) <- "RNA"
  
  # Subset cells based on metadata column and filter features
  if(!is.null(celltype_column)){
    cells_to_keep <- colnames(seur)[grepl(paste0("^",celltypes, "$"), seur[[celltype_column]][,1])]
    features_to_keep <- rownames(seur@assays[[assay_name]])[!grepl(feature_pattern_exclude, rownames(seur@assays[[assay_name]]))]
    seur <- subset(seur, cells = cells_to_keep, features = features_to_keep)
  } else {
    features_to_keep <- rownames(seur@assays[[assay_name]])[!grepl(feature_pattern_exclude, rownames(seur@assays[[assay_name]]))]
    seur <- subset(seur, features = features_to_keep)
  }
  
  if(!is.null(dataset)){
    cells_to_keep <- colnames(seur)[grepl(paste0("^",dataset, "$"), seur[["dataset"]][,1])]
    seur <- subset(seur, cells = cells_to_keep)
  }
  
  if(!is.null(filter_column)){
    cells_to_keep <- colnames(seur)[grepl(paste0("^",filter_value, "$"), seur[[filter_column]][,1])]
    seur <- subset(seur, cells = cells_to_keep)
  }
  
  y <- Seurat2PB(seur, sample = "hash.mcl.ID", cluster = trt_group)
  keep.samples <- y$samples$lib.size > lib.size
  y <- y[, keep.samples]
  meta <- seur[[]] %>% distinct(hash.mcl.ID, .keep_all = T) 
  meta <- meta[match(y$samples$sample, meta$hash.mcl.ID),]
  y$samples$group <- NULL
  y$samples <- bind_cols(meta, y$samples)
  
  keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count = min.count, min.total.count = min.total.count)
  table(keep.genes)
  y <- y[keep.genes, , keep=FALSE]
  y <- normLibSizes(y)
  
  return(y)
  
}

#' wrapper to run scDist
#' @title
#' @param obj seurat object
#' @param fixed.effects character column with groups of interest (maximum two levels)
#' @param random.effects character vector of columns with batch information
#' @param assay seurat assay to use to calculate distance
#' @param ident_column character name of column of clusters
#' @param subset_cells regex format for clusters to keep
#' @param d numeric number of PC dimensions used to calculate distance
#' @param nfeats numeric number of genes used to calculate PCA
#' @return
#' @author dylanmr
#' @export
run_scdist <- function(obj, fixed.effects = "geno", assay="SCT", 
                            ident_column = "predicted.celltype", 
                            random.effects = c("hash.mcl.ID", "orig.ident","sex_predicted","time", "treatment"),
                            d=20,
                            nfeats=5000) {

  DefaultAssay(obj) <- "RNA"
  obj <- process_seurat(obj, method = "qpoisson", cluster=F, nfeats = nfeats)
  mat <- GetAssayData(obj, slot="scale.data", assay=assay) %>%  as.matrix()
  scd_res <- scDist(mat, meta.data=obj@meta.data, min.counts.per.cell = 5,
         fixed.effects = fixed.effects, 
         random.effects = random.effects, 
         clusters = ident_column, d = d)
  return(scd_res)

}

#' Calculate DEG using edgeR
#'
#' @title
#' @param y pseudobulk obj generated from edger_prep
#' @param mm formula for differential expression
#' @param contrast vector to construct comparisons
#' @return fitted model
#' @author dylanmr
#' @export
run_edger <- function(y, mm = "~ 0 + cluster + sex_predicted + seq_pool") {

  design <- model.matrix(as.formula(mm), data = y$samples)
  colnames(design) <- gsub("/","_",colnames(design))
  y <- estimateGLMRobustDisp(y, design)
  fit <- glmFit(y, design)
  return(fit)

}


#' Create Milo object
#' @title
#' @param obj seurat object
#' @param pca_reduction_name name of where PCA is stored in seurat object
#' @param umap_reduction_name name of where UMAP is stored in seurat object
#' @param k numeric number of neighbors
#' @param d numeric number of dimensions
#' @param prop 
#' @param refinement_scheme
#' @return
#' @author dylanmr
#' @export

seurat_to_milo <- function(
  seurat_obj,
  pca_reduction_name = "pca",
  umap_reduction_name = "umap",
  k = 30,
  d = 30,
  prop = 0.1,
  refinement_scheme = "graph",
  sample_column = "orig.ident",
  additional_id_column = "hash.mcl.ID"
) {
  # Convert Seurat object to SingleCellExperiment object
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Check if specified reductions exist and add them
  if (pca_reduction_name %in% names(seurat_obj@reductions) && 
      umap_reduction_name %in% names(seurat_obj@reductions)) {
    reducedDim(sce, "PCA", withDimnames = TRUE) <- 
      seurat_obj[[pca_reduction_name]]@cell.embeddings
    reducedDim(sce, "UMAP", withDimnames = TRUE) <- 
      seurat_obj[[umap_reduction_name]]@cell.embeddings
  } else {
    stop("Specified reduction names not found in Seurat object")
  }
  
  # Create Milo object
  milo_obj <- Milo(sce)
  milo_obj <- buildGraph(milo_obj, k = k, d = d)
  milo_obj <- makeNhoods(milo_obj, prop = prop, k = k, d = d, 
                         refined = TRUE, refinement_scheme = refinement_scheme)
  
  # Add metadata and count cells
  if (!is.null(additional_id_column)) {
    colData(milo_obj)$ObsID <- paste(colData(milo_obj)[[sample_column]], 
                                     colData(milo_obj)[[additional_id_column]], 
                                     sep = "_")
  } else {
    colData(milo_obj)$ObsID <- colData(milo_obj)[[sample_column]]
  }
  
  milo_obj <- countCells(milo_obj, meta.data = data.frame(colData(milo_obj)), 
                         samples = "ObsID")
  
  return(milo_obj)
}
#' Create dds object
#' @title
#' @param pb pseudobulk object
#' @param design formula to construct design matrix
#' @author dylanmr
#' @export

create_dds <- function(pb, design) {
    rownames(pb$samples) <- colnames(pb$counts)
    dds <- DESeqDataSetFromMatrix(countData = getCounts(pb), colData = pb$samples, 
                            design = as.formula(design))
    return(dds)
}

create_vsd <- function(dds, design, batch1, batch2) {
    vsd <- vst(dds, blind=F)
    if(!is.null(batch1)) {
       assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch =  colData(vsd)[[batch1]], design = model.matrix(as.formula(design), data = colData(vsd)))
    } else if(!is.null(batch2)) {
      assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch = colData(vsd)[[batch1]], batch2= colData(vsd)[[batch2]], design = model.matrix(as.formula(design), data = colData(vsd)))
    }
    return(vsd)
}()


#'test reduced model
#' @title
#' @param dds DESeq2 object
#' @param reduced reduced formula to test against full design
#' @return
#' @author dylanmr
#' @export

lrt_dds <- function(dds, reduced =  ~ treatment + time + hash_pool) {
    DESeq(dds, test="LRT", reduced = reduced)
}

lrt_edger <- function(pb_lep, term_to_match) {
    coefs <- grep(term_to_match, colnames(pb_lep))
    res <- topTags(glmLRT(pb_lep, coef=coefs), n=Inf)  %>%  data.frame()
    return(res)
}

#' create training and test sets for splsda
#' @title
#' @param obj seurat object
#' @param folds number of random folds to create from the dataset
#' @param repeats how many times to create the number of folds
#' @param group grouping of cells
#' @return
#' @author dylanmr
#' @export

create_model_data <- function(obj, folds = 10, repeats = 2, group = "group", nfeats = 10000, method = "qpoisson") {

    out <- vfold_cv(data.frame(Cells(obj)), v = folds, repeats = repeats)
    training_sets <- purrr::map(out$splits, function(x) {
      train <- subset(obj, cells = colnames(obj)[x$in_id])
      min_cells <- min(table(train[[group]][,1]))
      DefaultAssay(train) <- "RNA"
      train <- process_seurat(train, method=method, cluster=F, nfeats = nfeats, features = "^Gm|Rik$")
      equalnums <- as.character(unlist(tapply(colnames(train), train[[group]][,1], function(x) sample(x, min_cells))))
      train <- subset(train, cells = equalnums)
      return(train)
    })
    
    testing_sets <- purrr::map(out$splits, function(x) {
      test <- subset(obj, cells = colnames(obj)[-x$in_id])
      DefaultAssay(test) <- "RNA"
      test <- process_seurat(test, method=method, cluster=F, nfeats = nfeats)
      return(test)
    })

    return(list("training_data" = training_sets, "testing_data" = testing_sets))
}

#' train splsda models on cross-fold data
#' @title
#' @param cv_data crossfold data
#' @param ext_datasets other datasets to ensure overlapping gene names 
#' @param ncomp number of components to retain in model
#' @param keepX number of genes to use in spls-DA mdodel
#' @return
#' @author dylanmr
#' @export
#' 
train_models <- function(cv_data, ext_datasets = camp, ncomp = 2, keepX = 200) {
  purrr::map2(cv_data$training_data,cv_data$testing_data, function(x,y) {
    features <- intersect(rownames(x@assays$SCT@data), rownames(y@assays$SCT@data))
    features <- intersect(features, rownames(ext_datasets@assays$SCT@data))
    train_x <- Matrix::t(x@assays$SCT@data[features,])
    train_y <- factor(x$treatment)
    final.splsda <- splsda(X = train_x, Y = train_y,  ncomp = ncomp, keepX = keepX)
    return(final.splsda)
  })
}


LoadXeniumNew <- function(data.dir, fov = 'fov', assay = 'Xenium') {
  data <- ReadXenium(
    data.dir = data.dir, 
    type = c("centroids", "segmentations"),
  )

  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )

  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  if("Blank Codeword" %in% names(data$matrix))
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  else
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])

  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}



project_pca <- function(ref_data, new_data, celltype, n_pcs = 10) {
  
  # Subset the reference and new data based on provided cluster and cell type
  ref_subset <- subset(ref_data, subset = seurat_clusters == celltype)
  new_subset <- subset(new_data, subset = predicted.celltype == celltype)
  
  # Find shared features
  shared_feats <- intersect(rownames(ref_subset), rownames(new_subset))
  
  # Process the reference data
  ref_subset <- process_seurat(subset(ref_subset, features = shared_feats), method = "integrate", batch = "hash_pool", cluster = FALSE, nfeats = 5000)
  
  # Perform PCA on the reference data
  ref_mat <- ref_subset@assays$integrated@scale.data
  ref_pca <- prcomp(t(ref_mat), center = TRUE, scale. = FALSE)

  
  # Normalize and scale the new data using the reference's variable features
  new_subset <- NormalizeData(new_subset) %>% ScaleData(., features = VariableFeatures(ref_subset))
  scale_mat <- new_subset@assays$RNA@scale.data
  
  # Align the rows of the new data's scaled matrix with the reference matrix
  scale_mat <- scale_mat[match(rownames(ref_mat), rownames(scale_mat)),]
  
  # Project the new data into the PCA space of the reference data
  new_data_pca_scores <- as.matrix(t(scale_mat)) %*% ref_pca$rotation  %>% data.frame() 
  
  return(list(ref_pca = data.frame(ref_pca$x)  %>%  bind_cols(ref_subset[[]]), new_data_pca_scores = new_data_pca_scores  %>%  bind_cols(new_subset[[]])))
}
