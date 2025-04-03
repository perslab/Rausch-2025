#' Create Seurat object from 10X h5 output files and 
#' cellbender h5 output files
#' @title generate seurat object
#' @return seurat object
#' @import Seurat scCustomize
#' @author dylanmr
#' @export

gen_seurat <- function(denoisedfile, rawpath, pattern = "_out.h5", min.cells=10, min.features=1000, cellbender = T) {
  
  name <- gsub(pattern,"",gsub(".*/","", denoisedfile))
  raw_file <- paste0(rawpath, "/", name, "_raw_feature_bc_matrix.h5")
  rawmat <- Read10X_h5(raw_file)
  colnames(rawmat) <- gsub("-1", replacement = "", colnames(rawmat))
  
  if(isTRUE(cellbender)){
    cellbender_mat <- Read_CellBender_h5_Mat(file_name = denoisedfile)
    colnames(cellbender_mat) <- gsub("-1", replacement = "", colnames(cellbender_mat))
    seur <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cellbender_mat, 
                                                   raw_counts_matrix = rawmat,
                                                   raw_assay_name = "RAW", min_cells = min.cells, 
                                                   min_features = min.features, project = as.character(name))
    seur <- Add_CellBender_Diff(seurat_object = seur, raw_assay_name = "RAW", cell_bender_assay_name = "RNA")
  } else {
    seur <- CreateSeuratObject(rawmat, min.cells = min.cells, min.features = min.feature, project = as.character(name))
  }

  return(seur)
  
}

#' Add HTO matrix to seurat object
#' Specifically reads in format from kb count (https://www.kallistobus.tools/kb_usage/kb_count/)
#'
#' @title Add HTO
#' @importFrom BUSpaRse read_count_output
#' @return seurat object with HTO matrix
#' @author dylanmr
#' @export
add_hto <- function(seur, path = "sc_preprocessing/hto_counts/") {
  
  hpool <- unique(seur$orig.ident)
  print(hpool)
  lane <- paste0(here::here(path), "/", hpool, "/", "counts_unfiltered")
  # Read HTO mat
  hto <- BUSpaRse::read_count_output(
    dir = lane, 
    name = "cells_x_features", 
    tcc = FALSE
  )
  
  bc_overlap <- intersect(colnames(hto), colnames(seur))
  seur <- subset(seur, cells = bc_overlap)
  seur[["HTO"]] <- CreateAssayObject(counts = hto[,bc_overlap])
  return(seur)
  
}

#' Functions for assigning cells to hash tags
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @import assertthat mclust diptest
#' @author dylanmr
#' @export

demux_hto <- function(obj, outs = "sc_preprocessing/hto_qc/", assay = "HTO") {
  outs <- here::here(outs)
  lane <- as.character(unique(obj$orig.ident))
  obj <- NormalizeData(obj, assay = assay, normalization.method = 'CLR', margin = 2, verbose = F)
  # Setting the threshold based on the quantiles of the negative and the positive clusters.
  q_l = 1
  q_h = 0.001
  obj <- HTODemux.mcl(obj, q_l =  q_l, q_h = q_h, assay=assay)
  table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Doublet']
  cat('#\t-\t\tDoublets:\t\t',table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Doublet'],'\n',
      '#\t-\t\tNegatives:\t\t',table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Negative'],'\n',
      '#\t-\t\tSinglets:\t\t',table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Singlet'],'\n',
      sep = '', file = paste0(outs, lane,".log"))
  HTODemux_mcl.visualization(obj, q_l =  q_l, q_h = q_h, outs = outs, lane=lane)
  return(obj)
  
}

is_multimodal <- function(x, p.cutoff = 1e-2) {
  # Test if the expression distribution violates unimodal distribution.
  p = diptest::dip.test(x)$p.value
  return(p < p.cutoff)
}

select_hash_cutoff_mcl <- function(x, q_l = 1, q_h = 0.01, seed = 42) {
  # Model HTO data as a mixture of two Gaussian distributions (for normalized [across cells] data)
  # And select HTO cutoff based on mclust (Model based clustering).
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  # q_l and q_h are the quantiles for negative and postive cluster, respectively. 
  cutoff <- max(quantile(x[cl == low_cl], q_l), quantile(x[cl == high_cl], q_h))
  # The higher the cut off, the less false positive (the more false negative).
  return(cutoff)
}

hash_mcl_p <- function(x, seed = 3030, q_l = 1, q_h = 0.001) {
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  p.high_cl <- km$z[,high_cl]
  # Correct assignment error from Mclust
  p.high_cl[which(x < max(quantile(cl_center[low_cl], q_l), quantile(cl_center[low_cl], q_h)))] = 0
  names(p.high_cl) = names(x)
  return(p.high_cl)
}

HTO_classifcation = function(discrete, hto_mcl.p, assay){
  # Based on HTODemux (Seurat)
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = discrete)
  hash.max <- apply(X = hto_mcl.p, MARGIN = 2, FUN = max) # This returns the probability of the most likely HashID (based on the Hashtag distribution among cells)
  hash.maxID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = which.max)])
  hash.second <- apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) sort(x,decreasing = T)[2])
  hash.secondID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) order(x,decreasing = T)[2])])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, 'mcl', c("maxID", "secondID", "margin", "classification", "classification.global"), sep = "_")
  return(classification.metadata)
}

HTODemux.mcl <- function(object, assay = "HTO", q_l = 1, q_h = 0.005, seed = 42){
  # A function to find the threshold for each hastag, the P, and singlet, doublets and negative.
  # The input is the HTO data matrix (normalized across cells).
  
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x, q_l = q_l, q_h = q_h), future.seed = T))
  hto_mcl.cutoff$Multi_modal = apply(data,1,function(x) is_multimodal(x))
  print(hto_mcl.cutoff)
  hto_mcl.p = t(apply(data,1,function(x) hash_mcl_p(x, seed = seed, q_l = q_l, q_h = q_h)))
  discrete <- data
  discrete[discrete > 0] <- 0
  for (iter in rownames(x = data)) {
    values <- data[iter, ]
    cutoff <- hto_mcl.cutoff[iter,'cut_off']
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
  }
  classification.metadata <- HTO_classifcation(discrete, hto_mcl.p, assay)
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste(assay, "mcl", "classification", sep = '_')
  doublets <- rownames(x = object[[]])[which(object[[paste(assay, "mcl","classification.global", sep = "_")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  Idents(object) = factor(Idents(object), levels = c('Doublet', 'Negative', rownames(object@assays$HTO)))
  object$hash.mcl.ID <- Idents(object = object)
  return(object)
}

HTODemux_mcl.visualization <- function(object, assay = "HTO", q_l = 1, q_h = 0.005, seed = 42, outs, lane){
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto.data.wide = data.frame(t(data))
  hto.data.long = data.table::melt(data.table::setDT(hto.data.wide,keep.rownames = T), id.vars = 'rn', value.name = 'Expression',variable.name = 'hto')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x, q_l = q_l, q_h = q_h), future.seed = T), hto = colnames(hto.data.wide)[-1])
  
  p <- ggplot(hto.data.long, aes(x = Expression)) +
    geom_histogram(bins = 100) +
    geom_vline(data = hto_mcl.cutoff, aes(xintercept = cut_off), col = 'red') +
    facet_wrap(~hto,scales = 'free',ncol = 2) +
    xlab('Expression') +
    ylab('Counts') +
    scale_y_sqrt() +
    ggtitle('Individual HTO distributions', subtitle = paste0('Q_negative = ', q_l, '; Q_positive = ', q_h)) +
    theme_minimal(base_size = 10)
  
  ggsave(filename = paste0(lane,'_HTO-hist.png'),
         plot = p,
         path = outs,
         width = 10,
         height = 12)
}

#' Add excel file to existing seurat object
#' @title add_metadata
#' @param joinby
#' @return seurat object
#' @author dylanmr
#' @export

add_metadata <- function(obj, file = here::here("sc_preprocessing/sequenced_mice_data_swapped_w_all.xlsx"), joinby = "hash") {
  
  data <- readxl::read_xlsx(file)  %>%  janitor::clean_names() 
  add_data <- left_join(obj[["hash.mcl.ID"]], data, by = c("hash.mcl.ID" = {{ joinby }})) 
  row.names(add_data) <- row.names(obj[[]])
  obj <- AddMetaData(obj, metadata = add_data)
  obj <- subset(obj, cells = colnames(obj)[which(!is.na(obj$dataset))])
  return(obj)
  
}

#' Assign cell annotations
#' @title CellAnnotator
#' @return seurat object with updated metadata column
#' @author dylanmr
#' @import CellAnnotatoR
#' @export

classify_cells <- function(obj, column_name = "annotation",
                           path = here::here("cell_markers/coarse_markers.txt"), 
                           clustering_res = clustering_res,
                           clustering_dims = clustering_dims) {
  
  obj <- process_seurat(obj, method = "log", res = clustering_res, dims = clustering_dims) 
  cm <- obj@assays$RNA@counts
  clusters <- setNames(obj@meta.data$seurat_clusters, rownames(obj@meta.data))
  clf_data <- getClassificationData(cm, markers = path)
  obj[[column_name]] <- assignCellsByScores(graph = NULL, clf_data, clusters=clusters)$annotation$l1
  return(obj)
  
}

#' process seurat objects
#'
#' @title
#' @param obj seurat or list of seurat objects
#' @param method log, glm, qpoisson, or integrate
#' @param res numeric value indicating resolution of clustering
#' @param features regex for features to remove
#' @param dims number of dimensions to use for umap
#' @param batch character matching to metadata to integrate on. if obj is a list will 
#' integrate by list entry
#' @param return_model logical, return umap model. useful for projection
#' @param cluster logical, perform RunUMAP, FindClusters, FindNeighbors
#' @param type character, 
#' @return
#' @author dylanmr
#' @export

process_seurat <- function(obj, method, ref_datasets = NULL, k.anchor=5,k.weight=100,
                           res=NULL, features=NULL, dims=NULL, 
                           batch=NULL, return_model = FALSE, cluster=TRUE, 
                           nfeats = 2000) {
  
  if (is(obj, "SingleCellExperiment")) {
    obj <- CreateSeuratObject(counts = counts(obj))
  }
  
  if(!is.null(features)) {
    features_to_remove <- rownames(obj)[grepl(features, rownames(obj))]
    obj <- subset(obj, features = setdiff(rownames(obj), features_to_remove))
  }
  
  obj <- switch(
    method,
    integrate = process_integrate(obj, batch, nfeats, ref_datasets, k.anchor, k.weight),
    log = process_log(obj, nfeats, batch),
    glm = process_glm(obj, nfeats, batch),
    qpoisson = process_qpoisson(obj, nfeats, batch),
    stop(paste("Unknown method:", method))
  )
  
  if (cluster) {
    obj <- cluster_obj(obj, dims, return_model, res)
  }
  
  return(obj)
}

process_integrate <- function(obj, batch, nfeats, ref_datasets, k.anchor, k.weight) {
  if(is.list(obj)) {
    obj_list <- obj
  } else {
    DefaultAssay(obj) <- "RNA"
    obj_list <- SplitObject(obj, split.by = batch)
  }
  
  obj_list <- lapply(X = obj_list, FUN = function(x) {
    if(sum("SCT" %in% names(x@assays))>0) {
      x[["SCT"]] <- NULL
    }
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeats)
  })
  
  features <- SelectIntegrationFeatures(object.list = obj_list)
  
  obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  anchors <- FindIntegrationAnchors(object.list = obj_list, reference = ref_datasets,
    k.anchor = k.anchor, anchor.features = features, reduction = "rpca")
  integrated <- IntegrateData(anchorset = anchors, k.weight=k.weight)
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated) %>% RunPCA(.)
  return(integrated)
  
}

process_log <- function(obj, nfeats, batch) {
  if (inherits(obj, "list")) {
    obj <- obj[[1]]
  }
  
  DefaultAssay(obj) <- "RNA"
  
  NormalizeData(obj) %>%
    FindVariableFeatures(., selection.method = "vst", nfeatures = nfeats) %>%
    ScaleData(., vars.to.regress=batch) %>%
    RunPCA(.)
}

process_glm <- function(obj, nfeats, batch) {
  if (inherits(obj, "list")) {
    obj <- obj[[1]]
  }
  
  DefaultAssay(obj) <- "RNA"

  SCTransform(obj, method = "glmGamPoi", batch_var=batch, variable.features.n = nfeats) %>%
    RunPCA(.)
}

process_qpoisson <- function(obj, nfeats, batch) {
  if (inherits(obj, "list")) {
    obj <- obj[[1]]
  }
  
  DefaultAssay(obj) <- "RNA"
  
  SCTransform(obj, method = "qpoisson", variable.features.n = nfeats, vars.to.regress = batch) %>%
    RunPCA(.)
}

cluster_obj <- function(obj, dims, return_model, res) {
  obj %>%
    RunUMAP(., dims = seq(dims), return.model=return_model) %>%
    FindNeighbors(., dims = seq(dims)) %>%
    FindClusters(., resolution = res)
}


#' Run gene expression qc to remove outlier cells
#'
#' @title
#' @param obj seurat object
#' @param times numeric, number of tests that must fail for a cell to be removed
#' @param batch character, indicates metadata column to split calculation of outlier on
#' @return
#' @author dylanmr
#' @export
run_gex_qc <- function(obj, times, batch) {
  
  # convert to single cell experiment
  ds <- as.SingleCellExperiment(obj)
  ds <- scuttle::addPerCellQC(ds, percent_top=c(20,50,100,200))
  
  # calculate qc columns
  ds$total_features <- ds$detected
  ds$log10_total_features <- log10(ds$detected)
  ds$total_counts <- ds$sum
  ds$log10_total_counts <- log10(ds$sum+1)
  ds$featcount_ratio <- ds$log10_total_counts/ds$log10_total_features
  ds$featcount_dist <- getFeatCountDist(ds)
  ds$pct_counts_top_20_features <- colData(ds)[[intersect(c("percent_top_20","pct_counts_in_top_20_features","percent.top_20"), colnames(colData(ds)))[[1]]]]
  ds$pct_counts_top_50_features <- colData(ds)[[intersect(c("percent_top_50","pct_counts_in_top_50_features","percent.top_50"), colnames(colData(ds)))[[1]]]]
  
  # build matrix to track qc
  discard <- matrix(c(scuttle::isOutlier(ds$log10_total_counts, nmads = 5, type = "lower", batch = ds[[batch]]) , 
                      scuttle::isOutlier(ds$log10_total_counts,nmads = 2.5,type = "higher", batch = ds[[batch]]) ,
                      scuttle::isOutlier(ds$log10_total_features,nmads = 5,type = "lower", batch = ds[[batch]]) , 
                      scuttle::isOutlier(ds$log10_total_features, nmads = 2.5,type = "higher", batch = ds[[batch]]) ,
                      scuttle::isOutlier(ds$pct_counts_top_20_features, nmads = 5,type = "both", batch = ds[[batch]]) ,
                      scuttle::isOutlier(ds$pct_counts_top_50_features, nmads = 5,type = "both", batch = ds[[batch]]),
                      scuttle::isOutlier(ds$featcount_dist, nmads = 5,type = "both", batch = ds[[batch]]), 
                      ds$HTO_mcl_margin < (mean(ds$HTO_mcl_margin) - sd(ds$HTO_mcl_margin)*3)),
                    nrow = length(ds$log10_total_counts))

  # filter cells
  discard <- rowSums(discard)>=times
  cells_to_keep <- colnames(ds)[!discard]
  obj <- subset(obj, cells = cells_to_keep)

}

getFeatCountDist <- function(df, do.plot=FALSE, linear=TRUE){
  if(is(df,"SingleCellExperiment")) df <- as.data.frame(colData(df))
  if(linear){
    mod <- lm(log10_total_features~log10_total_counts, data=df)
    pred <- predict.lm(mod, newdata=data.frame(log10_total_counts=df$log10_total_counts))
  }else{
    mod <- loess(log10_total_features~log10_total_counts, data=df)
    pred <- predict.loess(mod, newdata=data.frame(log10_total_counts=df$log10_total_counts))
  }
  df$diff <- df$log10_total_features - pred
  if(do.plot){
    library(ggplot2)
    ggplot(df, aes(x=total_counts, y=total_features, colour=diff)) + 
      geom_point() + geom_smooth(method = "loess", col="black")
  }
  df$diff
}

#' Remove doublets from seurat object
#'
#' @title
#' @param obj seurat object
#' @param k number of neighbors to search across for putative doublets
#' @return seurat object
#' @author dylanmr
#' @export

remove_doublets <- function(obj, k=30) {
  
  # Remove negative cells
  obj <- subset(obj, subset = HTO_mcl_classification.global == "Negative", invert=T)
  
  # Collect Data for scDblFinder
  dat <- GetAssayData(obj, slot = "data")
  dat <- dat[rownames(dat) %in% obj@assays$RNA@var.features,]
  doublets <- as.logical(ifelse(obj$HTO_mcl_classification.global == "Doublet", T, F))
  res <- scDblFinder::recoverDoublets(dat, doublets=doublets, samples=table(obj$hash.mcl.ID), k=k)
  print(paste(sum(res$predicted),"intra-sample doublets detected"))
  
  # Remove all doublets
  obj <- subset(obj, cells = Cells(obj)[which(res$predicted==T)], invert=T)
  obj <- subset(obj, subset= HTO_mcl_classification.global == "Singlet")
  
  return(obj)
  
}

#' generate xenium seurat objects
#' @title
#' @param path path to xenium files
#' @param min_umi minimum number of probes in a cell
#' @return
#' @author dylanmr
#' @export
gen_xenium <- function(path, min_umi = 10) {

  obj <- LoadXenium(path, fov="fov")
  pattern <- "__(R.*?)__"
  obj[["section"]] <- str_match_all(path, pattern)[[1]][,2]
  pattern <- "__(00.*?)__"
  obj[["slide"]] <- str_match_all(path, pattern)[[1]][,2]
  obj <- subset(obj, subset = nCount_Xenium > min_umi)
  return(obj)

}

#' Map campbell reference data
#' 
#' hardcoded to combine both agrp populations into a single annotation
#' 
#' @title
#' @param obj seurat object
#' @param dims number of dimensions for reference mapping
#' @param column_name metadata column to contain mapped information
#' @param ref campbell seurat object
#' @param dims
#' @return
#' @author dylanmr
#' @export

map_ref <- function(obj = obj, dims = 30,
                    ref = ref, ref_labels = "stablecamp7",
                    column_name = "final_labels") {
  
  anchors <- FindTransferAnchors(reference = ref, query = obj, dims = seq(dims), reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = ref[[ref_labels]][,1], dims = seq(dims))
  obj[[column_name]] <- predictions$predicted.id
  obj[["prediction.score.max"]] <- predictions$prediction.score.max
  return(obj)
  
}

#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param nameme1
#' @return
#' @author dylanmr
#' @export

prepCamp <- function(method) {
  
  camp <- scRNAseq::CampbellBrainData()
  seur <- process_seurat(camp, method = method, cluster=F, type="sce")
  seur <- AddMetaData(seur, metadata = data.frame(colData(camp)))
  return(seur)
  
}

#' prepare martin myers sc data

#' @title
#' @param path path to data
#' @param min.features filter cells based on the minimum number of features
#' @param min.cells filter genes based on the minimum number of cells in which they are identified
#' @param name identity to include in seurat object
#' @returns returns a processed seurat object
#' @author dylanmr

prep_myers <- function(path, min.features, min.cells, name = "leprgfp") {
  
  obj <- Read10X(here::here(path), gene.column = 1, strip.suffix = T)
  colnames(obj) <- stri_sub(stri_sub(colnames(obj),to =  -2), 2)
  meta <- read.csv(here::here(path, "metadata.csv.gz"), row.names = "ID")
  obj <- CreateSeuratObject(obj, min.features = min.features, min.cells = min.cells, project = name, meta.data = meta)
  obj <- subset(obj, subset = Cluster_region == "HY")
  obj <- subset(obj, subset = Lepr_cluster %in% c("L10.GABA","L8.Nts"), invert=T)
  clustokeep <- obj$Lepr_cluster[grep("^L", obj$Lepr_cluster)] %>% unique()
  obj <- subset(obj, subset = Lepr_cluster %in% clustokeep)
  obj <- process_seurat(obj, method="log", res=0.8, dims=30)
  return(obj)
  
}

#' prepare in house labeling for annotation

#' @title
#' @param path path to data
#' @param min.features filter cells based on the minimum number of features
#' @param min.cells filter genes based on the minimum number of cells in which they are identified
#' @param name identity to include in seurat object
#' @returns returns a processed seurat object
#' @author dylanmr

prep_internal_lab_transfer <- function(path = "external_data/20231027_ARC_filtered.qs", ncells = 20000) {
  
  ref <- qs::qread(here::here(path))
  ref_sub <- subset(ref, cells = sample(Cells(ref), ncells))
  ref_sub$stablecamp7 <- ifelse(grepl("Agrp",ref_sub$stablecamp7), "Agrp", ref_sub$stablecamp7)
  return(ref_sub)
  
}

prep_macaque <- function(path = "external_data") {
  genes <- read.table(here::here(path, "myers_macaque/features.tsv.gz"), sep="\t", header=T)
  bars <- read.table(here::here(path, "myers_macaque/barcodes.tsv.gz"), sep="\t", header=T)
  mtx <- Matrix::readMM(here::here(path,"myers_macaque/matrix.mtx.gz"))
  rownames(mtx) <- genes[,1]
  colnames(mtx) <- bars[,1]
  
  orth_genes <- read.table(here::here(path, "gene_orth.txt"), sep="\t", header=T) %>% janitor::clean_names() %>% mutate_all(na_if,"")
  orth_genes_filt <- orth_genes[!is.na(orth_genes$mouse_gene_name),]
  orth_genes_filt <- orth_genes_filt[!duplicated(orth_genes_filt$mouse_gene_name),]
  orth_genes_filt <- orth_genes_filt[!duplicated(orth_genes_filt$macaque_gene_stable_id),]
  overlap <- intersect(rownames(mtx), orth_genes_filt$macaque_gene_stable_id)
  orth_genes_final <- orth_genes_filt[orth_genes_filt$macaque_gene_stable_id %in% overlap,]
  mtx <- mtx[intersect(rownames(mtx), orth_genes_final$macaque_gene_stable_id),]

  so <- CreateSeuratObject(mtx, min.cells=10, min.features=1000)
  so[["sample"]] <-  gsub(".*_","", colnames(so))
}

corral_reduction <- function(sobj, norm_type="sqrt") {
  if(norm_type=="sqrt"){
    int <- rownames(sobj)[rowSums(sobj@assays$Xenium@counts)>7000]
    int <- int[!int=="Sparc"]
    reduc <- corral::corral(sqrt(sobj@assays$Xenium@counts[int,]))
    sobj[["corral_sqrt"]] <- CreateDimReducObject(embeddings = reduc$v, key="corral", assay="Xenium")
  } else if(norm_type=="ft") {
    int <- rownames(sobj)[rowSums(sobj@assays$Xenium@counts)>7000]
    int <- int[!int=="Sparc"]
    reduc <- corral::corral(sqrt(sobj@assays$Xenium@counts[int,]) + sqrt(sobj@assays$Xenium@counts[int,]+1))
    sobj[["corral_ft"]] <- CreateDimReducObject(embeddings = reduc$v, key="corral", assay="Xenium")
  }
  return(sobj)
}

annotate_spatial_data <- function(sobj) {
  sobj[["sample"]] <- gsub(".*-","",colnames(sobj))
  sobj[["slide"]] <- ifelse(sobj[["sample"]][,1] %in% c("1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8"),"slide1","slide2")
  sobj[["diet"]] <- ifelse(sobj[["sample"]][,1] %in% c("1_1","1_2","1_3","1_4","1_9","1_10","1_11","1_12"),"chow","dio")
  sobj[["animal"]] <- case_when(sobj[["sample"]][,1] %in% c("1_1","1_2","1_3") ~ "animal1",
                               sobj[["sample"]][,1] %in% c("1_4","1_12") ~ "animal2",
                               sobj[["sample"]][,1] %in% c("1_5","1_6","1_7") ~ "animal3",
                               sobj[["sample"]][,1] %in% c("1_8","1_16") ~ "animal4",
                               sobj[["sample"]][,1] %in% c("1_9","1_10","1_11") ~ "animal5",
                               sobj[["sample"]][,1] %in% c("1_13","1_14","1_15") ~ "animal6")
  return(sobj)
}

convert_mouse_to_human <- function(gene_list) { 
     output = c()
     mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

     for(gene in gene_list) {
          class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[["DB.Class.Key"]]
          if(!identical(class_key, integer(0)) ) {
               human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name == "human"))[,"Symbol"]
               for (human_gene in human_genes) {
                    output <- rbind(c(gene, human_gene), output)
               }
          }
     }
     return (output)
}
