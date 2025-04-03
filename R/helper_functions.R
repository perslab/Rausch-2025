#' Convert seurat to h5 object
#' @title
#' @import zellkonverter
#' @return 
#' @author dylanmr
#' @export
convert_to_h5 <- function(obj, filename) {
  
  DefaultAssay(obj) <- "RNA"
  sce <- as.SingleCellExperiment(obj)
  colData(sce) <- colData(sce)[colSums(is.na(colData(sce))) == 0]
  writeH5AD(sce, file = filename)
  
}

#' Convert seurat to loom file
#' @return
#' @import SeuratDisk
#' @author dylanmr
#' @export
create_loom <- function(obj, filepath, subset_value, subset_column) {
  
  if(!is.null(subset_value)) {
    obj <- filter_cells_by_column(obj, subset_column, subset_value)
  }
    
  DefaultAssay(obj) <- "RNA"
  loom_ma <- SeuratDisk::as.loom(obj, filename = filepath, verbose = FALSE,  overwrite=TRUE)
  loom_ma$close_all()
  
}

#' code to merge matrices
#'
#' .. content for \details{} ..
#'
#' @title
#' @param filtered_seurat
#' @return
#' @author dylanmr
#' @export
merge_mats <- function(filtered_seurat) {
  
  mat_list <- purrr::map(filtered_seurat, function(x) {
    cts <- x@assays$RNA@counts
    colnames(cts) <- paste0(colnames(cts), "_", x$orig.ident)
    return(cts)
  })  
  
  meta_df <- purrr::map_dfr(filtered_seurat, function(x) {
    meta <- x@meta.data
    rownames(meta) <- paste0(rownames(meta), "_", x$orig.ident)
    return(meta)
  })
  
  merged_mat <- RowMergeSparseMatrices(mat_list[[1]], mat_list[-1])
  seur <- CreateSeuratObject(merged_mat, meta.data = meta_df)
  
  return(seur)
  
}

#' Helper to create new column from existing columns in seurat object
#'
#' @title
#' @return
#' @author dylanmr
#' @export
create_new_column_seurat <- function(obj, selected_columns = c("treatment", "time", "geno"), new_column_name = "group") {
  
  df <- obj[[]]
  
  # Check if all selected columns are in the dataframe
  if(!all(selected_columns %in% names(df))) {
    stop("Not all selected columns are in the dataframe")
  }
  
  # Create the new column by pasting the selected columns together
  # The 'apply' function is used to iterate over rows
  obj[[new_column_name]] <- apply(df[, selected_columns, drop = FALSE], 1, paste, collapse = "_")
  
  # Return the modified dataframe
  return(obj)
  
}

#' Predict the sex of an animal by hash-tag
#'
#' @title
#' @param x seurat object
#' @param assay character, indicate assay from object to use to calculate expression
#' @param features character, vector of genes to calculate average expression
#' @return
#' @author dylanmr
#' @export
RunSexPrediction <- function(x, features = c("Xist","Tsix"), assay= "SCT") {

  sex_genes <- AverageExpression(x, features = features, assays = assay, group.by = "hash.mcl.ID")$RNA
  x@meta.data <- full_join(x@meta.data, enframe(factor(Mclust(colMeans(sex_genes),G=2)$classification), 
                                                name = "hash.mcl.ID", value="sex_predicted"), by=c("hash.mcl.ID"))
  rownames(x@meta.data) <- colnames(x)
  return(x)

}

#' Wrapper to subset cells 
#' easier to pass to a function
#' allows use of regex

#' @title
#' @param input seurat object
#' @param column_name column to filter
#' @param values character vector in regex format to match
#' @param invert logical if true keep match, if false keep non-match
#' @return
#' @author dylanmr
#' @export
#'

filter_cells_by_column <- function(input, column_name, values, invert = FALSE) {
  
  if (!is.vector(values) || is.null(values)) {
    stop("The 'values' parameter should be a non-null vector.")
  }
  
  if (inherits(input, "Seurat")) {
    seurat_list <- list(input)
  } else if (is.list(input) && all(sapply(input, inherits, "Seurat"))) {
    seurat_list <- input
  } else {
    stop("Input must be a Seurat object or a list of Seurat objects.")
  }
  
  filtered_objects <- lapply(seurat_list, function(s) {
    # Get cells that match any of the desired values in the specified column
    if (invert) {
      cells_to_keep <- !s@meta.data[[column_name]] %in% values
    } else {
      cells_to_keep <- s@meta.data[[column_name]] %in% values
    }
    
    # Subset the Seurat object to keep only the desired cells
    s_subset <- subset(s, cells = colnames(s)[cells_to_keep])
    
    return(s_subset)
  })
  
  # Check the length of the list and return appropriately
  if (length(filtered_objects) == 1) {
    return(filtered_objects[[1]])
  } else {
    return(filtered_objects)
  }
  
  return(filtered_objects)
  
}

#' project a query dataset into the reference embedding space of a second dataset
#' @title
#' @param ref reference dataset in the form of a seurat object. required to have run return_umap=T when processing the dataset
#' @param query query dataset in the form of a seurat object
#' @param label_to_transfer name of column in the reference object containing cluster identities
#' @param dims number of dimensions to use for label transfer
#' @param assay assay to use for integration
#' @return
#' @author dylanmr
#' @export
project_umap <- function(ref, query, label_to_transfer="predicted.id", dims=30, reference.assay = reference.assay, query.assay = query.assay) {
  
  anchors <- FindTransferAnchors(reference = ref, query = query, reference.assay = reference.assay, query.assay = query.assay, 
                                 dims = seq(dims), reference.reduction = "pca")
  
  query <- MapQuery(anchorset = anchors, reference = ref, query = query,
                    refdata = list(celltype = label_to_transfer), reference.reduction = "pca", reduction.model = "umap")
  
  return(query)
  
}

#' modify a seurat column
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
modify_column <- function(obj, column_to_update, original, update) {
  obj[[column_to_update]] <- ifelse(obj[[column_to_update]][,1] %in% original, update, obj[[column_to_update]][,1])
  return(obj)
}