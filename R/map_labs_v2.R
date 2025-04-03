#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param obj
#' @param type
#' @param dims
#' @return
#' @author dylanmr
#' @export


map_labs_v2 <- function(obj, type = "coarse", dims = 30) {
  # Validate inputs
  if (!is(obj, "Seurat")) {
    stop("The object provided is not a Seurat object.")
  }
  if (!type %in% c("neuron", "coarse")) {
    stop("Invalid type specified. Valid types are 'neuron' or 'coarse'.")
  }
  if (!is.numeric(dims) || dims <= 0) {
    stop("Parameter 'dims' must be a positive integer.")
  }
  
  # Path configurations (could be set externally)
  paths <- list(
    neuron = here::here("data/Arc_Neurons_33_Clusters.rds"),
    coarse = here::here("data/Arcuate_Full_Atlas.rds")
  )
  
  # Load the reference data
  camp <- readRDS(paths[[type]])
  
  # Set default assays
  DefaultAssay(camp) <- if (type == "neuron") "integrated" else "RNA"
  DefaultAssay(obj) <- if (type == "neuron") "integrated" else "RNA"
  
  # Normalize and find variable features for the coarse type
  if (type == "coarse") {
    obj <- NormalizeData(obj) %>% FindVariableFeatures()
  }
  
  # Perform the label transfer
  camp[["label"]] <- Idents(camp)
  anchors <- FindTransferAnchors(reference = camp, query = obj, dims = seq_len(dims), reduction = "cca")
  predictions <- TransferData(anchorset = anchors, weight.reduction = "cca", refdata = camp$label, dims = seq_len(dims))
  
  # Add metadata and assign new labels
  obj <- AddMetaData(obj, metadata = predictions[, c(1, ncol(predictions))])
  obj[["camp_labs"]] <- ifelse(grepl("Agrp", obj$predicted.id), "Agrp", obj$predicted.id)
  
  return(obj)
}
