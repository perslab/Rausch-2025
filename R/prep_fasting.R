#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
prep_fasting <- function(path=NULL) {

  fast <- qs::qread(file = here::here("data/hypomap/arc_hypomap_Dowsett.qs"))
  fast <- subset(fast, subset = Batch_ID == "Dowsett10xnuc_batch_1")
  fast <- process_seurat(fast, method = "qpoisson", res=0.8, dims = 40)
  return(fast)
  
}
