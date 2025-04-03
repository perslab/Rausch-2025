#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
predict_activity <- function(method, obj, ref_data = NULL, ngene) {

  if(method == "pca") {
    
    keep_genes <- intersect(rownames(assay(ref_data)), rownames(obj))
    vsd <- assay(ref_data)[keep_genes,]
    vsd <- vsd[head(order(-matrixStats::rowVars(vsd)),ngene),]
    pc <- prcomp(t(vsd))
    act_scores_pc <- predict(pc, t(obj@assays$alra@data)) %>% data.frame() %>% rownames_to_column('cell') %>% full_join(obj[[]] %>%  rownames_to_column('cell'))
    
  } else if(method == "IEG") {
    
    activity.gene.sets <- list(IEG = c('Fosb','Npas4','Fos','Junb','Nr4a1','Egr4','Arc','Egr2','Egr1','Maff','Ier2','Klf4','Dusp1','Gadd45g','Dusp5','Egr3','Btg2','Ppp1r15a','Amigo3'))
    res <- AddModuleScore(obj, features = activity.gene.sets)
    act_scores_ieg <- res[[]]
    
  }
}
