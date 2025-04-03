#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
infer_tf_activity <- function(deseq_obj) {

  # extract tf database
  dorothea_df <- dorothea_mm %>%
    dplyr::filter(confidence %in% c("A", "B", "C")) %>%
    dplyr::select(target, tf, mor) %>%
    as.data.frame()
  dorothea_df$likelihood <- 1
  
  #normalize and stabilize pseudobulk
  vsd <- DESeq2::vst(deseq_obj, blind=F)
  #print(unique(vsd$final_labels) %>%  as.character)
  
  results <- 
    decouple(
      mat = assay(vsd),
      network = dorothea_df,
      .source = "tf",
      .target = "target",
      statistics = c("wmean", "gsva", "ulm"),
      args = list(
        gsva = list(verbose = FALSE),
        wmean = list(.mor = "mor", .likelihood = "likelihood"),
        ulm = list(.mor = "mor", .likelihood = "likelihood")
      ),
      minsize = 5
    )
  
  # consensus_results <- run_consensus(results) %>%  
  #   inner_join(colData(vsd) %>%  as.data.frame() %>% 
  #                rownames_to_column("condition") %>%  
  #                dplyr::select(condition, geno, hash_pool, final_labels)) %>% 
  #   nest(-source) %>% 
  #   mutate(lm = purrr::map(data, ~lm(score ~ geno, data=.x) %>% broom::tidy(.)),
  #          celltype = purrr::map_chr(data, ~unique(as.character(.x$final_labels)))) %>% 
  #   unnest(lm) %>% 
  #   arrange(p.value) %>% 
  #   dplyr::filter(term == "genowt") %>%  
  #   arrange(-statistic)
  
  results <- results %>% 
    inner_join(
      colData(vsd) %>%
        as.data.frame() %>% 
        rownames_to_column("condition") %>%
        dplyr::select(condition, geno, hash_pool, final_labels)
    )
  
  return(results)
}
