#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param obwt_dds
#' @return
#' @author dylanmr
#' @export
infer_path_activity <- function(pb) {

  #normalize and stabilize pseudobulk
  vsd <- DESeq2::vst(pb, blind=F)
  #print(unique(vsd$final_labels) %>%  as.character)
  
  net <- get_progeny(organism = 'mouse', top = 100)
  
  results <- 
    decouple(
      mat = assay(vsd),
      network = net,
      .source = "source",
      .target = "target",
      statistics = c("wmean", "gsva", "ulm"),
      args = list(
        gsva = list(verbose = FALSE),
        wmean = list(.mor = "weight", .likelihood = "likelihood"),
        ulm = list(.mor = "weight", .likelihood = "likelihood")
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
        rownames_to_column("condition")
      )
  
  return(results)

}
