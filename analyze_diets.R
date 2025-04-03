# Seurat processing params
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 30

# targets options
tar_option_set(error = "null")

gen_diet <-
  list(
    tar_target(
      arc_diet,
      {
        classify_cells(diet_neurons, path = here::here("cell_markers/coarse_markers_vmh.txt"), column_name = "l2_arc_region",
                       clustering_res = clustering_res, clustering_dims = clustering_dims) %>% 
        filter_cells_by_column(., column_name = "l2_arc_region", value = "ARC_Neurons") %>%
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) %>%
         map_ref(., dims=30, column_name = "internal_labs", ref = inhouse_scrna)
      }
    ),
    tar_target(
      diet_to_remove,
      names(which(tapply(arc_diet$prediction.score.max, arc_diet$seurat_clusters, mean)<0.7))
    ),
    tar_target(
      arc_diet_filtered,
       subset(arc_diet, idents = diet_to_remove, invert=T) %>% 
          filter_cells_by_column(., column_name = "internal_labs", 
                                 values = c("Lpar1_oligo","Sim1/Rprm","Sim1/Ebf3"), invert = T) %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims, return_model=T) 
    ),
     tar_target(
      arc_diet_labeled,
      arc_diet_filtered %>% modify_column(., "seurat_clusters", c("1","14"),"Agrp")
    )
  )

analyze_diet <-
  list(    
    tar_target(
      celltypes_diet,
      names(table(arc_diet_labeled$seurat_clusters))[table(arc_diet_labeled$seurat_clusters) > 200]
    ),
    tar_target(
      hfdveh_celldist,
      {
        subset(arc_diet_labeled, subset = treatment %in% c("Chow", "HFD") & seurat_clusters %in% celltypes_diet) %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "hash_pool"), assay = "SCT", ident_column = "seurat_clusters", 
                   random.effects = c("orig.ident","hash.mcl.ID"), d = 20)
      }
    ),
    tar_target(
      fastveh_celldist,
      {
        subset(arc_diet_labeled, subset = treatment %in% c("Chow", "Fast") & seurat_clusters %in% celltypes_diet) %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "hash_pool"), assay = "SCT", ident_column = "seurat_clusters", 
                     random.effects = c("orig.ident","hash.mcl.ID"), d = 20)
      }
    ),
    tar_target(
      celltypes_refeed,
      names(table(arc_diet_labeled$seurat_clusters))[table(arc_diet_labeled$seurat_clusters[arc_diet_labeled$treatment=="Refeed"]) > 40]
    ),
    tar_target(
      refeed_celldist,
      {
        subset(arc_diet_labeled, subset = treatment %in% c("Refeed", "Fast") & seurat_clusters %in% celltypes_refeed) %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "hash_pool"), assay = "SCT", ident_column = "seurat_clusters", 
                     random.effects = c("orig.ident","hash.mcl.ID"), d = 20)
      }
    ),
    tar_target(
      pb_alldiet,
      edger_prep(arc_diet_labeled, celltype_column = "seurat_clusters", trt_group = "treatment", celltypes = celltypes_diet),
      pattern = map(celltypes_diet),
      iteration = "list"
    ),
    tar_target(
      milo_diet,
      seurat_to_milo(arc_diet_labeled, pca_reduction_name = "pca",
        umap_reduction_name = "umap",k = 30,d = 30,prop = 0.1,refinement_scheme = "graph", 
        sample_column = "orig.ident",additional_id_column = "hash.mcl.ID")
    )
  )