# Seurat processing params
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

#targets options
tar_option_set(error="null")

gen_lepip <- 
  list(
    tar_target(
      arc_ob_db,
      {
        classify_cells(ob_db_neuron, path = here::here("cell_markers/coarse_markers_vmh.txt"), column_name = "l2_arc_region", clustering_res = clustering_res,
                       clustering_dims = clustering_dims) %>% 
          filter_cells_by_column(., column_name = "l2_arc_region", value = "ARC_Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch,
                         res = clustering_res, dims = clustering_dims) %>% 
           project_umap(query = ., ref = arc_diet_labeled, dims = 30, label_to_transfer="seurat_clusters", reference.assay = "integrated", query.assay = "integrated")
      }
    ),
    tar_target(
      lepip_to_remove,
      union(names(which(prop.table(table(arc_ob_db$dataset, arc_ob_db$seurat_clusters), margin=2)[2,]>0.4)),
            names(which(tapply(arc_ob_db$predicted.celltype.score, arc_ob_db$seurat_clusters, mean)<0.7)))
    ),
    tar_target(
      arc_lepip,
      {
       subset(arc_ob_db, subset = seurat_clusters %in% lepip_to_remove, invert=T) %>%
          create_new_column_seurat(., selected_columns = c("geno","treatment"), new_column_name = "group") 
      }
    )
  )

analyze_lepip <- 
  list(
    tar_target(
      celltypes,
      names(table(arc_lepip$predicted.celltype))[table(arc_lepip$predicted.celltype) > 200]
    ),
    tar_target(
      pseudobulk_lepip, 
      edger_prep(arc_lepip, dataset = "lepip", celltype_column="predicted.celltype", trt_group = "group", celltypes = celltypes),
      pattern = map(celltypes),
      iteration="list"
    ),
    tar_target(
      obwt_celldist, 
      {
        filter_cells_by_column(arc_lepip, column_name = "treatment", values = "Sal") %>% 
          subset(., subset = predicted.celltype %in% celltypes) %>% 
          run_scdist(obj = ., fixed.effects = c("geno","dataset"),  
                     assay="SCT", ident_column = "predicted.celltype",
                     random.effects = c("orig.ident","hash.mcl.ID"), 
                     d = 20)
      }
    ),
    tar_target(
      lepip_celldist,
      {
        subset(arc_lepip, subset = time %in% c(1,3,6) & predicted.celltype %in% celltypes & dataset == "lepip") %>%
          create_new_column_seurat(., selected_columns = c("predicted.celltype","geno"), new_column_name = "scdist_group") %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "seq_pool","time"),
                     assay="SCT", ident_column = "scdist_group", 
                     random.effects = c("hash.mcl.ID","orig.ident"),
                     d = 20)
      }
    ),
    tar_target(
      myers_scrna.map,
      prep_myers(here::here("external_data/myers_sctrap/"), min.features = 1000, min.cells=10) %>%
        project_umap(arc_lepip_labeled, ., reference.assay="integrated", query.assay="RNA", label_to_transfer = "seurat_clusters")
    ),
    tar_target(
      campbell_scrna.map,
      readRDS(here::here("external_data/Arc_Neurons_33_Clusters.rds")) %>%
        project_umap(arc_lepip_labeled, ., reference.assay="integrated", query.assay="integrated", label_to_transfer = "seurat_clusters") %>%
        process_seurat(., method="qpoisson", cluster=F)
    ),
    tar_target(
      focus_object,
      subset(arc_lepip_labeled, subset = seurat_clusters %in% c("0", "9", "26"))
    )
  )