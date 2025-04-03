# Seurat processing params
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

#targets options
tar_option_set(error="null")

knockout_glp1 <- 
  list(
     tar_target(
      arc_glp1rko,
      {
        classify_cells(glp1rleprko_neurons, path = here::here("cell_markers/coarse_markers_vmh.txt"), column_name = "l2_arc_region",
                       clustering_res = clustering_res, clustering_dims = clustering_dims) %>% 
        filter_cells_by_column(., column_name = "l2_arc_region", value = "ARC_Neurons") %>%
        filter_cells_by_column(., column_name = "hash_pool", value = "H2303", invert=T) %>%
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) %>%
        project_umap(query = ., ref = arc_diet_labeled, dims = 30, label_to_transfer="seurat_clusters", reference.assay = "integrated", query.assay = "integrated") 
      }
    ),
    tar_target(
      glp1rko_to_remove,
      union(names(which(prop.table(table(arc_glp1rko$diet, arc_glp1rko$seurat_clusters), margin=2)[2,]<0.25)),
            names(which(tapply(arc_glp1rko$predicted.celltype.score, arc_glp1rko$seurat_clusters, mean)<0.7)))
    ),
    tar_target(
      arc_glp1rko_labeled,
      subset(arc_glp1rko, subset = seurat_clusters %in% glp1rko_to_remove, invert=T)
    ),
    tar_target(
      celltypes_glp1rko,
      names(table(arc_glp1rko_labeled$predicted.celltype))[table(arc_glp1rko_labeled$predicted.celltype) > 700]
    ),
    tar_target(
      glp1rko_celldist,
      {
        subset(arc_glp1rko_labeled, subset = predicted.celltype %in% celltypes_glp1rko) %>% 
        create_new_column_seurat(., selected_columns = c("predicted.celltype","diet"), new_column_name = "scdist_group") %>% 
          run_scdist(obj = ., fixed.effects = c("geno"), assay = "SCT", ident_column = "scdist_group", 
                   random.effects = c("orig.ident","hash.mcl.ID"), d = 20)
      }
    ),
    tar_target(
      pb_glp1rko,
      edger_prep(arc_glp1rko_labeled, celltype_column = "predicted.celltype", trt_group = "geno", celltypes = celltypes_glp1rko),
      pattern = map(celltypes_glp1rko),
      iteration = "list"
    ),
    tar_target(
      fit_glp1rko,
      run_edger(pb_glp1rko, mm = "~ geno + diet + geno:diet"),
      pattern = map(pb_glp1rko),
      iteration="list"
    ),
     tar_target(
      milo_glp1r,
      seurat_to_milo(arc_glp1rko_labeled, pca_reduction_name = "pca",
        umap_reduction_name = "umap",k = 30, d = 30, prop = 0.1, refinement_scheme = "graph", 
        sample_column = "orig.ident", additional_id_column = "hash.mcl.ID")
    )
  )