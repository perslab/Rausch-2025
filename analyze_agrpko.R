# Seurat processing params
ann_column <- "ann_neuron"
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

knockout_agrp <- 
  list(
    tar_target(
      arc_agrpko,
      classify_cells(agrpleprko_neurons, path = here::here("cell_markers/coarse_markers_vmh.txt"), column_name = "l2_arc_region",
                       clustering_res = clustering_res, clustering_dims = clustering_dims) %>% 
        filter_cells_by_column(., column_name = "genetic", values = "Cre neg; Flox WT", invert=T) %>% 
        filter_cells_by_column(., column_name = "l2_arc_region", value = "ARC_Neurons") %>%
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) %>%
        project_umap(query = ., ref = arc_diet_labeled, dims = 30, label_to_transfer="seurat_clusters", reference.assay = "integrated", query.assay = "integrated") 
    ),
    tar_target(
      agrpko_to_remove,
      union(names(which(colSums(prop.table(table(arc_agrpko$hash.mcl.ID, arc_agrpko$seurat_clusters), margin=2)<0.02)>2)),
            names(which(tapply(arc_agrpko$predicted.celltype.score, arc_agrpko$seurat_clusters, mean)<0.8)))
    ),
    tar_target(
      arc_agrpko_labeled,
      subset(arc_agrpko, subset = seurat_clusters %in% agrpko_to_remove, invert=T)
    ),
    tar_target(
      celltypes_agrpko,
      names(table(arc_agrpko_labeled$predicted.celltype))[table(arc_agrpko_labeled$predicted.celltype) > 300]
    ),
    tar_target(
      agrpko_celldist,
      {
        subset(arc_agrpko_labeled, subset = predicted.celltype %in% celltypes_agrpko) %>% 
        create_new_column_seurat(., selected_columns = c("predicted.celltype","diet"), new_column_name = "scdist_group") %>% 
          run_scdist(obj = ., fixed.effects = c("geno"), assay = "SCT", ident_column = "scdist_group", 
                   random.effects = c("orig.ident","hash.mcl.ID"), d = 10)
      }
    ),
    tar_target(
      pb_agrpko,
      edger_prep(arc_agrpko_labeled, celltype_column = "predicted.celltype", trt_group = "geno", celltypes = celltypes_agrpko),
      pattern = map(celltypes_agrpko),
      iteration = "list"
    ),
    tar_target(
      fit_agrpko,
      run_edger(pb_agrpko, mm = "~ geno + diet + geno:diet"),
      pattern = map(pb_agrpko),
      iteration="list"
    )
  )