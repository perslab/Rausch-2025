# File paths
xenium_path <- here::here("sc_preprocessing/xenium/")
roi_file <- here::here("sc_preprocessing/cells_to_keep.csv")

analyze_xenium <- 
  list(
    tarchetypes::tar_files_input(
      xenium_paths,
      list.files(xenium_path, full.names = T)
    ),
    tar_target(
      roi_file,
      command = roi_file,
      format = "file"
    ),
    tar_target(
      cells_to_keep,
      read.csv(roi_file) %>% 
        janitor::clean_names() %>% 
        mutate(sample = gsub(".*_","",name), cellname = paste0(dat_cell_id, "_", sample)) %>% 
        column_to_rownames("cellname")
    ),
    tar_target(
      xenium_obj,
      gen_xenium(path=xenium_paths),
      pattern = map(xenium_paths),
      iteration = "list"
    ),
    tar_target(
      xenium_merged,
      merge(xenium_obj[[1]], xenium_obj[-1]) %>% 
        AddMetaData(., cells_to_keep) %>% 
        subset(., reg == "toremove", invert=T) %>%
        SCTransform(., method = "qpoisson", assay="Xenium") %>%
        RunPCA(.) %>% 
        RunUMAP(., dims = 1:30) %>% 
        FindNeighbors(., reduction = "pca", dims = 1:30) %>% 
        FindClusters(., resolution = 0.8)
    ),
    tar_target(
      xenium_mbh,
      subset(xenium_merged, subset = reg == "no_assign", invert=T) %>% 
        SCTransform(., method = "qpoisson", assay="Xenium") %>%
        RunPCA(.) %>% 
        RunUMAP(., dims = 1:30) %>% 
        FindNeighbors(., reduction = "pca", dims = 1:30) %>% 
        FindClusters(., resolution = 0.8)
    ),
    tar_target(
      xenium_arc,
      subset(xenium_merged, subset = reg == "ARC") %>% 
        SCTransform(., method = "qpoisson", assay="Xenium") %>%
        RunPCA(.) %>% 
        RunUMAP(., dims = 1:30) %>% 
        FindNeighbors(., reduction = "pca", dims = 1:30) %>% 
        FindClusters(., resolution = 0.8)
    ),
    tar_target(
      xenium_neurons,
      classify_cells_v2(xenium_subset, path = "/projects/dylan/target_testing/cell_markers/spatial_markers.txt") %>% 
        filter_cells_by_column(., column_name = "ann_neuron_v2", value = "Neurons") %>% 
        SCTransform(., method = "qpoisson", assay="Xenium") %>%
        RunPCA(.) %>% 
        RunUMAP(., dims = 1:30) %>% 
        FindNeighbors(., reduction = "pca", dims = 1:30) %>% 
        FindClusters(., resolution = 0.8)
    )
  )
