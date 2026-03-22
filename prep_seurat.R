## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)
options(tidyverse.quiet = TRUE)

# File paths
raw_path <- here::here("sc_preprocessing/raw_h5_files")
denoised_path <- here::here("sc_preprocessing/denoised_counts")
spatial_folder <-  here::here("sc_preprocessing/xenium/diet")
metadata_file <- here::here("sc_preprocessing/sequenced_mice_data_swapped_w_all_160424.xlsx")
marker_file <- here::here("cell_markers/coarse_markers.txt")

# Seurat processing params
ann_column <- "ann_neuron"
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

# Step 2: Reference Parameters in Targets

seurat_processing <-
  list(
    tarchetypes::tar_files_input(
      paths,
      Sys.glob(file.path(denoised_path, "*_out.h5")), 
      cue = tar_cue(mode = "always")
    ),
    tarchetypes::tar_files_input(
      spatial,
      Sys.glob(file.path(spatial_folder, "*output*"))
    ),
    tar_target(
      metadata,
      command = metadata_file,
      format = "file"
    ),
    tar_target(
      cellmarkers,
      command = marker_file,
      format = "file"
    ),
    tar_target(
      spatial_objects,
      LoadXeniumNew(spatial),
      pattern = map(spatial),
      iteration = "list"
    ),
    tar_target(
      seurat_objects,
      {
        gen_seurat(denoisedfile = paths, rawpath = raw_path) %>% 
          add_hto() %>% 
          demux_hto()
      },
      pattern = map(paths),
      iteration = "list"
    ),
    tar_target(
      filtered_seurat,
      {
        seurat_objects %>%
          remove_doublets() %>%
          add_metadata(obj = ., file = metadata, joinby = "hash") %>%
          classify_cells(obj = ., path = marker_file, clustering_res = clustering_res,
                         clustering_dims = clustering_dims, column_name = "l1_neur") %>%
          run_gex_qc(., times=1, batch = "l1_neur") %>%
          process_seurat(., method = "log", res = clustering_res, dims = clustering_dims)
      },
      pattern = map(seurat_objects),
      iteration = "list"
    )
  )
