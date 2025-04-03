source("prep_seurat.R")
source("integrate_objects.R")
source("analyze_diets.R")
source("analyze_lepip.R")
source("analyze_agrpko.R")
source("analyze_xenium.R")
source("analyze_glp1rfl.R")

tar_option_set(
  packages = c("tidyverse", "Seurat", "miloR"), # packages that your targets need to run
  format = "qs", # default storage format,
  error = "null",
  retrieval = "worker",
  storage = "worker",
  priority = 0.1,
  memory = "transient",
  garbage_collection = TRUE,
  cue = tar_cue(
    mode = c("thorough", "always", "never"),
    command = TRUE,
    depend = TRUE,
    format = TRUE,
    repository = TRUE,
    iteration = TRUE,
    file = TRUE
  )
)

options(clustermq.scheduler = "multicore",
        clustermq.ssh.timeout=36000,
        clustermq.worker.timeout=36000,
        clustermq.error.timeout=36000,
        clustermq.ssh.log='clustermq_sshlog.log'
)


list(
  seurat_processing,
  integrate_objects,
  gen_diet,
  analyze_diet,
  gen_lepip,
  analyze_lepip,
  knockout_glp1,
  knockout_agrp
  )

#knockout_glp1,
#knockout_agrp