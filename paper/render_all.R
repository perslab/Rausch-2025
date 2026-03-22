#!/usr/bin/env Rscript
#
# Master render script for all paper figures and tables.
#
# Usage (standard environment):
#   module load gcc/11.2.0 R/4.3.3
#   Rscript paper/render_all.R
#
# Figure 1 (spatial) requires a separate environment:
#   module purge
#   module load gcc/13.2.0 openjdk/20.0.0 openssl/3.1.0 curl/8.6.0 sqlite3/3.44.2 \
#     python/3.11.3 tiff/4.5.0 geos/3.11.1 swig/4.2.1 hdf5/1.12.1 proj/7.0.1 \
#     gdal/3.5.3 udunits2/2.2.28 R/4.5.0
#   R_LIBS_USER=renv_spatial/lib Rscript -e "rmarkdown::render('paper/figure1.Rmd')"

library(rmarkdown)

setwd(here::here())
paper_dir <- "paper"
fi_dir    <- "paper/food_intake_analysis"

render_notebook <- function(path, label = basename(path)) {
    cat(sprintf("\n========== %s ==========\n", label))
    tryCatch({
        rmarkdown::render(path, quiet = TRUE)
        cat("  OK\n")
    }, error = function(e) {
        cat("  FAILED:", e$message, "\n")
    })
}

# --- Phase 1: LGS derivation (no dependencies) ---
render_notebook(file.path(paper_dir, "lgs.Rmd"))

# --- Phase 2: Main figures ---
# Figure 1 (spatial) — SKIPPED: requires R 4.5.0 + spatial env (see header)
cat("\n========== figure1.Rmd ==========\n")
cat("  SKIPPED: requires spatial environment (R 4.5.0). Run manually.\n")

render_notebook(file.path(paper_dir, "figure2.Rmd"))
render_notebook(file.path(paper_dir, "figure4.Rmd"))

# --- Phase 3: Food intake scripts ---
render_notebook(file.path(fi_dir, "fig3_dreadd.Rmd"))
render_notebook(file.path(fi_dir, "fig4_leptin_refeed.Rmd"))
render_notebook(file.path(fi_dir, "fig5_hfd_onset.Rmd"))
render_notebook(file.path(fi_dir, "fig5_desensitization.Rmd"))
render_notebook(file.path(fi_dir, "ed_body_composition.Rmd"))

# --- Phase 4: Extended data figures ---
render_notebook(file.path(paper_dir, "ed_figure2.Rmd"))
render_notebook(file.path(paper_dir, "ed_figure3.Rmd"))
render_notebook(file.path(paper_dir, "ed_figure4.Rmd"))
render_notebook(file.path(paper_dir, "ed_figure5.Rmd"))
render_notebook(file.path(paper_dir, "ed_lgs.Rmd"))
render_notebook(file.path(paper_dir, "ed_markers.Rmd"))
render_notebook(file.path(paper_dir, "ed_overlaps.Rmd"))

# --- Phase 5: Supplemental tables ---
render_notebook(file.path(paper_dir, "supp_tables.Rmd"))

cat("\n========== DONE ==========\n")
cat("All outputs in: output_figures/paper/\n")
cat("Spatial figure1 must be rendered separately in the R 4.5.0 environment.\n")
