# R environment lockfile (renv)

`renv.lock` records the R packages the shiny_DPDP app and the DPDP pipeline use, at the versions on the lab Linux PC.
- **R 4.6.1**, **Bioconductor 3.23**, CRAN `https://cloud.r-project.org`
- **363 packages** (325 CRAN, 37 Bioconductor, plus renv 1.2.2): the 66 packages the app and pipeline reference, plus their dependencies.
- Made on Ubuntu 24.04 from the shared library `/home/R_Packages` (+ the app's bundled `shiny_DPDP/Packages/pca3d`).

## Restore on a new machine
1. Install R 4.6.1 and the system libraries (Ubuntu 24.04; from `renv::sysreqs()`):
   ```bash
   sudo apt install -y cmake git gsfonts libcairo2-dev libcurl4-openssl-dev libfontconfig1-dev libfreetype6-dev \
     libgit2-dev libglpk-dev libicu-dev libmagick++-dev libnode-dev libpng-dev libssl-dev libx11-dev libxml2-dev \
     make pandoc perl python3
   ```
2. Restore into a library of your choice (most packages build from source; expect a long first run):
   ```r
   install.packages("renv")
   renv::restore(lockfile = "renv.lock", library = "/path/to/R_library", prompt = FALSE)
   ```
   Then point R at that library (for example `R_LIBS_USER=/path/to/R_library` in `~/.Renviron`). Run the pipeline with plain
   `Rscript`, never `--vanilla`, or the library setting is skipped.

Tested 2026-11-04: restoring `xtable`, `pca3d` (archived on CRAN) and `BiocGenerics` with their dependencies (33 packages)
into an empty library worked.

## Not in the lockfile
- `rmotifx`, `PTMphinder`: referenced by the app's Phos motif step only; not installed on the lab PC (GitHub packages).

## Update
After installing or updating packages on the lab PC: `Rscript tools/renv/snapshot.R` (from the project root), then rebuild the plugin
(`tools/build_plugin.sh`).
