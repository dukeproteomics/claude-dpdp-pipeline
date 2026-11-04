# Regenerate tools/renv/renv.lock: every package the shiny_DPDP app and the pipeline reference (renv::dependencies),
# plus their dependencies, at the versions installed on this machine. Does not create an renv project or touch the app.
# Run from the project root:  Rscript tools/renv/snapshot.R
d <- renv::dependencies(c("shiny_DPDP", "pipeline"), progress = FALSE, quiet = TRUE)
base <- rownames(installed.packages(priority = "base"))
libs <- c(.libPaths(), normalizePath("shiny_DPDP/Packages"))   # the app bundles pca3d in shiny_DPDP/Packages
wanted <- c(setdiff(sort(unique(d$Package)), c(base, "R")), "renv")
inst <- rownames(installed.packages(lib.loc = libs))
missing <- setdiff(wanted, inst)
if (length(missing)) message("Referenced but not installed (not in the lockfile): ", paste(missing, collapse = ", "))
renv::snapshot(project = "tools/renv", library = libs, lockfile = "tools/renv/renv.lock",
               packages = intersect(wanted, inst), prompt = FALSE, force = TRUE)
