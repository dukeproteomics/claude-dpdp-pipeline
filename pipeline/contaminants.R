# Contaminant and spike-in flags (lab rule, user 2026-11-04: flag contaminants, don't remove them).
# Same rule as the normalization exclude default: added material only.
#   contaminant: trypsin   digestion trypsin from pig or cow (not the sample organism's own trypsins or trypsin inhibitors)
#   contaminant: keratin   human keratins, when the sample organism is not human (in human samples they can't be told apart)
#   spike-in: <name>       accessions of the QC spikes (casein by level, ADH constant), from the QC step settings
# Returns a character vector ("" = not flagged), one per row.

sample_organism <- function(description) {
  org <- sub(".*org=([^;]*)$", "\\1", description[grepl("org=", description)])
  if (!length(org)) return(NA_character_)
  names(sort(table(org), decreasing = TRUE))[1]
}

flag_contaminants <- function(accession, description, spike_acc = c("P02662", "P02663", "P02666"),
                              constant_acc = c("P00330", "P00331"), sample_org = sample_organism(description)) {
  d <- tolower(description)
  flag <- rep("", length(d))
  human <- !is.na(sample_org) && grepl("homo sapiens", sample_org, ignore.case = TRUE)
  flag[grepl("trypsin, org=(sus scrofa|bos taurus)", d)] <- "contaminant: trypsin"
  if (!human) flag[grepl("keratin, type [^,]*, org=homo sapiens", d)] <- "contaminant: human keratin"
  first_acc <- sub(";.*", "", accession)
  flag[first_acc %in% spike_acc] <- "spike-in: casein"
  flag[first_acc %in% constant_acc] <- "spike-in: ADH"
  flag
}

# Spike accessions recorded by the QC step, else the lab defaults.
qc_spike_settings <- function(run_dir) {
  f <- file.path(run_dir, "steps", "09_qc", "summary.json")
  s <- if (file.exists(f)) jsonlite::fromJSON(f)$settings else NULL
  list(spike = union(unlist(s$spike_accession), c("P02662", "P02663", "P02666")),
       constant = unlist(s$constant_accession) %||% c("P00330", "P00331"))
}
