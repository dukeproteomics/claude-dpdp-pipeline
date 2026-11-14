# STRING (string-db.org) API client for the Pathway step (user, 2026-11-13: the STRING website is fine to use).
# Only accessions / STRING ids are sent (no intensities). Calls go to the version's stable address (e.g. version-12-5) so a run
# can be repeated against the same data. STRING asks for one call at a time with a short pause; caller_identity is generic.
# https://string-db.org/help/api/

STRING_CALLER <- "dpdp_pipeline"

string_version <- function() {
  r <- httr::GET("https://string-db.org/api/json/version", httr::timeout(30))
  httr::stop_for_status(r)
  v <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))
  list(version = v$string_version[1], base = v$stable_address[1])
}

string_post <- function(base, path, body, tries = 3) {
  body$caller_identity <- STRING_CALLER
  for (k in seq_len(tries)) {
    r <- tryCatch(httr::POST(paste0(base, "/api/json/", path), body = body, encode = "form", httr::timeout(180)), error = function(e) e)
    Sys.sleep(1)   # STRING: one request at a time, short pause
    if (!inherits(r, "error") && httr::status_code(r) == 200) {
      txt <- httr::content(r, "text", encoding = "UTF-8")
      return(if (nzchar(txt)) jsonlite::fromJSON(txt) else data.frame())
    }
    msg <- if (inherits(r, "error")) conditionMessage(r) else paste("HTTP", httr::status_code(r), substr(httr::content(r, "text", encoding = "UTF-8"), 1, 200))
    if (k == tries) stop("STRING ", path, " failed: ", msg, call. = FALSE)
    Sys.sleep(5 * k)
  }
}

# Accessions (first of a protein group) -> STRING ids for one species, in chunks. Returns data.frame(query, string_id, name).
string_map <- function(base, accessions, species, chunk = 1000) {
  acc <- unique(accessions[!is.na(accessions) & nzchar(accessions)])
  out <- list()
  for (s in split(acc, ceiling(seq_along(acc) / chunk))) {
    r <- string_post(base, "get_string_ids", list(identifiers = paste(s, collapse = "\r"), species = species, limit = 1, echo_query = 1))
    if (NROW(r)) out[[length(out) + 1]] <- data.frame(query = r$queryItem, string_id = r$stringId, name = r$preferredName, stringsAsFactors = FALSE)
  }
  m <- if (length(out)) do.call(rbind, out) else data.frame(query = character(0), string_id = character(0), name = character(0))
  m[!duplicated(m$query), ]
}

# Enrichment of ids against a background of ids (both STRING ids). Returns STRING's table (category, term, description, counts, p, fdr, genes).
string_enrich <- function(base, ids, background, species) {
  r <- string_post(base, "enrichment", list(identifiers = paste(ids, collapse = "\r"), background_string_identifiers = paste(background, collapse = "\r"),
                                            species = species))
  if (!NROW(r)) return(data.frame())
  r
}

string_network <- function(base, ids, species, required_score = 400, network_type = "functional") {
  r <- string_post(base, "network", list(identifiers = paste(ids, collapse = "\r"), species = species, required_score = required_score,
                                          network_type = network_type))
  if (!NROW(r)) return(data.frame(a = character(0), b = character(0), score = numeric(0)))
  data.frame(a = r$stringId_A, b = r$stringId_B, score = r$score, stringsAsFactors = FALSE)
}

string_link <- function(base, ids, species, required_score = 400, network_type = "functional") {
  r <- tryCatch(string_post(base, "get_link", list(identifiers = paste(ids, collapse = "\r"), species = species, required_score = required_score,
                                                   network_flavor = "confidence", network_type = network_type)), error = function(e) NULL)
  if (is.null(r) || !length(r)) return(NA_character_)
  as.character(unlist(r))[1]
}

# Organisms: from "org=..." (Spectronaut) or "OS=..." (UniProt) in the descriptions; NCBI taxon ids STRING uses.
STRING_TAXA <- c("Homo sapiens" = 9606, "Mus musculus" = 10090, "Rattus norvegicus" = 10116, "Saccharomyces cerevisiae" = 4932,
                 "Escherichia coli" = 511145, "Cryptococcus neoformans" = 214684, "Cryptococcus deneoformans" = 214684,
                 "Arabidopsis thaliana" = 3702, "Drosophila melanogaster" = 7227, "Danio rerio" = 7955, "Caenorhabditis elegans" = 6239,
                 "Bos taurus" = 9913, "Sus scrofa" = 9823, "Gallus gallus" = 9031, "Macaca mulatta" = 9544, "Xenopus laevis" = 8355,
                 "Schizosaccharomyces pombe" = 4896, "Bacillus subtilis" = 224308, "Pseudomonas aeruginosa" = 208964,
                 "Staphylococcus aureus" = 93061, "Mycobacterium tuberculosis" = 83332, "Candida albicans" = 237561, "Plasmodium falciparum" = 36329)
organism_of <- function(desc) {
  o <- ifelse(grepl("org=", desc), sub(".*org=([^,;]+).*", "\\1", desc), ifelse(grepl("OS=", desc), sub(".*OS=(.+?) OX=.*", "\\1", desc), NA))
  o <- trimws(sub(" \\(.*$", "", o))   # "Escherichia coli (strain K12)" -> "Escherichia coli"
  vapply(strsplit(ifelse(is.na(o), "", o), " "), function(w) paste(head(w, 2), collapse = " "), "")
}

# STRING's annotation of every protein of a species (the terms the enrichment uses), for the ranked analysis (GSEA). Downloaded once
# per species into cache_dir (shared by the runs of a data folder). Returns data.frame(id, category, term, description) with the
# category as the enrichment API names it (Process, Function, ... ). KEGG is not in STRING's download files (license).
STRING_FILE_CATEGORIES <- c("Biological Process (Gene Ontology)" = "Process", "Molecular Function (Gene Ontology)" = "Function",
  "Cellular Component (Gene Ontology)" = "Component", "Reactome Pathways" = "RCTM", "WikiPathways" = "WikiPathways",
  "Annotated Keywords (UniProt)" = "Keyword", "Local Network Cluster (STRING)" = "NetworkNeighborAL", "Protein Domains (Pfam)" = "Pfam",
  "Protein Domains and Features (InterPro)" = "InterPro", "Protein Domains (SMART)" = "SMART", "Subcellular localization (COMPARTMENTS)" = "COMPARTMENTS",
  "Tissue expression (TISSUES)" = "TISSUES", "Disease-gene associations (DISEASES)" = "DISEASES", "KEGG Pathways" = "KEGG")
string_terms <- function(species, cache_dir, categories, version = "12.0") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache_dir, sprintf("%s.protein.enrichment.terms.v%s.txt.gz", species, version))
  if (!file.exists(f) || file.size(f) < 1000) {
    url <- sprintf("https://stringdb-downloads.org/download/protein.enrichment.terms.v%s/%s.protein.enrichment.terms.v%s.txt.gz", version, species, version)
    cat("Downloading STRING annotation for species", species, "...\n")
    tmp <- paste0(f, ".part"); r <- httr::GET(url, httr::write_disk(tmp, overwrite = TRUE), httr::timeout(900))
    if (httr::status_code(r) != 200) { unlink(tmp); stop("STRING annotation download failed (HTTP ", httr::status_code(r), "): ", url, call. = FALSE) }
    file.rename(tmp, f)
  }
  d <- data.table::fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  names(d) <- c("id", "category", "term", "description")
  d$category <- unname(STRING_FILE_CATEGORIES[d$category])
  as.data.frame(d[!is.na(d$category) & d$category %in% categories, ])
}

# GO Consortium's generic GO slim (broad categories), downloaded once into cache_dir. Returns data.frame(term, name, ont) with
# ont BP / MF / CC; the three roots are left out.
goslim_generic <- function(cache_dir) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache_dir, "goslim_generic.obo")
  if (!file.exists(f) || file.size(f) < 1000) {
    cat("Downloading the generic GO slim...\n")
    r <- httr::GET("https://current.geneontology.org/ontology/subsets/goslim_generic.obo", httr::write_disk(paste0(f, ".part"), overwrite = TRUE), httr::timeout(120))
    if (httr::status_code(r) != 200) { unlink(paste0(f, ".part")); stop("GO slim download failed (HTTP ", httr::status_code(r), ")", call. = FALSE) }
    file.rename(paste0(f, ".part"), f)
  }
  x <- readLines(f, warn = FALSE)
  starts <- which(x == "[Term]"); ends <- c(starts[-1] - 1, length(x))
  t <- do.call(rbind, lapply(seq_along(starts), function(k) { b <- x[starts[k]:ends[k]]
    g <- function(tag) sub(paste0("^", tag, ": "), "", grep(paste0("^", tag, ": "), b, value = TRUE)[1])
    data.frame(term = g("id"), name = g("name"), ns = g("namespace"), obsolete = any(b == "is_obsolete: true"), stringsAsFactors = FALSE) }))
  t <- t[!t$obsolete & grepl("^GO:", t$term) & !t$term %in% c("GO:0008150", "GO:0003674", "GO:0005575"), ]
  t$ont <- c(biological_process = "BP", molecular_function = "MF", cellular_component = "CC")[t$ns]
  t[!is.na(t$ont), c("term", "name", "ont")]
}

# GO annotations (id, category, term) with every ancestor term added (GO.db), so a protein counts in each broader term.
go_propagate <- function(an) {
  if (!NROW(an)) return(an)
  anc <- list(Process = GO.db::GOBPANCESTOR, Function = GO.db::GOMFANCESTOR, Component = GO.db::GOCCANCESTOR)
  out <- lapply(names(anc), function(cat_) { a <- an[an$category == cat_, c("id", "term")]; if (!nrow(a)) return(NULL)
    u <- unique(a$term); m <- AnnotationDbi::mget(u, anc[[cat_]], ifnotfound = NA)
    lk <- data.frame(term = rep(u, lengths(m)), anc = unlist(m, use.names = FALSE), stringsAsFactors = FALSE)
    lk <- lk[!is.na(lk$anc) & lk$anc != "all", ]
    e <- merge(a, lk, by = "term")
    r <- unique(rbind(a, data.frame(id = e$id, term = e$anc, stringsAsFactors = FALSE)))
    cbind(r, category = cat_) })
  do.call(rbind, out)[, c("id", "category", "term")]
}
