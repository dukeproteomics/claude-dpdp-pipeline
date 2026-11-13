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
