# Step 12 Phospho: sites, sequence windows, enrichment logos and motif-x (pipeline; user 2026-11-13, replacing rmotifx / MoMo).
#
# Sites: every localized site of the tested phosphopeptides (Local2 = "Y"; Protein_PTM_Loc gives residue + protein position), one
# row per site; a site on several peptides takes the most significant one. Windows: the 15 residues around the site (+-7) from the
# UniProt sequence of the first accession of the protein group ('_' past the ends); the site must be that residue in the sequence.
# Background for logos and motif-x: every tested site of the comparison (what the experiment measured), not the proteome.

AA20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

# UniProt sequences by accession (rest.uniprot.org, batches of 400), cached in a FASTA file. Returns a named character vector.
uniprot_sequences <- function(accessions, cache) {
  have <- if (file.exists(cache)) read_fasta(cache) else character(0)
  need <- setdiff(unique(accessions[!is.na(accessions) & nzchar(accessions)]), names(have))
  for (b in split(need, ceiling(seq_along(need) / 400))) {
    url <- paste0("https://rest.uniprot.org/uniprotkb/accessions?format=fasta&size=500&accessions=", paste(b, collapse = ","))
    r <- NULL
    for (k in 1:3) { r <- tryCatch(httr::GET(url, httr::timeout(120)), error = function(e) NULL); if (!is.null(r) && httr::status_code(r) == 200) break; Sys.sleep(3 * k) }
    if (is.null(r) || httr::status_code(r) != 200) stop("UniProt request failed (", if (is.null(r)) "no answer" else httr::status_code(r), ")", call. = FALSE)
    cat(httr::content(r, "text", encoding = "UTF-8"), file = cache, append = TRUE)
    Sys.sleep(0.5)
  }
  read_fasta(cache)
}
read_fasta <- function(f) {
  l <- readLines(f, warn = FALSE); h <- grepl("^>", l)
  if (!any(h)) return(character(0))
  acc <- sub("^>[a-z]+\\|([^|]+)\\|.*", "\\1", l[h]); id <- cumsum(h)
  sq <- tapply(l[!h], id[!h], paste, collapse = "")
  s <- setNames(as.character(sq[as.character(seq_along(acc))]), acc); s[!duplicated(names(s)) & !is.na(s)]
}

# One row per site: acc, gene, site (e.g. S48), residue, pos, fc, p, call, peptides (rows carrying it), localized
phospho_sites <- function(st, nm, hits = "stats calls", p_type = "adjusted", p_cut = 0.05, fc_cut = 1.5, localized_only = TRUE) {
  keep <- !is.na(st$Protein_PTM_Loc) & nzchar(st$Protein_PTM_Loc)
  if (localized_only && !is.null(st$Local2)) keep <- keep & st$Local2 %in% "Y"
  st <- st[keep, , drop = FALSE]
  fc <- comp_get(st, nm, "_FC"); praw <- comp_get(st, nm, "_pval"); padj <- comp_get(st, nm, "_adjpval") %||% praw
  call <- if (identical(hits, "custom")) { pv <- if (identical(p_type, "raw")) praw else padj
    ifelse(!is.na(pv) & pv <= p_cut & fc >= fc_cut, "Up", ifelse(!is.na(pv) & pv <= p_cut & fc <= -fc_cut, "Down", "")) } else ifelse(is.na(st$Stats), "", st$Stats)
  locs <- strsplit(gsub("\\s", "", st$Protein_PTM_Loc), "[,;]")
  n <- lengths(locs)
  d <- data.frame(acc = rep(sub(";.*", "", st$Accession), n), gene = rep(ifelse(is.na(st$Genes), "", st$Genes), n), site = unlist(locs),
                  fc = rep(fc, n), p = rep(praw, n), call = rep(call, n), stringsAsFactors = FALSE)
  d <- d[grepl("^[STY][0-9]+$", d$site), ]
  d$residue <- substr(d$site, 1, 1); d$pos <- as.integer(substring(d$site, 2))
  d$key <- paste(d$acc, d$site)
  pep <- table(d$key)
  d <- d[order(is.na(d$p), d$p), ]; d <- d[!duplicated(d$key), ]
  d$peptides <- as.integer(pep[d$key])
  d
}

# 15-mer window, site in lower case, '_' past the ends; NA when the sequence is missing or the residue doesn't match
site_window <- function(seq, pos, residue, flank = 7) {
  vapply(seq_along(pos), function(i) { s <- seq[i]
    if (is.na(s) || pos[i] > nchar(s) || substr(s, pos[i], pos[i]) != residue[i]) return(NA_character_)
    w <- vapply((pos[i] - flank):(pos[i] + flank), function(j) if (j < 1 || j > nchar(s)) "_" else substr(s, j, j), "")
    w[flank + 1] <- tolower(w[flank + 1]); paste(w, collapse = "") }, "")
}

# Enrichment logo (pLogo-style, O'Shea et al. 2013): for each position and residue, the binomial significance of its count in the
# foreground given its frequency in the background; + over-represented, - under-represented (log10). Bonferroni line returned.
logo_stats <- function(fg, bg, flank = 7) {
  pos <- setdiff(seq_len(2 * flank + 1), flank + 1)
  M <- matrix(0, length(AA20), length(pos), dimnames = list(AA20, as.character(pos - flank - 1)))
  F <- M
  if (length(fg) < 3 || length(bg) < 10) return(NULL)
  for (j in seq_along(pos)) {
    a <- toupper(substr(fg, pos[j], pos[j])); b <- toupper(substr(bg, pos[j], pos[j]))
    for (r in AA20) { k <- sum(a == r); n <- length(a); q <- (sum(b == r) + 0.5) / (length(b) + 1)
      F[r, j] <- k / n
      up <- pbinom(k - 1, n, q, lower.tail = FALSE); dn <- pbinom(k, n, q)
      M[r, j] <- if (k / n >= q) -log10(max(up, 1e-300)) else log10(max(dn, 1e-300)) }
  }
  list(score = round(M, 2), freq = round(F, 3), positions = as.integer(colnames(M)), residues = AA20, n_fg = length(fg), n_bg = length(bg),
       threshold = round(-log10(0.05 / (length(AA20) * length(pos))), 2))
}

# motif-x (Schwartz & Gygi, Nat Biotechnol 2005), re-implemented: iteratively fix the position/residue pair whose count in the
# foreground is the least likely under the background frequency (binomial), while p < sig and the count >= min_occ; the fixed
# pattern is a motif; remove its matches from both sets and search again. Foreground / background: windows with the same central
# residue. Returns motifs with their counts and fold increase.
motif_x <- function(fg, bg, central, min_occ = 20, sig = 1e-6, flank = 7, max_motifs = 25) {
  up <- function(x) toupper(x)
  fg <- up(fg[substr(up(fg), flank + 1, flank + 1) == central]); bg <- up(bg[substr(up(bg), flank + 1, flank + 1) == central])
  n_fg0 <- length(fg); n_bg0 <- length(bg); out <- list()
  if (n_fg0 < min_occ || n_bg0 < 10) return(data.frame())
  chars <- function(s, j) substr(s, j, j)
  pos <- setdiff(seq_len(2 * flank + 1), flank + 1)
  while (length(fg) >= min_occ && length(out) < max_motifs) {
    fixed <- list(); f <- fg; b <- bg; score <- 0
    repeat {
      best <- NULL
      for (j in setdiff(pos, as.integer(names(fixed)))) {
        a <- chars(f, j); bb <- chars(b, j); ta <- table(a); tb <- table(bb)
        for (r in intersect(names(ta), AA20)) { k <- ta[[r]]; if (k < min_occ) next
          q <- (if (r %in% names(tb)) tb[[r]] else 0) / length(b); if (q <= 0) q <- 0.5 / length(b)
          pv <- pbinom(k - 1, length(f), q, lower.tail = FALSE)
          if (pv < sig && (is.null(best) || pv < best$p)) best <- list(j = j, r = r, p = pv) }
      }
      if (is.null(best)) break
      fixed[[as.character(best$j)]] <- best$r; score <- score - log10(max(best$p, 1e-16))
      f <- f[chars(f, best$j) == best$r]; b <- b[chars(b, best$j) == best$r]
      if (length(f) < min_occ) break
    }
    if (!length(fixed)) break
    m <- rep(".", 2 * flank + 1); m[flank + 1] <- tolower(central); for (j in names(fixed)) m[as.integer(j)] <- fixed[[j]]
    pat <- paste(m, collapse = ""); rx <- gsub("\\.", ".", toupper(pat))
    in_f <- grepl(paste0("^", rx, "$"), fg); in_b <- grepl(paste0("^", rx, "$"), bg)
    out[[length(out) + 1]] <- data.frame(motif = sub("^\\.+", "", sub("\\.+$", "", pat)), window = pat, score = round(score, 1),
      fg_matches = sum(in_f), fg_size = length(fg), bg_matches = sum(in_b), bg_size = length(bg),
      fold = round((sum(in_f) / length(fg)) / max(sum(in_b) / length(bg), 1e-9), 2), stringsAsFactors = FALSE)
    fg <- fg[!in_f]; bg <- bg[!in_b]
  }
  if (length(out)) do.call(rbind, out) else data.frame()
}

# The Python that has kinase_library: DPDP_KL_PYTHON, the web app's venv, or the repo's dev venv
kl_python <- function(root) {
  cand <- c(Sys.getenv("DPDP_KL_PYTHON"), "/srv/dpdp/venv/bin/python", file.path(root, ".venv-kl", "bin", "python"))
  for (py in cand[nzchar(cand) & file.exists(cand)])
    if (identical(suppressWarnings(system2(py, c("-c", shQuote("import kinase_library")), stdout = FALSE, stderr = FALSE)), 0L)) return(py)
  NA_character_
}
