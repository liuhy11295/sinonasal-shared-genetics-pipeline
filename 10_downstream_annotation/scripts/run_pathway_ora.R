suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(KEGGREST)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

root <- Sys.getenv("ENRICHMENT_ROOT")
if (!nzchar(root)) stop("Set ENRICHMENT_ROOT to the final enrichment analysis directory.")
domains <- c("asthma_lower_airway", "atopic_allergic", "autoimmune_IBD", "ENT_infection")

read_genes <- function(path) {
  if (!file.exists(path)) return(character())
  x <- readLines(path, warn = FALSE)
  x <- unique(trimws(x))
  x[nzchar(x)]
}

bh <- function(p) {
  if (length(p) == 0) return(numeric())
  p.adjust(p, method = "BH")
}

map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[nzchar(symbols)])
  if (!length(symbols)) return(data.frame(SYMBOL=character(), ENTREZID=character()))
  suppressMessages({
    m <- AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
                               columns = c("SYMBOL", "ENTREZID"))
  })
  m <- m[!is.na(m$ENTREZID) & nzchar(m$ENTREZID), , drop=FALSE]
  unique(m)
}

build_go_bp_sets <- function(background_symbols) {
  m <- map_symbols_to_entrez(background_symbols)
  if (!nrow(m)) return(list())
  suppressMessages({
    go <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(m$ENTREZID), keytype = "ENTREZID",
                                columns = c("ENTREZID", "GOALL", "ONTOLOGYALL"))
  })
  go <- go[go$ONTOLOGYALL == "BP" & !is.na(go$GOALL), , drop=FALSE]
  go <- merge(go, m, by = "ENTREZID", all.x = FALSE)
  term_names <- AnnotationDbi::select(GO.db, keys = unique(go$GOALL), keytype = "GOID", columns = c("GOID", "TERM"))
  names_map <- setNames(term_names$TERM, term_names$GOID)
  split(go$SYMBOL, go$GOALL) |>
    lapply(unique) |>
    Filter(function(x) length(x) > 0, x = _) |>
    structure(term_names = names_map)
}

build_reactome_sets <- function(background_symbols) {
  m <- map_symbols_to_entrez(background_symbols)
  if (!nrow(m)) return(list())
  suppressMessages({
    er <- AnnotationDbi::select(reactome.db::reactome.db, keys = unique(m$ENTREZID), keytype = "ENTREZID", columns = c("REACTOMEID", "PATHNAME"))
  })
  er <- er[!is.na(er$REACTOMEID), , drop=FALSE]
  if (!nrow(er)) return(list())
  er <- merge(er, m, by = "ENTREZID", all.x = FALSE)
  names_map <- setNames(er$PATHNAME, er$REACTOMEID)
  split(er$SYMBOL, er$REACTOMEID) |>
    lapply(unique) |>
    Filter(function(x) length(x) > 0, x = _) |>
    structure(term_names = names_map)
}

build_kegg_sets <- function(background_symbols) {
  m <- map_symbols_to_entrez(background_symbols)
  if (!nrow(m)) return(list())
  links <- tryCatch(KEGGREST::keggLink("pathway", "hsa"), error = function(e) character())
  plist <- tryCatch(KEGGREST::keggList("pathway", "hsa"), error = function(e) character())
  if (!length(links)) return(list())
  df <- data.frame(ENTREZID = sub("^hsa:", "", names(links)),
                   term_id = sub("^path:", "", as.character(links)),
                   stringsAsFactors = FALSE)
  df <- merge(df, m, by = "ENTREZID", all.x = FALSE)
  names_map <- setNames(sub(" - Homo sapiens \\(human\\)$", "", as.character(plist)), sub("^path:", "", names(plist)))
  split(df$SYMBOL, df$term_id) |>
    lapply(unique) |>
    Filter(function(x) length(x) > 0, x = _) |>
    structure(term_names = names_map)
}

ora <- function(input_symbols, background_symbols, term_sets, database, analysis_name) {
  input_symbols <- unique(input_symbols[input_symbols %in% background_symbols])
  background_symbols <- unique(background_symbols)
  N <- length(background_symbols)
  n <- length(input_symbols)
  term_names <- attr(term_sets, "term_names")
  rows <- list()
  i <- 1
  for (term_id in names(term_sets)) {
    term_genes <- unique(term_sets[[term_id]])
    term_genes <- term_genes[term_genes %in% background_symbols]
    K <- length(term_genes)
    if (K < 3 || K > N - 1 || n == 0) next
    overlap <- intersect(input_symbols, term_genes)
    k <- length(overlap)
    if (k == 0) next
    mat <- matrix(c(k, n-k, K-k, N-K-n+k), nrow=2, byrow=TRUE)
    ft <- fisher.test(mat, alternative = "greater")
    rows[[i]] <- data.frame(
      analysis_name = analysis_name,
      database = database,
      term_id = term_id,
      term_name = ifelse(!is.null(term_names[[term_id]]) && !is.na(term_names[[term_id]]), term_names[[term_id]], term_id),
      input_gene_count = n,
      background_gene_count = N,
      overlap_count = k,
      term_size_in_background = K,
      gene_ratio = paste0(k, "/", n),
      background_ratio = paste0(K, "/", N),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      FDR = NA_real_,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
    i <- i + 1
  }
  out <- if (length(rows)) bind_rows(rows) else data.frame(
    analysis_name=analysis_name, database=database, term_id=character(), term_name=character(),
    input_gene_count=integer(), background_gene_count=integer(), overlap_count=integer(),
    term_size_in_background=integer(), gene_ratio=character(), background_ratio=character(),
    odds_ratio=numeric(), p_value=numeric(), FDR=numeric(), overlap_genes=character()
  )
  if (nrow(out)) out$FDR <- bh(out$p_value)
  out |> arrange(FDR, p_value, desc(overlap_count))
}

run_set <- function(input_file, background_file, out_prefix, out_dir) {
  input <- read_genes(input_file)
  bg <- read_genes(background_file)
  go <- build_go_bp_sets(bg)
  react <- build_reactome_sets(bg)
  kegg <- build_kegg_sets(bg)
  dbs <- list(GO_BP=go, Reactome=react, KEGG=kegg)
  for (db in names(dbs)) {
    res <- ora(input, bg, dbs[[db]], db, out_prefix)
    write_tsv(res, file.path(out_dir, paste0(out_prefix, "_", db, ".tsv")))
  }
}

global_dir <- file.path(root, "enrichment_results/global")
domain_dir <- file.path(root, "enrichment_results/domain")
dir.create(global_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(domain_dir, recursive=TRUE, showWarnings=FALSE)

run_set(file.path(root, "enrichment_inputs/global/global_primary_123_genes.txt"),
        file.path(root, "enrichment_inputs/global/global_background_AB_genes.txt"),
        "enrichment_global_primary_123_vs_AB3782", global_dir)
run_set(file.path(root, "enrichment_inputs/global/global_strict_genes.txt"),
        file.path(root, "enrichment_inputs/global/global_background_AB_genes.txt"),
        "enrichment_global_strict_vs_AB3782", global_dir)
run_set(file.path(root, "enrichment_inputs/global/global_primary_123_genes.txt"),
        file.path(root, "enrichment_inputs/global/global_background_MAGMA_candidate_genes.txt"),
        "enrichment_global_conditional_123_vs_MAGMA205", global_dir)

for (d in domains) {
  input <- file.path(root, "enrichment_inputs/domain", paste0(d, "_primary_genes.txt"))
  if (length(read_genes(input)) < 10) next
  run_set(input, file.path(root, "enrichment_inputs/domain", paste0(d, "_AB_background_genes.txt")),
          paste0("enrichment_domain_", d, "_vs_domain_AB_background"), domain_dir)
  run_set(input, file.path(root, "enrichment_inputs/global/global_background_AB_genes.txt"),
          paste0("enrichment_domain_", d, "_vs_global_AB3782"), domain_dir)
  run_set(input, file.path(root, "enrichment_inputs/domain", paste0(d, "_MAGMA_candidate_background_genes.txt")),
          paste0("enrichment_domain_", d, "_vs_domain_MAGMA_candidate"), domain_dir)
}

rename_outputs <- function() {
  files <- list.files(global_dir, pattern="\\.tsv$", full.names=TRUE)
  files <- c(files, list.files(domain_dir, pattern="\\.tsv$", full.names=TRUE))
  for (f in files) {
    nf <- sub("_GO_BP.tsv$", "_GO_BP.tsv", f)
    nf <- sub("_Reactome.tsv$", "_Reactome.tsv", nf)
    nf <- sub("_KEGG.tsv$", "_KEGG.tsv", nf)
    if (nf != f) file.rename(f, nf)
  }
}
rename_outputs()

top_terms <- function(dir_path) {
  files <- list.files(dir_path, pattern="\\.tsv$", full.names=TRUE)
  rows <- lapply(files, function(f) {
    x <- read_tsv(f, show_col_types=FALSE)
    if (!nrow(x)) return(NULL)
    x |> arrange(FDR, p_value) |> head(10)
  })
  bind_rows(rows)
}
write_tsv(top_terms(domain_dir), file.path(domain_dir, "enrichment_domain_summary_top_terms.tsv"))
write_tsv(bind_rows(top_terms(global_dir), top_terms(domain_dir)), file.path(domain_dir, "enrichment_global_and_domain_summary_top_terms.tsv"))

plot_dot <- function(path, out_pdf, title) {
  x <- read_tsv(path, show_col_types=FALSE)
  if (!nrow(x)) return(FALSE)
  x <- x |> arrange(FDR, p_value) |> head(15)
  x$term_name <- factor(x$term_name, levels=rev(x$term_name))
  p <- ggplot(x, aes(x=-log10(pmax(FDR, 1e-300)), y=term_name, size=overlap_count, color=odds_ratio)) +
    geom_point() + scale_color_viridis_c(option="C") +
    labs(x="-log10(FDR)", y=NULL, title=title, size="Overlap", color="OR") +
    theme_bw(base_size=9)
  ggsave(out_pdf, p, width=8, height=5)
  TRUE
}

fig_dir <- file.path(root, "figures")
dir.create(fig_dir, showWarnings=FALSE)
warnings <- list()
figs <- list(
  list(file.path(global_dir, "enrichment_global_primary_123_vs_AB3782_GO_BP.tsv"), "dotplot_global_primary_123_vs_AB3782_GO_BP.pdf", "Global primary GO BP"),
  list(file.path(global_dir, "enrichment_global_primary_123_vs_AB3782_Reactome.tsv"), "dotplot_global_primary_123_vs_AB3782_Reactome.pdf", "Global primary Reactome"),
  list(file.path(global_dir, "enrichment_global_primary_123_vs_AB3782_KEGG.tsv"), "dotplot_global_primary_123_vs_AB3782_KEGG.pdf", "Global primary KEGG")
)
for (d in domains) {
  figs[[length(figs)+1]] <- list(file.path(domain_dir, paste0("enrichment_domain_", d, "_vs_domain_AB_background_GO_BP.tsv")),
                                 paste0("dotplot_domain_", d, "_GO_BP.pdf"), paste0("Domain ", d, " GO BP"))
}
for (x in figs) {
  ok <- tryCatch(plot_dot(x[[1]], file.path(fig_dir, x[[2]]), x[[3]]), error=function(e) FALSE)
  if (!ok) warnings[[length(warnings)+1]] <- data.frame(figure=x[[2]], reason="no_terms_or_plot_error")
}
if (length(warnings)) write_tsv(bind_rows(warnings), file.path(fig_dir, "figure_generation_warnings.tsv")) else write_tsv(data.frame(figure=character(), reason=character()), file.path(fig_dir, "figure_generation_warnings.tsv"))

cat("ORA enrichment complete\n")
