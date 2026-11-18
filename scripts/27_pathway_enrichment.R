# 27_pathway_enrichment.R - plan step 9: pathway enrichment of the inflamed vs non-inflamed changes per cell unit.
# Input : de_pseudobulk_mrvi_all.csv.gz (26_de_calls.py)   Gene sets: MSigDB Hallmark (H) and Reactome (C2:CP:REACTOME) from the msigdbr package (offline, human symbols).
# (1) ORA on CALLED genes (both methods agree), up and down separately. Background = genes tested in the unit AND among the 3,000 MrVI model genes (only those can be called).
# (2) ORA on the BROAD robust list = robust calls + pseudobulk-only robust genes outside the MrVI gene set; background = all genes tested in the unit.
# (3) Ranked-list GSEA (fgsea multilevel) on signed sqrt(F) of all tested genes of the unit, excluding guarded (Ig, epithelial-dominant, mito, ribo, Hb, sex) and ambient-risk genes.
# Gene sets are restricted to the background and must contain 10-500 genes there; ORA needs >= 3 overlapping genes. BH within (unit, list, collection).
suppressMessages({library(data.table); library(fgsea); library(msigdbr)})
de <- fread(cmd = "gzip -dc de_pseudobulk_mrvi_all.csv.gz"); mv <- as.character(packageVersion("msigdbr"))
get_sets <- function(cat, sub = NULL) { d <- tryCatch(if (is.null(sub)) msigdbr(species = "Homo sapiens", category = cat) else msigdbr(species = "Homo sapiens", category = cat, subcategory = sub), error = function(e) NULL)
  if (is.null(d)) d <- if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = cat) else msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub); split(d$gene_symbol, d$gs_name) }
sets <- list(Hallmark = get_sets("H"), Reactome = get_sets("C2", "CP:REACTOME")); cat("msigdbr", mv, "| Hallmark sets:", length(sets$Hallmark), "| Reactome sets:", length(sets$Reactome), "\n")
ora <- function(genes, bg, ps, minsz = 10, maxsz = 500, minov = 3) {
  genes <- intersect(genes, bg); N <- length(bg); K <- length(genes); out <- list()
  for (nm in names(ps)) { s <- intersect(unique(ps[[nm]]), bg); m <- length(s); if (m < minsz || m > maxsz) next; ov <- intersect(genes, s); k <- length(ov); if (k < minov) next
    p <- phyper(k - 1, m, N - m, K, lower.tail = FALSE); out[[length(out) + 1]] <- data.table(pathway = nm, overlap = k, set_size = m, n_genes = K, bg_size = N, fold_enrichment = (k / K) / (m / N), p = p, genes = paste(sort(ov), collapse = ";")) }
  if (!length(out)) return(NULL); o <- rbindlist(out); o[, padj := p.adjust(p, "BH")]; o[order(p)] }
res_ora <- list(); res_gsea <- list()
for (u in unique(de$unit)) {
  d <- de[unit == u]; tested <- d$gene; bg_m <- d[in_mrvi_genes == TRUE]$gene
  lists <- list(called_up = list(d[called == TRUE & logFC > 0]$gene, bg_m), called_down = list(d[called == TRUE & logFC < 0]$gene, bg_m),
                broad_up = list(d[(robust == TRUE | pb_robust_not_mrvi_covered == TRUE) & logFC > 0]$gene, tested), broad_down = list(d[(robust == TRUE | pb_robust_not_mrvi_covered == TRUE) & logFC < 0]$gene, tested))
  for (ln in names(lists)) { g <- lists[[ln]][[1]]; bg <- lists[[ln]][[2]]; if (length(g) < 8) next
    for (cl in names(sets)) { o <- ora(g, bg, sets[[cl]]); if (!is.null(o)) { o[, `:=`(unit = u, list = ln, collection = cl)]; res_ora[[length(res_ora) + 1]] <- o } } }
  k <- d[guarded == FALSE & ambient_risk == FALSE]; st <- sign(k$logFC) * sqrt(pmax(k$F, 0)); names(st) <- k$gene; st <- st[is.finite(st)]; st <- st[!duplicated(names(st))]
  if (length(st) >= 1000) for (cl in names(sets)) { set.seed(1); fg <- suppressWarnings(fgseaMultilevel(sets[[cl]], st, minSize = 10, maxSize = 500, eps = 0)); fg <- as.data.table(fg); fg[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x, 12), collapse = ";"))]; fg[, `:=`(unit = u, collection = cl, n_ranked = length(st))]; fg[, padj_unit := p.adjust(pval, "BH")]; res_gsea[[length(res_gsea) + 1]] <- fg }
  cat(u, "done\n") }
# (4) ranked-list GSEA on MrVI mean log fold change (status effect; composition-free, all 3,000 model genes, no pseudobulk information); statistics have no biopsy-level error model, so use only as a consistency check
rk <- fread("mrvi_rank_stats.csv"); res_m <- list()
for (u in unique(rk$unit)) { k <- rk[unit == u]; st <- k$stat; names(st) <- k$gene; st <- st[is.finite(st)]; st <- st[!duplicated(names(st))]; if (length(st) < 500) next
  for (cl in names(sets)) { set.seed(1); fg <- suppressWarnings(fgseaMultilevel(sets[[cl]], st, minSize = 10, maxSize = 500, eps = 0)); fg <- as.data.table(fg); fg[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x, 12), collapse = ";"))]; fg[, `:=`(unit = u, collection = cl, n_ranked = length(st))]; fg[, padj_unit := p.adjust(pval, "BH")]; res_m[[length(res_m) + 1]] <- fg } }
fwrite(rbindlist(res_m, fill = TRUE), "pathway_gsea_mrvi_ranked.csv")
fwrite(rbindlist(res_ora, fill = TRUE), "pathway_ora.csv"); fwrite(rbindlist(res_gsea, fill = TRUE), "pathway_gsea.csv"); cat("saved\n")
