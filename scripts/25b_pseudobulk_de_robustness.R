# 25b_pseudobulk_de_robustness.R - plan step 10d: robustness of the pseudobulk DE to (a) dropping the extreme biopsy HA98AC_inf and (b) a fuller covariate model (site + biologic + age + ileocolonic disease location).
# Same inputs and edgeR quasi-likelihood pipeline as 25_pseudobulk_de.R. Output: pseudobulk_de_robustness.csv (unit, gene, logFC/FDR for each variant).
suppressMessages({library(edgeR); library(limma)}); wd <- "pbde"; units <- read.csv(file.path(wd, "pb_units.csv"), stringsAsFactors = FALSE); bm <- read.csv("biopsy_metadata_mrvi.csv", stringsAsFactors = FALSE); bm$ileocolonic <- as.numeric(grepl("^ileocolonic", bm$cd_history)); rownames(bm) <- bm$sample_id
out <- list()
for (i in seq_len(nrow(units))) {
  u <- units$unit[i]; meta <- read.csv(file.path(wd, sprintf("pb_%02d_meta.csv", i)), stringsAsFactors = FALSE); cnt <- as.matrix(read.csv(gzfile(file.path(wd, sprintf("pb_%02d_counts.csv.gz", i))), row.names = 1, check.names = FALSE))
  meta$status <- factor(meta$status, levels = c("non_inflamed", "inflamed")); meta$site <- factor(meta$site, levels = c("TI", "AC")); meta$biologic <- factor(meta$biologic, levels = c("none", "anti_TNF", "anti_IL12_IL23")); meta$agez <- as.numeric(scale(meta$age)); meta$ileo <- bm[meta$sample_id, "ileocolonic"]; rownames(meta) <- meta$sample_id; cnt <- cnt[, meta$sample_id, drop = FALSE]
  if (min(table(meta$status)) < 5) next
  run <- function(m, cn, design) { y <- DGEList(cn); keep <- filterByExpr(y, design); y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y); y <- estimateDisp(y, design, robust = TRUE); f <- glmQLFit(y, design, robust = TRUE); q <- glmQLFTest(f, coef = "statusinflamed"); topTags(q, n = Inf, sort.by = "none")$table }
  base <- run(meta, cnt, model.matrix(~ site + status, meta)); res <- data.frame(unit = u, gene = rownames(base), logFC_base = base$logFC, FDR_base = base$FDR, stringsAsFactors = FALSE)
  k <- meta$sample_id != "HA98AC_inf"; if (!all(k) && min(table(meta$status[k])) >= 5) { m2 <- droplevels(meta[k, ]); t2 <- run(m2, cnt[, k, drop = FALSE], model.matrix(~ site + status, m2)); res$logFC_drop_HA98 <- t2[res$gene, "logFC"]; res$FDR_drop_HA98 <- t2[res$gene, "FDR"] } else { res$logFC_drop_HA98 <- NA; res$FDR_drop_HA98 <- NA }
  d3 <- model.matrix(~ site + biologic + agez + ileo + status, meta); if (qr(d3)$rank == ncol(d3)) { t3 <- run(meta, cnt, d3); res$logFC_full <- t3[res$gene, "logFC"]; res$FDR_full <- t3[res$gene, "FDR"] } else { res$logFC_full <- NA; res$FDR_full <- NA }
  out[[length(out) + 1]] <- res; cat(u, nrow(res), "\n") }
write.csv(do.call(rbind, out), "pseudobulk_de_robustness.csv", row.names = FALSE)
