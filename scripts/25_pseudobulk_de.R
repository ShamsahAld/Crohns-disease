# 25_pseudobulk_de.R - plan step 8c: pseudobulk differential expression, inflamed vs non-inflamed, per lineage and per state.
# Inputs (written by 25a_prepare_pseudobulk.py): pb_units.csv, and for each unit i: pb_<i>_counts.csv.gz (genes x kept biopsies, raw summed UMI) and pb_<i>_meta.csv (sample_id, status, site, patient, biologic, age, n_cells).
# Biopsies enter a unit only if they contribute >= 20 cells to it. Genes are filtered with edgeR::filterByExpr on the design (status + site).
# Primary  : edgeR quasi-likelihood, TMM, design ~ site + status, robust dispersion.  logFC is log2 (inflamed / non-inflamed), BH FDR within unit.
# Patient-aware sensitivity: limma-voom + duplicateCorrelation(block = patient), same design.
# Therapy/age sensitivity   : edgeR QL with ~ site + biologic + age(z) + status.
suppressMessages({library(edgeR); library(limma)})
args <- commandArgs(trailingOnly = TRUE); wd <- if (length(args) >= 1) args[1] else "."; units <- read.csv(file.path(wd, "pb_units.csv"), stringsAsFactors = FALSE); out <- list(); info <- list()
for (i in seq_len(nrow(units))) {
  u <- units[i, , drop = FALSE]; meta <- read.csv(file.path(wd, sprintf("pb_%02d_meta.csv", i)), stringsAsFactors = FALSE); cnt <- as.matrix(read.csv(gzfile(file.path(wd, sprintf("pb_%02d_counts.csv.gz", i))), row.names = 1, check.names = FALSE))
  meta$status <- factor(meta$status, levels = c("non_inflamed", "inflamed")); meta$site <- factor(meta$site, levels = c("TI", "AC")); meta$biologic <- factor(meta$biologic, levels = c("none", "anti_TNF", "anti_IL12_IL23")); meta$agez <- as.numeric(scale(meta$age)); rownames(meta) <- meta$sample_id; cnt <- cnt[, meta$sample_id, drop = FALSE]
  tab <- table(meta$status, meta$site); n_min <- min(table(meta$status))
  if (n_min < 5 || nlevels(droplevels(meta$site)) < 2) { info[[i]] <- data.frame(unit = u$unit, n_biopsies = nrow(meta), n_inflamed = sum(meta$status == "inflamed"), n_non = sum(meta$status == "non_inflamed"), n_genes = NA, note = "skipped: <5 biopsies in a status group"); next }
  design <- model.matrix(~ site + status, meta); y <- DGEList(cnt); keep <- filterByExpr(y, design); y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE); fit <- glmQLFit(y, design, robust = TRUE); qlf <- glmQLFTest(fit, coef = "statusinflamed"); tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  res <- data.frame(unit = u$unit, gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM, F = tt$F, PValue = tt$PValue, FDR = tt$FDR, stringsAsFactors = FALSE)
  # patient-aware: voom + duplicateCorrelation
  cons <- 0; if (any(duplicated(meta$patient))) { v0 <- voom(y, design); dc0 <- duplicateCorrelation(v0, design, block = meta$patient); if (is.finite(dc0$consensus)) cons <- dc0$consensus }; dc <- list(consensus = cons); v <- voom(y, design, block = meta$patient, correlation = cons); fit2 <- lmFit(v, design, block = meta$patient, correlation = cons); fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE); t2 <- topTable(fit2, coef = "statusinflamed", number = Inf, sort.by = "none")
  res$logFC_patient <- t2[res$gene, "logFC"]; res$P_patient <- t2[res$gene, "P.Value"]; res$FDR_patient <- t2[res$gene, "adj.P.Val"]
  # therapy + age adjusted
  des3 <- model.matrix(~ site + biologic + agez + status, meta); ok3 <- qr(des3)$rank == ncol(des3)
  if (ok3) { y3 <- estimateDisp(y, des3, robust = TRUE); f3 <- glmQLFit(y3, des3, robust = TRUE); q3 <- glmQLFTest(f3, coef = "statusinflamed"); t3 <- topTags(q3, n = Inf, sort.by = "none")$table; res$logFC_adj <- t3[res$gene, "logFC"]; res$FDR_adj <- t3[res$gene, "FDR"] } else { res$logFC_adj <- NA; res$FDR_adj <- NA }
  out[[i]] <- res; info[[i]] <- data.frame(unit = u$unit, n_biopsies = nrow(meta), n_inflamed = sum(meta$status == "inflamed"), n_non = sum(meta$status == "non_inflamed"), n_genes = nrow(res), note = sprintf("dup.cor=%.3f; fdr<0.05: %d; fdr_patient<0.05: %d; therapy-adjusted full rank: %s", dc$consensus, sum(res$FDR < 0.05), sum(res$FDR_patient < 0.05, na.rm = TRUE), ok3))
  cat(sprintf("%s: %d biopsies (%d infl / %d non), %d genes, FDR<0.05: %d\n", u$unit, nrow(meta), sum(meta$status == "inflamed"), sum(meta$status == "non_inflamed"), nrow(res), sum(res$FDR < 0.05)))
}
write.csv(do.call(rbind, out), file.path(wd, "pseudobulk_de_results.csv"), row.names = FALSE); write.csv(do.call(rbind, info), file.path(wd, "pseudobulk_de_unit_info.csv"), row.names = FALSE)
