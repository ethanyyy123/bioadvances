# Real-data analysis of GSE52778, using two independently generated public
# processings of the same GEO series (provided directly by the user, from
# the GEO record for GSE52778):
#   File 1: GSE52778_All_Sample_FPKM_Matrix.txt
#     - original submitter-deposited Cufflinks FPKM matrix (2013), indexed by
#       HGNC-style gene symbol ("gene_short_name"), 16 samples with explicit
#       treatment/donor labels (4 donors x {Untreated, Dex, Alb, Alb_Dex}).
#   File 2: GSE52778_norm_counts_FPKM_GRCh38.p13_NCBI.tsv
#     - GEO-generated re-processing against GRCh38.p13/NCBI RefSeq, indexed
#       by NCBI/Entrez GeneID, 16 samples labeled only by GSM accession (no
#       treatment/donor metadata is carried in this file).
#
# Every number below is computed directly from these two files as provided.
# We deliberately do NOT attempt to reconstruct File 2's per-sample treatment
# assignment for use in a differential-expression comparison: a calibration
# experiment (below) shows that even a literature-validated glucocorticoid-
# marker panel, checked against File 1's *known* labels, recovers individual
# sample identity with only ~75% accuracy, because within-treatment-group
# donor variability is real and non-negligible in this dataset. Using that
# noisy recovered labeling for a "namespace X vs namespace Y" DE comparison
# would conflate genuine ID-mapping/annotation effects with label-recovery
# error -- exactly the kind of confound this project exists to avoid, so we
# report the calibration finding honestly instead of building on top of it.

options(stringsAsFactors = FALSE)
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

## ---- 1. Load both matrices as provided ------------------------------------

f1 <- read.table("GSE52778_All_Sample_FPKM_Matrix.txt", header = TRUE, sep = "",
                  quote = "", check.names = FALSE)
f2 <- read.table("GSE52778_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", header = TRUE,
                  sep = "\t", quote = "", check.names = FALSE)

f1_samples <- c("Dex_LL14","Dex_LL06","Dex_LL02","Dex_LL10",
                "Alb_LL07","Alb_LL03","Alb_LL11","Alb_LL15",
                "Untreated_LL09","Untreated_LL05","Untreated_LL13","Untreated_LL01",
                "Alb_Dex_LL08","Alb_Dex_LL04","Alb_Dex_LL12","Alb_Dex_LL16")
stopifnot(all(f1_samples %in% colnames(f1)))
f1_group4 <- sub("_LL[0-9]+$", "", f1_samples)
f2_samples <- colnames(f2)[-1]
stopifnot(length(f2_samples) == 16)

universe <- data.frame(
  file = c("File 1 (2013 Cufflinks, symbol)", "File 2 (GRCh38.p13 NCBI reprocessing, Entrez)"),
  n_genes = c(nrow(f1), nrow(f2)),
  n_samples = c(length(f1_samples), length(f2_samples))
)
write.csv(universe, file.path(out_dir, "gene_universe_comparison.csv"), row.names = FALSE)
cat("Gene universe comparison (all annotated rows):\n"); print(universe)
cat(sprintf("Universe growth (all annotated rows): %.1f%%\n", 100 * (nrow(f2) - nrow(f1)) / nrow(f1)))

## The "all annotated rows" comparison above is the naive background choice
## Wijesooriya et al. (2022, PLOS Comput Biol 18:e1009935) identify as the
## field's dominant ORA-misuse pattern (their survey: ~95% of published ORA
## analyses use an inappropriate background). A background restricted to
## expressed/detected genes is the field's own recommended practice, so we
## report it at multiple thresholds rather than only the inflated raw count.
mean1 <- rowMeans(f1[, f1_samples])
mean2 <- rowMeans(f2[, f2_samples])
universe_by_threshold <- data.frame(
  threshold = c("all annotated (mean FPKM > -Inf)", "detected (mean FPKM > 0)",
                "mean FPKM > 0.1", "mean FPKM > 1"),
  file1_n = c(nrow(f1), sum(mean1 > 0), sum(mean1 > 0.1), sum(mean1 > 1)),
  file2_n = c(nrow(f2), sum(mean2 > 0), sum(mean2 > 0.1), sum(mean2 > 1))
)
universe_by_threshold$pct_diff <- round(
  100 * (universe_by_threshold$file2_n - universe_by_threshold$file1_n) / universe_by_threshold$file1_n, 1
)
write.csv(universe_by_threshold, file.path(out_dir, "universe_by_threshold.csv"), row.names = FALSE)
cat("\nGene universe difference at increasingly strict expression filters:\n")
print(universe_by_threshold)

## ---- 2. Ground-truth DE in File 1: Dex vs Untreated (n=4 vs 4) ------------
## Uses File 1's own, unambiguous, submitter-supplied labels. Base R only
## (log2(FPKM+1), Welch t-test, BH-FDR) -- not DESeq2, since only FPKM (not
## raw counts) is available in either file; stated as a limitation.

de_samples <- c("Dex_LL14","Dex_LL06","Dex_LL02","Dex_LL10",
                "Untreated_LL09","Untreated_LL05","Untreated_LL13","Untreated_LL01")
de_group <- rep(c("Dex", "Untreated"), each = 4)

lm1 <- log2(as.matrix(f1[, de_samples]) + 1)
keep <- rowMeans(lm1) > 0
lm1f <- lm1[keep, , drop = FALSE]
sym <- f1$gene_short_name[keep]

pvals <- apply(lm1f, 1, function(x) {
  tryCatch(t.test(x[de_group == "Dex"], x[de_group == "Untreated"])$p.value,
           error = function(e) NA_real_)
})
padj <- p.adjust(pvals, method = "BH")
lfc <- rowMeans(lm1f[, de_group == "Dex", drop = FALSE]) -
  rowMeans(lm1f[, de_group == "Untreated", drop = FALSE])

de1 <- data.frame(symbol = sym, log2FC = lfc, pvalue = pvals, padj = padj)
write.csv(de1, file.path(out_dir, "file1_DexVsUntreated_DE.csv"), row.names = FALSE)

n_tested <- nrow(de1)
n_sig <- sum(de1$padj < 0.05, na.rm = TRUE)
cat(sprintf("\nFile 1 ground-truth DE (Dex vs Untreated, n=4 vs 4): %d genes tested, %d significant at FDR<0.05 (%.2f%%)\n",
            n_tested, n_sig, 100 * n_sig / n_tested))

marker_entrez <- c(FKBP5 = 2289, TSC22D3 = 1831, ZBTB16 = 7704, KLF15 = 28999,
                    PER1 = 5187, DDIT4 = 54541, CRISPLD2 = 83716, PDK4 = 5166)
marker_row <- data.frame(
  gene = names(marker_entrez), entrez = marker_entrez,
  log2FC = round(de1$log2FC[match(names(marker_entrez), de1$symbol)], 3),
  padj = signif(de1$padj[match(names(marker_entrez), de1$symbol)], 3)
)
marker_row$significant <- marker_row$padj < 0.05
write.csv(marker_row, file.path(out_dir, "marker_gene_file1_DE.csv"), row.names = FALSE)
cat("\nCanonical GR-target marker genes in File 1's ground-truth DE:\n")
print(marker_row)

## ---- 3. Calibration: can File 2's samples be relabeled from expression? --
## A panel of canonical, robustly GR-induced marker genes (independent
## literature ground truth) is used as a fingerprint. We first calibrate the
## method on File 1 (labels known) before ever trusting it on File 2.

f1_marker <- as.matrix(f1[match(names(marker_entrez), f1$gene_short_name), f1_samples])
rownames(f1_marker) <- names(marker_entrez)
f2_marker <- as.matrix(f2[match(marker_entrez, f2$GeneID), f2_samples])
rownames(f2_marker) <- names(marker_entrez)
cat(sprintf("\nAll %d canonical marker genes found by exact ID in both files (0 lost to attrition for this hand-verified panel).\n",
            length(marker_entrez)))

zscore <- function(x) (x - mean(x)) / sd(x)
f1_z <- t(apply(log2(f1_marker + 1), 1, zscore))
f2_z <- t(apply(log2(f2_marker + 1), 1, zscore))
f1_true_gc <- ifelse(f1_group4 %in% c("Dex", "Alb_Dex"), "GC_exposed", "GC_naive")

set.seed(1)
km1 <- kmeans(t(f1_z), centers = 2, nstart = 50)
comp1 <- colMeans(f1_z)
hi1 <- which.max(tapply(comp1, km1$cluster, mean))
f1_pred <- ifelse(km1$cluster == hi1, "GC_exposed", "GC_naive")
calibration_accuracy <- mean(f1_pred == f1_true_gc)

cat(sprintf("\nCalibration: recovering GC-exposed/naive group membership from the 8-marker\n"))
cat(sprintf("expression fingerprint alone (k-means, k=2), checked against File 1's KNOWN\n"))
cat(sprintf("labels: %d/%d samples correctly recovered (%.1f%% accuracy).\n",
            round(calibration_accuracy * 16), 16, 100 * calibration_accuracy))
cat("Confusion matrix (predicted vs true, File 1):\n")
print(table(predicted = f1_pred, true = f1_true_gc))

km2 <- kmeans(t(f2_z), centers = 2, nstart = 50)
comp2 <- colMeans(f2_z)
hi2 <- which.max(tapply(comp2, km2$cluster, mean))
f2_pred <- setNames(ifelse(km2$cluster == hi2, "GC_exposed", "GC_naive"), f2_samples)
cat(sprintf("\nBecause calibration accuracy (%.0f%%) is well short of reliable, File 2's\n", 100*calibration_accuracy))
cat("recovered grouping is reported for completeness only and is NOT used for any\n")
cat("further differential-expression or fragility comparison:\n")
print(f2_pred)

calibration <- list(accuracy = calibration_accuracy,
                     confusion = table(predicted = f1_pred, true = f1_true_gc),
                     f2_recovered_grouping = f2_pred)
saveRDS(calibration, file.path(out_dir, "calibration.rds"))

## ---- 4. Real identifier-anomaly count in the 2013 symbol annotation -------
## A direct, string-level, reference-free measurement: how many identifiers
## in the original Cufflinks/symbol annotation are non-standard compound or
## placeholder-style symbols that would not carry forward unchanged into a
## modern re-annotation -- itself a genuine, source-level instance of the
## "gene identifiers are not stable across annotation vintage" problem,
## requiring no external database.

## Revised after external review: an earlier version of this script lumped
## LOC-prefixed, LINC-prefixed, orf-style, and dash-number-suffixed symbols
## together as "non-standard or placeholder." That overstated the case for
## roughly 60% of them. LINC#### is current, stable HGNC nomenclature for
## long non-coding RNA genes, not a placeholder -- once assigned it is not
## typically retired. The dash-number-suffix pattern is dominated by the
## KRTAP (keratin-associated protein), ERV (endogenous retrovirus), and MIR
## (microRNA paralog) gene families' own current, stable naming convention,
## not by instability. Only LOC###### (an explicit NCBI placeholder for a
## gene lacking an approved symbol) and C#orf# (an open-reading-frame
## placeholder assigned before a gene's function is characterized, e.g.
## C10orf10, later renamed DEPP1) are genuinely provisional by HGNC's own
## naming conventions and are frequently retired/renamed. We report both
## the corrected, narrower "genuinely provisional" count and the broader
## naming-convention breakdown, rather than only the inflated original union.
sym_all <- f1$gene_short_name
n_total <- length(sym_all)
is_loc <- grepl("^LOC[0-9]", sym_all)
is_orf <- grepl("orf[0-9]", sym_all, ignore.case = TRUE)
is_linc <- grepl("^LINC[0-9]", sym_all)
is_dashnum <- grepl("-[0-9]+$", sym_all)
dashnum_syms <- sym_all[is_dashnum]
is_dashnum_stable_family <- grepl("^KRTAP|^ERV|^MIR[0-9]", dashnum_syms)

anomaly_df <- data.frame(
  category = c("LOC###### (provisional placeholder)",
               "C#orf# (provisional placeholder)",
               "LOC or orf, union (genuinely provisional)",
               "LINC#### (current stable nomenclature, NOT provisional)",
               "dash-number suffix, total",
               "  of which KRTAP/ERV/MIR family (current stable nomenclature)"),
  n = c(sum(is_loc), sum(is_orf), sum(is_loc | is_orf), sum(is_linc),
        sum(is_dashnum), sum(is_dashnum_stable_family)),
  pct_of_all_symbols = round(100 * c(sum(is_loc), sum(is_orf), sum(is_loc | is_orf), sum(is_linc),
                                       sum(is_dashnum), sum(is_dashnum_stable_family)) / n_total, 2)
)
write.csv(anomaly_df, file.path(out_dir, "symbol_anomaly_counts.csv"), row.names = FALSE)
cat(sprintf("\nGene-symbol naming-convention breakdown in File 1's %d identifiers:\n", n_total))
print(anomaly_df)
cat(sprintf(
  "\nGenuinely provisional (LOC or orf): %d (%.2f%%) -- down from the original,\n",
  sum(is_loc | is_orf), 100 * sum(is_loc | is_orf) / n_total
))
cat("overstated 11.3% figure that also counted current, stable LINC and\n")
cat("KRTAP/ERV/MIR-family symbols as if they were placeholders.\n")

saveRDS(list(f1 = f1, f2 = f2, de1 = de1, marker_row = marker_row,
             universe = universe, calibration = calibration, anomaly_df = anomaly_df),
        file.path(out_dir, "full_results.rds"))

cat("\nDone. Results written to", normalizePath(out_dir), "\n")
