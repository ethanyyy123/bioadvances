# Fixed differential-expression step (plan §2, §5.2): run exactly once, in
# native Ensembl-ID space, before any ID-mapping branching happens. This is
# the sole input every mapping branch downstream must share, so that
# ID-mapping choice is the only manipulated variable.

#' Load the airway prototype dataset and run the standard DESeq2 contrast
#'
#' @return a list with `dds` (the fitted `DESeqDataSet`), `res` (the
#'   `DESeqResults` for the default Dex-vs-untreated contrast, standard
#'   apeglm/normal shrinkage), `sig_genes` (Ensembl IDs at
#'   `padj < fdr_cutoff`), `background_genes` (Ensembl IDs passing the
#'   expression filter), and `ranked_stat` (named numeric vector of the
#'   Wald statistic, for GSEA ranking).
#' @param fdr_cutoff significance threshold for the DE step itself (default
#'   0.05). This is independent of, and upstream of, the enrichment FDR
#'   cutoff used later in the pipeline.
#' @export
run_airway_de <- function(fdr_cutoff = 0.05) {
  for (pkg in c("airway", "DESeq2", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for run_airway_de(); see idmapaudit Suggests.")
    }
  }
  utils::data("airway", package = "airway", envir = environment())
  se <- get("airway", envir = environment())

  # native, unversioned Ensembl gene IDs, as documented in plan §4.1
  ens_versioned <- rownames(se)
  rownames(se) <- sub("\\.[0-9]+$", "", ens_versioned)

  SummarizedExperiment::colData(se)$dex <- relevel(
    SummarizedExperiment::colData(se)$dex, ref = "untrt"
  )
  dds <- DESeq2::DESeqDataSet(se, design = ~cell + dex)

  # standard count-threshold filter, filterByExpr-equivalent (plan §5.2)
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]

  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("dex", "trt", "untrt"), alpha = fdr_cutoff)

  sig_genes <- rownames(res)[!is.na(res$padj) & res$padj < fdr_cutoff]
  background_genes <- rownames(res)

  ranked_stat <- stats::setNames(res$stat, rownames(res))
  ranked_stat <- ranked_stat[!is.na(ranked_stat)]

  list(
    dds = dds, res = res,
    sig_genes = sig_genes, background_genes = background_genes,
    ranked_stat = ranked_stat
  )
}
