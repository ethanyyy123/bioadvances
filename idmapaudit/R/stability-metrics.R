#' Jaccard index of two significant-pathway sets
#'
#' @param sig_a,sig_b character vectors of pathway IDs called significant
#'   in two ID-mapping branches.
#' @return numeric in \[0, 1\]. Defined as 1 when both sets are empty.
#' @export
jaccard_index <- function(sig_a, sig_b) {
  sig_a <- unique(sig_a)
  sig_b <- unique(sig_b)
  union_n <- length(union(sig_a, sig_b))
  if (union_n == 0) return(1)
  length(intersect(sig_a, sig_b)) / union_n
}

#' Overlap coefficient (Szymkiewicz-Simpson) of two significant-pathway sets
#'
#' @inheritParams jaccard_index
#' @return numeric in \[0, 1\]. Defined as 1 when both sets are empty, 0 when
#'   exactly one is empty.
#' @export
overlap_coefficient <- function(sig_a, sig_b) {
  sig_a <- unique(sig_a)
  sig_b <- unique(sig_b)
  if (length(sig_a) == 0 && length(sig_b) == 0) return(1)
  if (length(sig_a) == 0 || length(sig_b) == 0) return(0)
  length(intersect(sig_a, sig_b)) / min(length(sig_a), length(sig_b))
}

#' Pairwise set-level concordance across all ID-mapping branches
#'
#' @param sig_by_branch a named list of character vectors: one significant
#'   pathway-ID vector per ID-mapping branch.
#' @return a data.frame with one row per unordered branch pair and columns
#'   `branch_a`, `branch_b`, `jaccard`, `overlap`.
#' @export
pairwise_concordance <- function(sig_by_branch) {
  stopifnot(is.list(sig_by_branch), !is.null(names(sig_by_branch)))
  branches <- names(sig_by_branch)
  pairs <- utils::combn(branches, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    a <- sig_by_branch[[pr[1]]]
    b <- sig_by_branch[[pr[2]]]
    data.frame(
      branch_a = pr[1],
      branch_b = pr[2],
      jaccard = jaccard_index(a, b),
      overlap = overlap_coefficient(a, b),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Rank concordance of an enrichment statistic across two branches
#'
#' Restricted to pathways testable (non-`NA` statistic) in both branches, per
#' the methods note §3 — untestable pathways are excluded, not imputed.
#'
#' @param stat_a,stat_b named numeric vectors (names = pathway IDs) of the
#'   ranking statistic (`-log10(p_adj)` for ORA, `NES` for GSEA) in each
#'   branch. Entries absent from a vector, or `NA`, are treated as untestable.
#' @return a list with `rho` (Spearman correlation) and `n` (number of
#'   pathways used).
#' @export
rank_concordance <- function(stat_a, stat_b) {
  common <- intersect(names(stat_a)[!is.na(stat_a)], names(stat_b)[!is.na(stat_b)])
  if (length(common) < 3) {
    return(list(rho = NA_real_, n = length(common)))
  }
  list(
    rho = stats::cor(stat_a[common], stat_b[common], method = "spearman"),
    n = length(common)
  )
}

#' Audit / fragility score for a single pathway
#'
#' Implements methods note §4: the fraction of testable branches in which a
#' pathway is called significant, plus the entropy- and variance-based
#' fragility indices that treat "always in" and "always out" symmetrically as
#' stable, and mid-range as fragile.
#'
#' @param sig logical vector: one entry per ID-mapping branch, `TRUE` if the
#'   pathway was called significant in that branch.
#' @param testable logical vector, same length as `sig`: `TRUE` if the
#'   pathway had at least one mapped gene in that branch. Branches where
#'   `testable` is `FALSE` are excluded from the denominator entirely (they
#'   are not counted as "not significant"). Defaults to `!is.na(sig)`.
#' @return a list with `F` (fragility proportion, `NA` if no branch is
#'   testable), `H` (binary entropy of `F`), `V` (Bernoulli variance of `F`),
#'   and `k_p` (number of testable branches).
#' @export
fragility_score <- function(sig, testable = !is.na(sig)) {
  stopifnot(length(sig) == length(testable))
  sig <- sig[testable]
  k_p <- length(sig)
  if (k_p == 0) {
    return(list(F = NA_real_, H = NA_real_, V = NA_real_, k_p = 0L))
  }
  F_p <- sum(sig) / k_p
  list(F = F_p, H = binary_entropy(F_p), V = F_p * (1 - F_p), k_p = k_p)
}

#' Binary (Shannon) entropy in bits, with the 0 log 0 := 0 convention
#'
#' @param p numeric vector of probabilities in \[0, 1\].
#' @return numeric vector, same length as `p`, in \[0, 1\].
#' @export
binary_entropy <- function(p) {
  stopifnot(all(p >= 0 & p <= 1, na.rm = TRUE))
  term <- function(x) ifelse(x <= 0, 0, x * log2(x))
  -(term(p) + term(1 - p))
}

#' Classify a pathway as stable or mapping-fragile from its fragility score
#'
#' @param F_p numeric fragility proportion (methods note §4).
#' @param tau_lo,tau_hi lower/upper stability thresholds. Default 0.1 / 0.9
#'   per methods note §4.2.
#' @return character: `"stable"`, `"fragile"`, or `NA_character_` if `F_p` is
#'   `NA`.
#' @export
classify_fragility <- function(F_p, tau_lo = 0.1, tau_hi = 0.9) {
  stopifnot(tau_lo < tau_hi)
  vapply(F_p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x <= tau_lo || x >= tau_hi) "stable" else "fragile"
  }, character(1))
}

#' Gene-level attrition rate for one ID-mapping branch
#'
#' Methods note §6: the fraction of a native-ID gene list that has *no*
#' target under the branch's mapping.
#'
#' This is deliberately defined on native keys with zero targets, not on the
#' size of the image set, because the image-set-ratio version conflates
#' attrition with two things that are not attrition at all: (i) many-to-one
#' collapse (two native IDs mapping to the same target both succeeded; the
#' image is merely smaller) would be reported as attrition though nothing
#' was lost, and (ii) one-to-many expansion under the `"list"`/`"all"`
#' resolution policy (plan §5.1) can inflate the image past the native set
#' size, driving the image-ratio version negative. Counting native keys with
#' at least one target is well-defined and non-negative under every
#' resolution policy.
#'
#' @param native_genes character vector of gene IDs in the native
#'   (pre-mapping) namespace, e.g. the DESeq2-significant Ensembl gene list.
#' @param mapped_native_genes character vector: the subset of `native_genes`
#'   that has at least one target under the branch's mapping (i.e. the
#'   `from`-keys of the resolved mapping table, restricted to `native_genes`
#'   -- not the mapped-to IDs themselves).
#' @return numeric in \[0, 1\].
#' @export
attrition_rate <- function(native_genes, mapped_native_genes) {
  native_genes <- unique(native_genes)
  n <- length(native_genes)
  if (n == 0) return(NA_real_)
  1 - length(intersect(native_genes, mapped_native_genes)) / n
}

#' Cross-dataset agreement of stable/fragile classification
#'
#' Methods note §5: simple agreement rate and Fleiss' kappa across datasets
#' (treated as raters) for a shared set of pathways, each classified
#' `"stable"` or `"fragile"` independently within its own dataset.
#'
#' @param class_matrix a pathway-by-dataset character matrix/data.frame of
#'   `"stable"`/`"fragile"` calls (rows = pathways, columns = datasets). Rows
#'   with any `NA` are dropped before computing agreement.
#' @return a list with `agreement_rate` (fraction of pathways where all
#'   datasets agree) and `fleiss_kappa`. Both are `NA` if fewer than two
#'   pathways remain after dropping incomplete rows, or fewer than two
#'   raters (datasets/columns) are supplied -- Fleiss' kappa is undefined
#'   for a single rater (its `raters * (raters - 1)` denominator is 0), and
#'   "agreement across datasets" is not a meaningful question for one.
#' @export
cross_dataset_agreement <- function(class_matrix) {
  m <- as.matrix(class_matrix)
  complete <- stats::complete.cases(m)
  m <- m[complete, , drop = FALSE]
  n <- nrow(m)
  if (n == 0 || ncol(m) < 2) {
    return(list(agreement_rate = NA_real_, fleiss_kappa = NA_real_))
  }
  agreement_rate <- mean(apply(m, 1, function(r) length(unique(r)) == 1))

  categories <- c("stable", "fragile")
  raters <- ncol(m)
  n_ic <- t(apply(m, 1, function(r) table(factor(r, levels = categories))))
  p_j <- colSums(n_ic) / (n * raters)
  P_i <- rowSums(n_ic * (n_ic - 1)) / (raters * (raters - 1))
  P_bar <- mean(P_i)
  P_e <- sum(p_j^2)
  kappa <- if (isTRUE(all.equal(P_e, 1))) NA_real_ else (P_bar - P_e) / (1 - P_e)

  list(agreement_rate = agreement_rate, fleiss_kappa = kappa)
}

#' Excess instability against a mapping-noise null distribution
#'
#' Methods note §7 step 4: compares an observed stability statistic (e.g. a
#' pairwise Jaccard index, or a pathway's `H_p`) against the distribution of
#' that same statistic under pure random gene dropout at matched attrition
#' rates.
#'
#' @param x_obs numeric, the observed statistic.
#' @param x_null numeric vector, the statistic recomputed over `R`
#'   noise-null replicates.
#' @return a list with `excess` (`x_obs` minus the null median), `z_score`,
#'   and `percentile` (fraction of null replicates at or below `x_obs`).
#' @export
excess_instability <- function(x_obs, x_null) {
  x_null <- x_null[!is.na(x_null)]
  if (length(x_null) == 0) {
    return(list(excess = NA_real_, z_score = NA_real_, percentile = NA_real_))
  }
  sd_null <- stats::sd(x_null)
  list(
    excess = x_obs - stats::median(x_null),
    z_score = if (sd_null == 0) NA_real_ else (x_obs - mean(x_null)) / sd_null,
    percentile = mean(x_null <= x_obs)
  )
}
