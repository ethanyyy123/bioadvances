#' Build the reusable per-pathway fragility lookup table
#'
#' The tool's headline community-facing artifact (plan §6, §9): one row per
#' `(pathway, dataset, db, mode)` with its fragility score, entropy index,
#' and stable/fragile classification, exported flat so it can be inspected
#' directly (CSV/TSV) rather than only read out of a figure.
#'
#' @param sig_table a data.frame with one row per `(pathway_id, branch)`
#'   combination and columns `pathway_id`, `branch`, `testable` (logical),
#'   `significant` (logical). Callers assemble this from
#'   [run_enrichment_battery()] output across all branches for one fixed
#'   `(dataset, db, mode)`.
#' @param dataset,db,mode labels stamped onto every output row, so tables
#'   from multiple runs can be row-bound into the full lookup table.
#' @param tau_lo,tau_hi passed to [classify_fragility()].
#' @return a data.frame with columns `dataset`, `db`, `mode`, `pathway_id`,
#'   `k_p`, `F`, `H`, `V`, `class`.
#' @export
build_fragility_table <- function(sig_table, dataset, db, mode,
                                   tau_lo = 0.1, tau_hi = 0.9) {
  stopifnot(all(c("pathway_id", "branch", "testable", "significant") %in% names(sig_table)))
  split_by_pathway <- split(sig_table, sig_table$pathway_id)
  rows <- lapply(names(split_by_pathway), function(pid) {
    sub <- split_by_pathway[[pid]]
    fs <- fragility_score(sig = sub$significant, testable = sub$testable)
    data.frame(
      dataset = dataset, db = db, mode = mode, pathway_id = pid,
      k_p = fs$k_p, F = fs$F, H = fs$H, V = fs$V,
      class = classify_fragility(fs$F, tau_lo = tau_lo, tau_hi = tau_hi),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Write the fragility lookup table to a standalone TSV
#'
#' Plan §6 repo-hygiene requirement: the fragility lookup table must be
#' exported as a clean, standalone file from the start, not just embedded in
#' a figure.
#'
#' @param fragility_table output of [build_fragility_table()] (optionally
#'   row-bound across multiple dataset/db/mode runs).
#' @param path output file path.
#' @return invisibly, `path`.
#' @export
write_fragility_table <- function(fragility_table, path) {
  utils::write.table(
    fragility_table, file = path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  invisible(path)
}
