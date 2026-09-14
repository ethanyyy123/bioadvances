#' Resolve a one-to-many/many-to-one ID mapping under a resolution policy
#'
#' Pure-logic core of the ID-mapping step (plan §5.1 secondary factor):
#' given a raw mapping table (possibly containing duplicate keys, e.g. one
#' Ensembl gene mapping to several Entrez IDs), collapse it under one of
#' three policies. This function takes a plain mapping table rather than
#' calling `org.Hs.eg.db`/`biomaRt` directly, so the resolution-policy logic
#' can be unit-tested without those (large) Bioconductor dependencies
#' installed.
#'
#' @param mapping a data.frame/list with columns `from` and `to` (character),
#'   one row per raw mapping edge. A `from` key with zero rows means it
#'   failed to map at all (attrition, not ambiguity).
#' @param policy one of `"drop"` (remove any `from` key with more than one
#'   distinct `to` value), `"first"` (keep the first `to` value per `from`
#'   key, in input row order), or `"all"` (keep every edge, exploding
#'   one-to-many keys into duplicate rows).
#' @return a data.frame with columns `from`, `to`, deduplicated per policy.
#' @export
resolve_mapping <- function(mapping, policy = c("drop", "first", "all")) {
  policy <- match.arg(policy)
  mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
  stopifnot(all(c("from", "to") %in% names(mapping)))
  mapping <- mapping[!is.na(mapping$to) & mapping$to != "", , drop = FALSE]

  if (policy == "all") {
    return(unique(mapping[, c("from", "to")]))
  }

  split_by_from <- split(mapping$to, mapping$from)
  if (policy == "first") {
    resolved <- vapply(split_by_from, function(x) x[1], character(1))
  } else {
    resolved <- vapply(split_by_from, function(x) {
      u <- unique(x)
      if (length(u) == 1) u else NA_character_
    }, character(1))
  }
  out <- data.frame(from = names(resolved), to = unname(resolved), stringsAsFactors = FALSE)
  out[!is.na(out$to), , drop = FALSE]
}

#' Apply an ID-mapping branch to a query gene list and its background,
#' keeping the two redefinitions separable
#'
#' Plan §5.3: a gene that fails to map is invisible to both the query list
#' and the background, and these two effects must be reported separately
#' rather than conflated. This function returns both mapped sets plus the
#' attrition rate of each, computed independently.
#'
#' @param query_genes,background_genes character vectors of native-namespace
#'   (pre-mapping) gene IDs.
#' @param mapping a data.frame with columns `from`, `to` — the *resolved*
#'   mapping (i.e. already passed through [resolve_mapping()]) covering both
#'   `query_genes` and `background_genes`.
#' @return a list with `query` (mapped query gene IDs), `background` (mapped
#'   background gene IDs), `query_attrition`, `background_attrition`.
#' @export
apply_mapping_branch <- function(query_genes, background_genes, mapping) {
  # Deliberately a row filter, not a named-vector lookup: `mapping$to[key]`
  # with duplicate `from` keys (the "all"/"list" resolution policy, plan
  # §5.1) would silently keep only the first match per key under R's `[`
  # semantics. Filtering rows preserves every edge of a one-to-many mapping.
  mapped_query <- mapping$to[mapping$from %in% query_genes]
  mapped_background <- mapping$to[mapping$from %in% background_genes]
  # attrition_rate() takes the *native* keys that found a target, not the
  # mapped-to values themselves (see its docs for why: an image-size ratio
  # is wrong under many-to-one collapse and under the "list" policy).
  list(
    query = mapped_query,
    background = mapped_background,
    query_attrition = attrition_rate(query_genes, intersect(query_genes, mapping$from)),
    background_attrition = attrition_rate(background_genes, intersect(background_genes, mapping$from))
  )
}

#' Build the full ID-mapping matrix for a DE result
#'
#' Orchestrates the primary factor (namespace) x secondary factor
#' (resolution policy, applied only where the plan calls for it — Entrez and
#' Symbol) described in plan §5.1. Namespace resolvers that require
#' Bioconductor packages (`org.Hs.eg.db`, `biomaRt`) are supplied by the
#' caller as `resolver` functions so this function, and its branch-selection
#' logic, stay testable without those packages installed; see
#' `vignettes/airway-prototype.R` for the concrete resolvers used against the
#' airway dataset.
#'
#' @param query_genes,background_genes character vectors of native
#'   (unversioned Ensembl) gene IDs from the fixed DE step.
#' @param resolvers a named list of resolver functions. Each resolver is
#'   called as `resolver(genes)` (the union of `query_genes` and
#'   `background_genes`) and must return a two-column `from`/`to`
#'   data.frame in that gene's *raw, unresolved* mapping (i.e. it may
#'   contain duplicate `from` keys) if the resolver's namespace has a
#'   resolution policy applied afterward, or an already-resolved table
#'   otherwise (e.g. plain Ensembl version-stripping has no ambiguity to
#'   resolve).
#' @param resolution_policies a named list, subset of `names(resolvers)`,
#'   mapping resolver name to the policy/policies to apply, e.g.
#'   `list(entrez = c("first", "list"), symbol_orgdb = c("first", "list"))`.
#'   A resolver not named here is used as-is (already resolved). `"list"`
#'   here is shorthand for [resolve_mapping()]'s `"all"` policy (kept as a
#'   distinct exploded-rows branch, per plan §5.1 item 4).
#' @return a named list, one entry per resulting branch (resolver name, or
#'   `"<resolver>__<policy>"` when a policy was applied), each the output of
#'   [apply_mapping_branch()].
#' @export
run_mapping_matrix <- function(query_genes, background_genes, resolvers,
                                resolution_policies = list()) {
  all_genes <- union(query_genes, background_genes)
  branches <- list()
  for (name in names(resolvers)) {
    raw <- resolvers[[name]](all_genes)
    policies <- resolution_policies[[name]]
    if (is.null(policies)) {
      branches[[name]] <- apply_mapping_branch(query_genes, background_genes, raw)
    } else {
      for (pol in policies) {
        internal_pol <- if (identical(pol, "list")) "all" else pol
        resolved <- resolve_mapping(raw, policy = internal_pol)
        branches[[paste0(name, "__", pol)]] <-
          apply_mapping_branch(query_genes, background_genes, resolved)
      }
    }
  }
  branches
}
