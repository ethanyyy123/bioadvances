# Concrete ID-mapping resolvers for the airway (GSE52778) prototype.
#
# Each resolver is a function(genes) -> data.frame(from, to) consumed by
# idmapaudit::run_mapping_matrix(); see plan §5.1 for the six namespace
# branches and docs/methods-note-fragility-score.md for the definitions
# these feed. All of them require Bioconductor packages that are declared
# in idmapaudit's Suggests and are NOT installed in the lightweight sandbox
# this repo was scaffolded in — they are exercised by the CI smoke test
# (.github/workflows/R-CMD-check.yaml), which installs the full
# Bioconductor stack, not in a local `Rscript` without that stack.

#' Ensembl gene ID, unversioned — the native namespace, identity mapping
#' @export
resolver_ensembl_unversioned <- function(genes) {
  data.frame(from = genes, to = genes, stringsAsFactors = FALSE)
}

#' Ensembl gene ID, versioned (e.g. ENSG00000141510.16)
#'
#' Deliberately included as a built-in positive control (plan §5.1): failing
#' to strip version suffixes is a common, silent real-world bug, and this
#' branch should show dramatic, almost-trivial instability against every
#' pathway database that indexes on the unversioned ID.
#' @export
resolver_ensembl_versioned <- function(genes) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt is required for resolver_ensembl_versioned(); see idmapaudit Suggests.")
  }
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  bm <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "version"),
    filters = "ensembl_gene_id", values = genes, mart = mart
  )
  data.frame(
    from = bm$ensembl_gene_id,
    to = paste0(bm$ensembl_gene_id, ".", bm$version),
    stringsAsFactors = FALSE
  )
}

#' Entrez Gene ID via org.Hs.eg.db — raw (possibly one-to-many) mapping
#'
#' Resolution policy ("first" vs "list") is applied afterward by
#' [idmapaudit::run_mapping_matrix()], not here — this returns every edge.
#' @export
resolver_entrez_orgdb <- function(genes) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("AnnotationDbi and org.Hs.eg.db are required for resolver_entrez_orgdb().")
  }
  tbl <- suppressWarnings(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", columns = "ENTREZID"
  ))
  data.frame(from = tbl$ENSEMBL, to = tbl$ENTREZID, stringsAsFactors = FALSE)
}

#' HGNC/official symbol via org.Hs.eg.db — raw (possibly one-to-many) mapping
#' @export
resolver_symbol_orgdb <- function(genes) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("AnnotationDbi and org.Hs.eg.db are required for resolver_symbol_orgdb().")
  }
  tbl <- suppressWarnings(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", columns = "SYMBOL"
  ))
  data.frame(from = tbl$ENSEMBL, to = tbl$SYMBOL, stringsAsFactors = FALSE)
}

#' HGNC/official symbol via biomaRt
#'
#' Deliberately a second resolver for the "same" namespace as
#' [resolver_symbol_orgdb()] (plan §5.1 item 6): the two can legitimately
#' disagree, which is itself a reportable finding, not a bug to reconcile.
#' Collapsed to first-hit internally, since the plan explicitly keeps the
#' resolution-policy axis a small separate sub-experiment restricted to
#' Entrez/org.Hs.eg.db-symbol rather than a second full factorial.
#' @export
resolver_symbol_biomart <- function(genes) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt is required for resolver_symbol_biomart().")
  }
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  bm <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id", values = genes, mart = mart
  )
  raw <- data.frame(from = bm$ensembl_gene_id, to = bm$hgnc_symbol, stringsAsFactors = FALSE)
  idmapaudit::resolve_mapping(raw, policy = "first")
}

#' All six primary-factor resolvers, named as in plan §5.1
#' @export
airway_resolvers <- function() {
  list(
    ens_unversioned = resolver_ensembl_unversioned,
    ens_versioned = resolver_ensembl_versioned,
    entrez_orgdb = resolver_entrez_orgdb,
    symbol_orgdb = resolver_symbol_orgdb,
    symbol_biomart = resolver_symbol_biomart
  )
}

#' The secondary-factor resolution policies, applied only where plan §5.1
#' says they bite (Entrez and org.Hs.eg.db-symbol)
#' @export
airway_resolution_policies <- function() {
  list(
    entrez_orgdb = c("first", "list"),
    symbol_orgdb = c("first", "list")
  )
}
