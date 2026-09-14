#' Run the standard ORA/GSEA battery for one ID-mapping branch
#'
#' Thin, swappable wrapper around clusterProfiler/ReactomePA/fgsea (plan
#' §5.3). Requires the relevant Bioconductor packages, which are declared in
#' `Suggests` rather than `Imports` because they are heavy and not needed for
#' the pure stability-metric logic in this package; callers get an explicit
#' error naming the missing package rather than a cryptic failure.
#'
#' @param query character vector of mapped, branch-specific significant gene
#'   IDs.
#' @param background character vector of mapped, branch-specific background
#'   gene IDs (plan §5.3: redefined per branch, not reused from another
#'   branch's namespace).
#' @param ranked_stat named numeric vector (names = gene IDs in this
#'   branch's namespace) of the DESeq2 test statistic, for GSEA mode. May be
#'   `NULL` if `modes` excludes `"GSEA"`.
#' @param dbs character vector, subset of `c("GO_BP", "GO_MF", "GO_CC",
#'   "KEGG", "Reactome", "WikiPathways")`.
#' @param modes character vector, subset of `c("ORA", "GSEA")`.
#' @param org_db the annotation package/object to pass through to
#'   clusterProfiler (e.g. `org.Hs.eg.db::org.Hs.eg.db`), or a string such as
#'   `"org.Hs.eg.db"`.
#' @param key_type the clusterProfiler/OrgDb key type of `query`/`background`
#'   in *this* branch, e.g. `"ENSEMBL"`, `"SYMBOL"`, or `"ENTREZID"` (see
#'   `AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)` for the full set
#'   `enrichGO()` accepts). This is not optional: `clusterProfiler::enrichGO()`
#'   silently assumes `"ENTREZID"` if it is not supplied, which for any
#'   Ensembl- or Symbol-space branch returns an empty or nonsensical result
#'   rather than an error -- exactly the kind of silent per-branch failure
#'   this package exists to catch, so it must not also be present in its own
#'   enrichment call. There is no equivalent for KEGG/Reactome/WikiPathways
#'   (see `dbs` below): those three functions accept only Entrez-family IDs.
#' @param fdr_cutoff FDR threshold used only to label results; full result
#'   tables are always returned unfiltered so callers can recompute at a
#'   different threshold (plan §8 threshold-sensitivity check).
#' @param dbs character vector, subset of `c("GO_BP", "GO_MF", "GO_CC",
#'   "KEGG", "Reactome", "WikiPathways")`. Only `enrichGO` (GO_*) genuinely
#'   supports arbitrary `key_type` via the OrgDb; `clusterProfiler::enrichKEGG()`
#'   accepts only `"kegg"`/`"ncbi-geneid"`/`"ncbi-proteinid"`/`"uniprot"`, and
#'   `ReactomePA::enrichPathway()`/`clusterProfiler::enrichWP()` accept only
#'   Entrez Gene IDs with no key-type argument at all. Requesting KEGG,
#'   Reactome, or WikiPathways with `key_type != "ENTREZID"` therefore errors
#'   explicitly (see `.run_one_enrichment`) rather than silently converting or
#'   returning an empty/wrong result: the primary namespace factor (plan
#'   §5.1) can only be fully crossed against GO; KEGG/Reactome/WikiPathways
#'   fragility is measured only across the Entrez-space branches (native
#'   Entrez, entrez_first, entrez_list).
#' @return a named list keyed `"<mode>.<db>"`, each element a data.frame with
#'   at least columns `pathway_id`, `p_adjust`, `stat` (`-log10(p_adjust)`
#'   for ORA, `NES` for GSEA), and `significant` (logical, at `fdr_cutoff`).
#' @export
run_enrichment_battery <- function(query, background, ranked_stat = NULL,
                                    dbs = c("GO_BP", "KEGG", "Reactome", "WikiPathways"),
                                    modes = c("ORA", "GSEA"),
                                    org_db = "org.Hs.eg.db",
                                    key_type = "ENTREZID",
                                    fdr_cutoff = 0.05) {
  if ("GSEA" %in% modes && is.null(ranked_stat)) {
    stop("`ranked_stat` is required when `modes` includes \"GSEA\".")
  }
  results <- list()
  for (mode in modes) {
    for (db in dbs) {
      key <- paste(mode, db, sep = ".")
      results[[key]] <- .run_one_enrichment(
        mode = mode, db = db, query = query, background = background,
        ranked_stat = ranked_stat, org_db = org_db, key_type = key_type,
        fdr_cutoff = fdr_cutoff
      )
    }
  }
  results
}

.run_one_enrichment <- function(mode, db, query, background, ranked_stat, org_db, key_type, fdr_cutoff) {
  if (mode == "ORA") {
    .require_pkg("clusterProfiler")
    if (db %in% c("KEGG", "Reactome", "WikiPathways") && !identical(key_type, "ENTREZID")) {
      stop(
        db, " enrichment requires Entrez Gene IDs (clusterProfiler::enrichKEGG(), ",
        "ReactomePA::enrichPathway(), and clusterProfiler::enrichWP() have no general ",
        "key-type argument), but this branch supplies key_type = \"", key_type, "\". ",
        "Restrict `dbs` to GO_* for non-Entrez branches, or run this branch only for ",
        "the Entrez-space arm of the mapping matrix.", call. = FALSE
      )
    }
    res <- switch(db,
      GO_BP = ,
      GO_MF = ,
      GO_CC = clusterProfiler::enrichGO(
        gene = query, universe = background, OrgDb = org_db, keyType = key_type,
        ont = sub("GO_", "", db), pAdjustMethod = "BH"
      ),
      KEGG = clusterProfiler::enrichKEGG(
        gene = query, universe = background, keyType = "kegg", pAdjustMethod = "BH"
      ),
      Reactome = {
        .require_pkg("ReactomePA")
        ReactomePA::enrichPathway(gene = query, universe = background, pAdjustMethod = "BH")
      },
      WikiPathways = clusterProfiler::enrichWP(gene = query, universe = background),
      stop("Unknown ORA database: ", db)
    )
    tbl <- as.data.frame(res)
    if (nrow(tbl) == 0) {
      return(.empty_result())
    }
    data.frame(
      pathway_id = tbl$ID,
      p_adjust = tbl$p.adjust,
      stat = -log10(pmax(tbl$p.adjust, .Machine$double.xmin)),
      significant = tbl$p.adjust < fdr_cutoff,
      stringsAsFactors = FALSE
    )
  } else if (mode == "GSEA") {
    .require_pkg("fgsea")
    stop(
      "GSEA requires a gene-set list for `db`, which must be supplied via ",
      "a pathway-database loader (plan section 5.3: MSigDB Hallmark, ",
      "C2:CP:KEGG, Reactome). Wire the concrete `fgsea::fgsea(pathways = ..., ",
      "stats = ranked_stat, ...)` call here once the gene-set source for ",
      "database '", db, "' is finalized; see docs/methods-note-fragility-score.md."
    )
  } else {
    stop("Unknown mode: ", mode)
  }
}

.empty_result <- function() {
  data.frame(
    pathway_id = character(0), p_adjust = numeric(0),
    stat = numeric(0), significant = logical(0), stringsAsFactors = FALSE
  )
}

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required for this enrichment call but is not ",
      "installed. It is declared in Suggests (Bioconductor); install it with ",
      "BiocManager::install(\"", pkg, "\") to run the real enrichment battery. ",
      "Pure stability-metric functions (see R/stability-metrics.R) do not ",
      "require it.",
      call. = FALSE
    )
  }
}

#' Simulate the mapping-noise null for one branch
#'
#' Methods note §7: draws `R` replicates of pure random gene dropout at the
#' branch's *observed* attrition rate (independent of which genes real
#' mapping actually drops), to give a null distribution against which
#' observed instability can be judged "in excess of chance."
#'
#' @param query_genes,background_genes the *native*-namespace gene lists
#'   before mapping (attrition rate is applied directly to these; this
#'   function does not perform ID translation).
#' @param attrition_rate the branch's observed attrition rate (methods note
#'   §6), applied identically to `query_genes` and `background_genes`.
#' @param R number of replicates (default 1000; methods note §7).
#' @param seed optional RNG seed for reproducibility.
#' @return a list of length `R`, each element a list with `query` and
#'   `background` (the noise-perturbed gene lists), ready to pass into
#'   [run_enrichment_battery()].
#' @export
simulate_mapping_noise_null <- function(query_genes, background_genes,
                                         attrition_rate, R = 1000, seed = NULL) {
  stopifnot(attrition_rate >= 0, attrition_rate <= 1)
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    })
    set.seed(seed)
  }
  n_drop_query <- floor(attrition_rate * length(query_genes))
  n_drop_bg <- floor(attrition_rate * length(background_genes))
  lapply(seq_len(R), function(i) {
    list(
      query = setdiff(query_genes, sample(query_genes, n_drop_query)),
      background = setdiff(background_genes, sample(background_genes, n_drop_bg))
    )
  })
}
