# targets pipeline for the Phase 0 airway prototype (plan §6, §14).
#
# Orchestrates: fixed DE step -> ID-mapping matrix -> enrichment battery
# (per branch) -> stability metrics -> fragility lookup table. Requires the
# Bioconductor stack declared in idmapaudit's DESCRIPTION (Suggests); run
# with `targets::tar_make()` from this directory in an environment that has
# them installed (see .github/workflows/R-CMD-check.yaml for the exact set,
# or `renv::restore()` against the repo's renv.lock).
library(targets)

tar_option_set(packages = c("idmapaudit"))

for (f in list.files("R", full.names = TRUE)) source(f)

FDR_CUTOFF <- 0.05
ORA_DBS <- c("GO_BP", "KEGG", "Reactome", "WikiPathways")

list(
  tar_target(de, run_airway_de(fdr_cutoff = FDR_CUTOFF)),

  tar_target(
    mapping_matrix,
    run_mapping_matrix(
      query_genes = de$sig_genes,
      background_genes = de$background_genes,
      resolvers = airway_resolvers(),
      resolution_policies = airway_resolution_policies()
    )
  ),

  # One enrichment battery per branch; kept as a single target (not a
  # tar_target per branch) for Phase 0 simplicity — split into
  # tar_target(..., pattern = map(branch)) once branch count grows with the
  # full three-dataset study (plan §Phase 1) and re-computation cost bites.
  tar_target(enrichment_by_branch, {
    stats::setNames(
      lapply(names(mapping_matrix), function(branch) {
        b <- mapping_matrix[[branch]]
        run_enrichment_battery(
          query = b$query, background = b$background,
          modes = "ORA", dbs = ORA_DBS, fdr_cutoff = FDR_CUTOFF
        )
      }),
      names(mapping_matrix)
    )
  }),

  # Long-format significance table per DB, one row per (pathway, branch),
  # feeding compute_stability_metrics()/build_fragility_table() (methods
  # note §4). A pathway absent from a branch's result table is untestable
  # in that branch, not "not significant" (methods note §1).
  tar_target(sig_tables_by_db, {
    stats::setNames(lapply(ORA_DBS, function(db) {
      key <- paste0("ORA.", db)
      rows <- lapply(names(enrichment_by_branch), function(branch) {
        tbl <- enrichment_by_branch[[branch]][[key]]
        if (is.null(tbl) || nrow(tbl) == 0) {
          return(data.frame(
            pathway_id = character(0), branch = character(0),
            testable = logical(0), significant = logical(0)
          ))
        }
        data.frame(
          pathway_id = tbl$pathway_id, branch = branch,
          testable = TRUE, significant = tbl$significant,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    }), ORA_DBS)
  }),

  tar_target(fragility_table, {
    rows <- lapply(ORA_DBS, function(db) {
      build_fragility_table(
        sig_tables_by_db[[db]], dataset = "airway", db = db, mode = "ORA"
      )
    })
    do.call(rbind, rows)
  }),

  tar_target(pairwise_concordance_by_db, {
    stats::setNames(lapply(ORA_DBS, function(db) {
      key <- paste0("ORA.", db)
      sig_by_branch <- lapply(enrichment_by_branch, function(x) {
        tbl <- x[[key]]
        tbl$pathway_id[tbl$significant]
      })
      pairwise_concordance(sig_by_branch)
    }), ORA_DBS)
  }),

  tar_target(
    fragility_table_tsv,
    write_fragility_table(fragility_table, "results/airway_fragility_table.tsv"),
    format = "file"
  )
)
