test_that("simulate_mapping_noise_null drops the expected count and is reproducible with a seed", {
  query <- letters[1:10]
  background <- LETTERS[1:20]
  reps <- simulate_mapping_noise_null(query, background, attrition_rate = 0.3, R = 5, seed = 42)

  expect_length(reps, 5)
  for (r in reps) {
    expect_equal(length(r$query), 10 - floor(0.3 * 10))
    expect_equal(length(r$background), 20 - floor(0.3 * 20))
    expect_true(all(r$query %in% query))
    expect_true(all(r$background %in% background))
  }

  reps_again <- simulate_mapping_noise_null(query, background, attrition_rate = 0.3, R = 5, seed = 42)
  expect_identical(reps, reps_again)
})

test_that("simulate_mapping_noise_null handles zero attrition (no dropout)", {
  query <- letters[1:5]
  reps <- simulate_mapping_noise_null(query, query, attrition_rate = 0, R = 3, seed = 1)
  for (r in reps) {
    expect_setequal(r$query, query)
  }
})

test_that("simulate_mapping_noise_null rejects an out-of-range attrition rate", {
  expect_error(simulate_mapping_noise_null(letters[1:3], letters[1:3], attrition_rate = 1.5))
  expect_error(simulate_mapping_noise_null(letters[1:3], letters[1:3], attrition_rate = -0.1))
})

test_that("run_enrichment_battery fails fast and clearly when clusterProfiler is unavailable", {
  skip_if(requireNamespace("clusterProfiler", quietly = TRUE),
          "clusterProfiler is installed; the missing-package path is not exercised")
  expect_error(
    run_enrichment_battery(query = c("g1"), background = c("g1", "g2"), modes = "ORA", dbs = "GO_BP"),
    "clusterProfiler"
  )
})

test_that("run_enrichment_battery requires ranked_stat when GSEA is requested", {
  expect_error(
    run_enrichment_battery(query = c("g1"), background = c("g1", "g2"), modes = "GSEA"),
    "ranked_stat"
  )
})

test_that("run_enrichment_battery errors on KEGG/Reactome/WikiPathways for a non-Entrez key_type", {
  skip_if_not(requireNamespace("clusterProfiler", quietly = TRUE))
  expect_error(
    run_enrichment_battery(query = c("g1"), background = c("g1", "g2"), modes = "ORA",
                            dbs = "KEGG", key_type = "ENSEMBL"),
    "Entrez"
  )
})

test_that("run_enrichment_battery defaults organism to \"Homo sapiens\" for enrichWP (regression)", {
  # This is the exact condition that crashed the real airway targets
  # pipeline in CI *after* the .is_total_mapping_failure fix: DE and the
  # ID-mapping matrix completed, enrichment_by_branch reached WikiPathways,
  # and clusterProfiler::enrichWP() -- unlike enrichKEGG() ("hsa") and
  # enrichPathway() ("human"), which both default their organism -- has no
  # default at all, and errored with `argument "organism" is missing, with
  # no default`. Verified here at the argument-default level (no
  # Bioconductor call), since exercising the real error requires the full
  # WikiPathways data download that only the CI smoke test performs.
  expect_equal(formals(run_enrichment_battery)$organism, "Homo sapiens")
})

test_that(".is_total_mapping_failure recognizes AnnotationDbi's zero-valid-keys error (regression)", {
  # This is the exact condition that crashed the real airway targets
  # pipeline in CI: enrichGO()'s internal AnnotationDbi::select() call
  # throws this hard error -- rather than returning NAs with a warning --
  # when *none* of the query genes match the branch's keyType, which is
  # exactly what the versioned-Ensembl positive-control branch produces.
  e <- tryCatch(stop("None of the keys entered are valid keys for 'ENSEMBL'."),
                error = function(e) e)
  expect_true(.is_total_mapping_failure(e))

  other <- tryCatch(stop("some unrelated failure"), error = function(e) e)
  expect_false(.is_total_mapping_failure(other))
})
