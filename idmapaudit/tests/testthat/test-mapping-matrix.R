test_that("resolve_mapping 'drop' removes ambiguous keys only", {
  mapping <- data.frame(
    from = c("g1", "g1", "g2", "g3"),
    to = c("e1", "e2", "e3", "e4"),
    stringsAsFactors = FALSE
  )
  out <- resolve_mapping(mapping, "drop")
  expect_setequal(out$from, c("g2", "g3"))
  expect_equal(out$to[out$from == "g2"], "e3")
  expect_equal(out$to[out$from == "g3"], "e4")
})

test_that("resolve_mapping 'first' keeps the first-seen value per key", {
  mapping <- data.frame(
    from = c("g1", "g1", "g2", "g3"),
    to = c("e1", "e2", "e3", "e4"),
    stringsAsFactors = FALSE
  )
  out <- resolve_mapping(mapping, "first")
  expect_equal(nrow(out), 3)
  expect_equal(out$to[out$from == "g1"], "e1")
  expect_equal(out$to[out$from == "g2"], "e3")
})

test_that("resolve_mapping 'all' explodes one-to-many keys and dedups exact duplicates", {
  mapping <- data.frame(
    from = c("g1", "g1", "g2", "g3", "g3"),
    to = c("e1", "e2", "e3", "e4", "e4"),
    stringsAsFactors = FALSE
  )
  out <- resolve_mapping(mapping, "all")
  expect_equal(nrow(out), 4)
  expect_setequal(out$to[out$from == "g1"], c("e1", "e2"))
})

test_that("resolve_mapping drops rows with missing/blank targets before resolving", {
  mapping <- data.frame(from = c("g1", "g2"), to = c(NA, ""), stringsAsFactors = FALSE)
  out <- resolve_mapping(mapping, "first")
  expect_equal(nrow(out), 0)
})

test_that("apply_mapping_branch computes mapped sets and per-side attrition independently", {
  mapping <- data.frame(
    from = c("g1", "g2", "g3"),
    to = c("e1", "e2", "e3"),
    stringsAsFactors = FALSE
  )
  out <- apply_mapping_branch(
    query_genes = c("g1", "g2", "gX"),
    background_genes = c("g1", "g2", "g3", "gY"),
    mapping = mapping
  )
  expect_setequal(out$query, c("e1", "e2"))
  expect_equal(out$query_attrition, 1 / 3)
  expect_setequal(out$background, c("e1", "e2", "e3"))
  expect_equal(out$background_attrition, 0.25)
})

test_that("apply_mapping_branch keeps every edge of a one-to-many mapping (regression)", {
  # Regression test: an earlier implementation used named-vector indexing
  # (`mapping$to[mapping$from %in% ...]` vs `key[genes]`), which silently
  # drops all but the first match per duplicated `from` key under R's `[`
  # semantics. This must not happen for the "list"/"all" resolution policy.
  mapping <- data.frame(
    from = c("g1", "g1", "g2"),
    to = c("e1", "e2", "e3"),
    stringsAsFactors = FALSE
  )
  out <- apply_mapping_branch(c("g1", "g2"), c("g1", "g2"), mapping)
  expect_setequal(out$query, c("e1", "e2", "e3"))
  expect_length(out$query, 3)
})

test_that("run_mapping_matrix builds one branch per policy and leaves unlisted resolvers as-is", {
  resolvers <- list(
    ens = function(genes) data.frame(from = genes, to = genes, stringsAsFactors = FALSE),
    entrez = function(genes) {
      data.frame(
        from = c("g1", "g1", "g2", "g3"),
        to = c("e1", "e2", "e3", "e4"),
        stringsAsFactors = FALSE
      )
    }
  )
  out <- run_mapping_matrix(
    query_genes = c("g1", "g2"),
    background_genes = c("g1", "g2", "g3"),
    resolvers = resolvers,
    resolution_policies = list(entrez = c("first", "list"))
  )
  expect_named(out, c("ens", "entrez__first", "entrez__list"))
  expect_setequal(out$ens$query, c("g1", "g2"))
  expect_equal(out$ens$query_attrition, 0)
  expect_setequal(out$entrez__first$query, c("e1", "e3"))
  expect_setequal(out$entrez__list$query, c("e1", "e2", "e3"))
})
