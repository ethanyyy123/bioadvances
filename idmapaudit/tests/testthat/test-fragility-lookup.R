test_that("build_fragility_table reproduces per-pathway fragility_score", {
  sig_table <- data.frame(
    pathway_id = rep(c("P1", "P2"), each = 3),
    branch = rep(c("b1", "b2", "b3"), 2),
    testable = TRUE,
    significant = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  out <- build_fragility_table(sig_table, dataset = "airway", db = "GO_BP", mode = "ORA")
  out <- out[order(out$pathway_id), ]

  expect_equal(out$F, c(2 / 3, 0))
  expect_equal(out$k_p, c(3L, 3L))
  expect_equal(out$class, c("fragile", "stable"))
  expect_equal(out$dataset, c("airway", "airway"))
  expect_equal(out$H[out$pathway_id == "P1"], binary_entropy(2 / 3))
})

test_that("build_fragility_table respects per-pathway testability", {
  sig_table <- data.frame(
    pathway_id = c("P1", "P1", "P1"),
    branch = c("b1", "b2", "b3"),
    testable = c(TRUE, FALSE, TRUE),
    significant = c(TRUE, NA, TRUE),
    stringsAsFactors = FALSE
  )
  out <- build_fragility_table(sig_table, dataset = "d1", db = "KEGG", mode = "ORA")
  expect_equal(out$k_p, 2L)
  expect_equal(out$F, 1)
  expect_equal(out$class, "stable")
})

test_that("write_fragility_table round-trips through a TSV file", {
  tbl <- build_fragility_table(
    data.frame(
      pathway_id = c("P1", "P1"), branch = c("b1", "b2"),
      testable = TRUE, significant = c(TRUE, FALSE), stringsAsFactors = FALSE
    ),
    dataset = "d1", db = "GO_BP", mode = "ORA"
  )
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp))
  write_fragility_table(tbl, tmp)
  reloaded <- utils::read.table(tmp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(reloaded$pathway_id, tbl$pathway_id)
  expect_equal(reloaded$F, tbl$F)
})
