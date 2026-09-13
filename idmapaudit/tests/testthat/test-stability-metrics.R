test_that("jaccard_index matches hand computation", {
  a <- c("a", "b", "c")
  b <- c("b", "c", "d")
  expect_equal(jaccard_index(a, b), 2 / 4)
  expect_equal(jaccard_index(character(0), character(0)), 1)
  expect_equal(jaccard_index(a, character(0)), 0)
  expect_equal(jaccard_index(a, a), 1)
})

test_that("overlap_coefficient matches hand computation", {
  a <- c("a", "b", "c")
  b <- c("b", "c", "d")
  expect_equal(overlap_coefficient(a, b), 2 / 3)
  expect_equal(overlap_coefficient(character(0), character(0)), 1)
  expect_equal(overlap_coefficient(a, character(0)), 0)
})

test_that("pairwise_concordance computes every branch pair correctly", {
  branches <- list(
    br1 = c("p1", "p2", "p3"),
    br2 = c("p2", "p3", "p4"),
    br3 = c("p1", "p2", "p3", "p4")
  )
  out <- pairwise_concordance(branches)
  expect_equal(nrow(out), 3)
  row12 <- out[out$branch_a == "br1" & out$branch_b == "br2", ]
  expect_equal(row12$jaccard, 2 / 4)
  expect_equal(row12$overlap, 2 / 3)
  row13 <- out[out$branch_a == "br1" & out$branch_b == "br3", ]
  expect_equal(row13$jaccard, 3 / 4)
  expect_equal(row13$overlap, 1)
})

test_that("rank_concordance restricts to jointly testable pathways", {
  a <- c(p1 = 1, p2 = 2, p3 = 3, p4 = 4)
  b <- c(p1 = 4, p2 = 3, p3 = 2, p4 = 1)
  out <- rank_concordance(a, b)
  expect_equal(out$rho, -1)
  expect_equal(out$n, 4)

  a2 <- c(p1 = 1, p2 = NA, p3 = 3, p4 = 4)
  b2 <- c(p1 = 4, p3 = 2, p4 = 1)
  out2 <- rank_concordance(a2, b2)
  expect_equal(out2$n, 3)
  expect_equal(out2$rho, -1)

  out3 <- rank_concordance(c(p1 = 1, p2 = 2), c(p1 = 1))
  expect_true(is.na(out3$rho))
  expect_equal(out3$n, 1)
})

test_that("binary_entropy hits known values", {
  expect_equal(binary_entropy(0), 0)
  expect_equal(binary_entropy(1), 0)
  expect_equal(binary_entropy(0.5), 1)
  expect_equal(binary_entropy(5 / 6), 0.6500224, tolerance = 1e-6)
})

test_that("fragility_score reproduces the methods-note worked example", {
  # 3 of 6 branches significant: F = 0.5, maximally fragile (H = 1)
  out <- fragility_score(rep(c(TRUE, FALSE), each = 3))
  expect_equal(out$F, 0.5)
  expect_equal(out$H, 1)
  expect_equal(out$V, 0.25)
  expect_equal(out$k_p, 6)

  # always significant: stable at F = 1, H = 0
  out_all <- fragility_score(rep(TRUE, 6))
  expect_equal(out_all$F, 1)
  expect_equal(out_all$H, 0)

  # untestable branches are excluded from the denominator, not counted
  # against the pathway: only the 1 testable branch (significant) counts
  out_partial <- fragility_score(
    sig = c(TRUE, FALSE),
    testable = c(TRUE, FALSE)
  )
  expect_equal(out_partial$k_p, 1)
  expect_equal(out_partial$F, 1)

  out_none <- fragility_score(c(NA, NA), testable = c(FALSE, FALSE))
  expect_equal(out_none$k_p, 0L)
  expect_true(is.na(out_none$F))
})

test_that("classify_fragility applies default thresholds correctly", {
  expect_equal(classify_fragility(0.05), "stable")
  expect_equal(classify_fragility(0.1), "stable")
  expect_equal(classify_fragility(0.9), "stable")
  expect_equal(classify_fragility(0.95), "stable")
  expect_equal(classify_fragility(0.5), "fragile")
  expect_equal(classify_fragility(0.10001), "fragile")
  expect_true(is.na(classify_fragility(NA_real_)))
  expect_equal(
    classify_fragility(c(0.05, 0.5, 0.95)),
    c("stable", "fragile", "stable")
  )
})

test_that("attrition_rate matches hand computation, including duplicate genes", {
  expect_equal(attrition_rate(c("a", "b", "c", "d"), c("a", "b")), 0.5)
  expect_equal(attrition_rate(c("a", "a", "b"), c("a")), 0.5)
  expect_equal(attrition_rate(c("a", "b"), c("a", "b")), 0)
  expect_equal(attrition_rate(c("a", "b"), character(0)), 1)
})

test_that("cross_dataset_agreement matches the hand-derived Fleiss' kappa", {
  m <- rbind(
    c("stable", "stable", "stable"),
    c("stable", "fragile", "stable"),
    c("fragile", "fragile", "fragile")
  )
  out <- cross_dataset_agreement(m)
  expect_equal(out$agreement_rate, 2 / 3)
  # hand derivation in the accompanying methods note review: kappa = 22/40
  expect_equal(out$fleiss_kappa, 0.55, tolerance = 1e-8)
})

test_that("cross_dataset_agreement drops incomplete rows", {
  m <- rbind(
    c("stable", "stable", "stable"),
    c("stable", NA, "fragile")
  )
  out <- cross_dataset_agreement(m)
  expect_equal(out$agreement_rate, 1)
})

test_that("excess_instability matches hand computation", {
  x_null <- c(1, 2, 3, 4, 5)
  out <- excess_instability(5, x_null)
  expect_equal(out$excess, 2)
  expect_equal(out$z_score, 2 / sd(x_null), tolerance = 1e-8)
  expect_equal(out$percentile, 1)

  out_mid <- excess_instability(3, x_null)
  expect_equal(out_mid$excess, 0)
  expect_equal(out_mid$percentile, 0.6)
})
