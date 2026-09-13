# Methods note: audit/fragility score and mapping-noise null

Internal spec, written before implementation, per plan §14 ("write the precise
mathematical definitions ... as a short internal methods note before coding
them"). This is the contract the functions in `idmapaudit` must satisfy;
`tests/testthat/test-stability-metrics.R` checks it against hand-computed
cases.

## 1. Inputs

Fix a dataset `d`, a pathway database `DB`, and an enrichment mode
`mode ∈ {ORA, GSEA}`. Let `b = 1, …, k` index the ID-mapping branches (the
primary factor: unversioned Ensembl, versioned Ensembl, Entrez/first,
Entrez/list, Symbol via org.Hs.eg.db, Symbol via biomaRt).

For a pathway `p`, running the enrichment battery on branch `b` yields:

- `testable_p,b ∈ {0,1}` — 1 iff at least one gene of `p`'s gene set is
  present in branch `b`'s mapped background (an untestable pathway
  contributes no evidence and must not silently count as "not
  significant").
- `sig_p,b ∈ {0,1}` — 1 iff `testable_p,b = 1` and `p` is called
  significant at the chosen FDR threshold `α` (default 0.05) in branch `b`.
- `stat_p,b` — the ranking statistic for branch `b` (`-log10(p_adj)` for
  ORA, `NES` for GSEA), defined only where `testable_p,b = 1`.

Let `T_p = {b : testable_p,b = 1}` and `k_p = |T_p|`. A pathway with
`k_p = 0` is excluded from all stability metrics and reported separately
as "untestable in every branch" — it is a data point about coverage, not
about fragility.

## 2. Set-level concordance (pairwise, per dataset/DB/mode)

For branches `b, b'` with significant-pathway sets `Sig_b`, `Sig_b'`
(restricted to pathways testable in both):

```
Jaccard(b, b')  = |Sig_b ∩ Sig_b'| / |Sig_b ∪ Sig_b'|      (define as 1 if both sets empty)
Overlap(b, b')  = |Sig_b ∩ Sig_b'| / min(|Sig_b|, |Sig_b'|) (define as 1 if either set empty and both empty; 0 if one empty and the other isn't)
```

## 3. Rank concordance

Spearman ρ of `stat_·,b` vs `stat_·,b'` over `{p : b, b' ∈ T_p}` only —
pathways untestable in either branch are excluded, not imputed.

## 4. Audit / fragility score (primary novel metric)

Per pathway `p`, per dataset `d` (one DB/mode at a time — do not pool
across DB or mode):

```
F_p = (1 / k_p) * Σ_{b ∈ T_p} sig_p,b        F_p ∈ [0, 1]
```

`F_p` near 0 → stably absent. `F_p` near 1 → stably present. Both are
"stable." Mid-range `F_p` is the fragile case.

### 4.1 Entropy version (captures "flips" symmetrically)

```
H_p = -F_p * log2(F_p) - (1 - F_p) * log2(1 - F_p)     (0 * log2(0) := 0)
```

`H_p ∈ [0, 1]`. `H_p = 0` at `F_p ∈ {0, 1}` (fully stable); `H_p = 1` at
`F_p = 0.5` (maximally fragile — as likely in as out). `H_p` is strictly
decreasing in `|F_p - 0.5|`, so it already satisfies the plan's
requirement that a 5/6 pathway (`F=0.833`, `H≈0.650`) reads as *less*
fragile than a 3/6 pathway (`F=0.5`, `H=1.0`), unlike treating "any
mid-range hit" as equally unstable.

Bernoulli variance `V_p = F_p * (1 - F_p)` is a monotone-equivalent,
cheaper alternative reported alongside `H_p` for readers who prefer a
variance framing; it never changes rank order relative to `H_p`.

### 4.2 Classification

Fixed thresholds `τ_lo = 0.1`, `τ_hi = 0.9` (configurable):

- **stable**: `F_p ≤ τ_lo` or `F_p ≥ τ_hi`
- **mapping-fragile**: `τ_lo < F_p < τ_hi`, *and* the gene-attrition
  regression (§6) assigns it a non-trivial attrition-explained share —
  otherwise it is flagged `borderline-significance` instead (see §4.3)
- **dataset-fragile**: computed at the cross-dataset step (§5), not here

### 4.3 Threshold-proximity control

Because ORA/GSEA calls near `α` flip branch-to-branch for reasons that
have nothing to do with ID mapping, every pathway also gets a
threshold-proximity flag: `near_boundary_p,b = 1` if
`|stat_p,b - stat_threshold| < δ` for a small `δ` (default: within a
2-fold p-value band of `α`). Report `F_p`/`H_p` **both** including and
excluding branches flagged `near_boundary`; a pathway whose fragility
vanishes once boundary-flagged branches are dropped is reclassified
`borderline-significance`, not `mapping-fragile` (plan §8, risk register
item 2).

## 5. Cross-dataset generalization (RQ2)

For a pathway `p` present (testable in ≥1 branch) in all three datasets,
let `class_p,1, class_p,2, class_p,3 ∈ {stable, fragile}` be its binary
classification (mapping-fragile and dataset-fragile both count as
"fragile" here) in each dataset. Report:

- simple agreement rate: fraction of pathways where all three classes
  agree
- Fleiss' κ across the 3 datasets (raters) × 2 categories

A pathway is **dataset-fragile** if its per-dataset `F_p` values disagree
enough that `class_p,·` is not unanimous, *after* already being resolved
as stable-or-fragile within each dataset — i.e. dataset-fragility is
measured on the classification, not by re-averaging `F_p` across
datasets (averaging would hide a pathway that is stably-hit in one
tissue and stably-absent in another behind a misleading mid-range mean).

## 6. Gene-level attrition → pathway instability (RQ3)

Per branch `b`, per dataset `d`:

```
Attrition_b = 1 - |map_b(G_sig) ∩ background_b| / |G_sig|
```

where `G_sig` is the native-Ensembl significant-gene list from the fixed
DE step and `map_b(G_sig)` is that list's image under branch `b`'s
mapping. Regress a per-pathway instability outcome (`H_p`, or the range
of `stat_p,·` across branches) against `Attrition_b`, controlling for
pathway gene-set size, separately for query-list attrition and
background-list attrition (plan §5.3 — these are conflated by default
and must be reported separately).

## 7. Mapping-noise null (the one that must be *built*, not assumed)

Purpose: establish how much set-level instability (§2) is expected from
finite-list-size sampling alone, so excess instability can be attributed
to the *structured* (non-random) nature of real mapping failures.

Procedure, per dataset `d`, per branch `b` with observed attrition rate
`π_b` (§6):

1. Let `n = |G_sig|` (native Ensembl significant genes) and
   `m = |background_d|`.
2. For `r = 1, …, R` (default `R = 1000`):
   - Draw `G_sig^(r,b)` by removing `⌊π_b · n⌋` genes selected **uniformly
     at random without replacement** from `G_sig` (independent of which
     genes real mapping actually drops).
   - Draw `background^(r,b)` analogously from `background_d` at the same
     rate `π_b`.
3. Run the identical enrichment battery on each `(G_sig^(r,b),
   background^(r,b))` and recompute Jaccard/`F_p`/`H_p` exactly as in
   §2–§4, giving a null distribution of each statistic per pathway/pair
   over the `R` replicates.
4. **Excess instability** for a real observed statistic `x_obs` (e.g.
   observed pairwise Jaccard, or observed `H_p`):

```
excess = x_obs - median(x_null^(1..R))
z_score = (x_obs - mean(x_null)) / sd(x_null)
percentile = rank of x_obs within {x_null^(1..R)} / R
```

Report `excess`/`percentile`, not just the raw observed statistic — a
pathway's raw `H_p` alone does not distinguish "unstable because ID
mapping structurally drops its genes" from "unstable because any random
~10% gene dropout would flip a list this size." Only the excess is the
mapping-mechanism claim.

Competitive null (separate, simpler): same procedure but drawing
uniformly random gene sets matched in size to each real pathway (no
relation to `G_sig`), used only as a sanity check that the enrichment
tools' base sensitivity/specificity is well behaved — not used for the
excess-instability calculation above.

## 8. Multiple testing

BH-FDR is computed independently within each `(branch, DB, dataset,
mode)` combination, exactly as an end user would run it in isolation.
p-values are never pooled across branches before correction.
