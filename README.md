# idmapaudit — a gene-ID-mapping audit tool

Every RNA-seq → pathway enrichment pipeline makes an invisible choice: which
gene identifier namespace to run enrichment in (Ensembl, Entrez, HGNC
symbol), and how to resolve genes that don't map cleanly. This repo holds a
tool that fixes everything else in a differential-expression → enrichment
pipeline, varies only that choice, and reports a per-pathway **fragility
score**: how much of a pathway's "significant" call survives the
ID-mapping decision, versus being an artifact of it.

Prototype dataset: `airway` (GSE52778, Himes et al. 2014), a glucocorticoid
(dexamethasone) response study in airway smooth muscle. Full study design,
research questions, and statistical framework are in the project plan; the
precise metric definitions this code implements are in
[`docs/methods-note-fragility-score.md`](docs/methods-note-fragility-score.md)
— read that first if you're modifying `R/stability-metrics.R`.

## Repository layout

```
idmapaudit/            R package: the reusable audit tool
  R/
    stability-metrics.R    Jaccard/overlap/rank-concordance, fragility score,
                            entropy/variance index, cross-dataset agreement,
                            mapping-noise-null excess-instability — pure R,
                            no Bioconductor dependency, fully unit tested
    mapping-matrix.R        ID-mapping resolution policies (drop/first/list)
                             and the query/background mapping-branch builder
    enrichment-battery.R    ORA/GSEA wrapper (clusterProfiler/ReactomePA/
                             fgsea) and the mapping-noise null generator
    fragility-lookup.R      builds and exports the per-pathway fragility
                             lookup table
  tests/testthat/           unit tests, including hand-computed known-answer
                             cases for every stability metric

analysis/airway-prototype/  the Phase 0 prototype run (plan §14)
  R/resolvers.R              the 5 concrete ID-mapping resolvers (Ensembl
                              unversioned/versioned, Entrez via org.Hs.eg.db,
                              Symbol via org.Hs.eg.db, Symbol via biomaRt)
  R/de-step.R                 the single, fixed DESeq2 run (plan §5.2 —
                              ID-mapping is the only manipulated variable)
  _targets.R                  orchestrates DE -> mapping matrix -> enrichment
                              -> stability metrics -> fragility table

docs/methods-note-fragility-score.md   precise math spec, written before
                                        implementation (plan §14)
data/, results/              raw data (gitignored) vs. generated output
                              (tracked) — plan §6 repo hygiene
```

## Reproducing the Phase 0 prototype

Requires R ≥ 4.2 and the Bioconductor packages listed in
`idmapaudit/DESCRIPTION` under `Suggests` (`airway`, `DESeq2`,
`SummarizedExperiment`, `org.Hs.eg.db`, `AnnotationDbi`, `biomaRt`,
`clusterProfiler`, `ReactomePA`, `fgsea`) plus `targets`. These are heavy
Bioconductor packages and are **not** installed in every environment this
repo is checked out in — see `.github/workflows/R-CMD-check.yaml`'s
`airway-smoke-test` job for the exact install recipe (`r-lib/actions`,
which resolves Bioconductor packages automatically from `DESCRIPTION`).

```r
# from the repo root
BiocManager::install(c(
  "airway", "DESeq2", "SummarizedExperiment", "org.Hs.eg.db",
  "AnnotationDbi", "biomaRt", "clusterProfiler", "ReactomePA", "fgsea"
))
install.packages(c("targets", "testthat", "roxygen2"))

R CMD INSTALL idmapaudit

setwd("analysis/airway-prototype")
targets::tar_make()
```

This writes `analysis/airway-prototype/results/airway_fragility_table.tsv`
— the reusable per-pathway fragility lookup table for the airway dataset,
across GO (BP), KEGG, Reactome, and WikiPathways, under the five primary
ID-mapping branches plus the Entrez/Symbol resolution-policy sub-branches.

`renv.lock` at the repo root is a scaffold (`renv::init(bare = TRUE)`
followed by `renv::snapshot()` in an environment that only had base R,
`renv`, and `roxygen2` installed) — it correctly pins the R version but not
yet the Bioconductor stack, since that stack isn't installed in the
environment this repo was scaffolded in. Regenerate a complete lockfile by
running `renv::snapshot()` again after installing the packages above.

## Running just the unit tests

The stability-metric and mapping-resolution-policy logic
(`R/stability-metrics.R`, `R/mapping-matrix.R`, `R/fragility-lookup.R`) has
no Bioconductor dependency and can be tested standalone:

```r
setwd("idmapaudit")
for (f in list.files("R", full.names = TRUE)) source(f)
testthat::test_dir("tests/testthat")
```

This is what the `unit-tests` CI job runs on every push, deliberately kept
separate from the slower `airway-smoke-test` job (plan §6: CI should catch
breakage in the core logic without every push waiting on a full
Bioconductor install).

## Scope of this Phase 0 scaffold

This commit implements the "immediate next steps" from the project plan:
the package skeleton with CI, the fully unit-tested stability-metric core,
the ID-mapping matrix and enrichment-battery wiring, and the airway
`targets` pipeline. It does **not** yet include: the two additional GEO
datasets (GSE159094, GSE80651) for the cross-dataset generalization study,
the GSEA gene-set loading (the `fgsea` call in
`R/enrichment-battery.R::.run_one_enrichment` is stubbed pending the
MSigDB/Reactome gene-set source decision), the GR/NR3C1 ChEA 2022 ground
truth anchor, or the manuscript. See the project plan for the full 12-week
timeline.
