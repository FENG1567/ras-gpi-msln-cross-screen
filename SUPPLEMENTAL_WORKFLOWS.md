# Scope of the supporting analyses

The figure-generation and summary-check commands run without downloading external inputs. The primary AsPC-1 and MSLN analysis and the alternative AsPC-1/UM53 coverage and ranking analysis use the public files in [SOURCES.md](SOURCES.md). The included ranking-bias simulation needs no inputs.

## DepMap 24Q2 baseline analysis

Obtain `CRISPRGeneEffect.csv` and `Model.csv` from the same DepMap Public 24Q2 release. The first command retains model rows with `Pancrea` in `OncotreePrimaryDisease`, preserving every gene column. The second computes model medians per gene and the across-gene distribution for the GPI and comparator sets, with a PAAD-only sensitivity.

```bash
python code/ancillary/project_depmap24q2_pancreas.py inputs/CRISPRGeneEffect.csv inputs/Model.csv generated/pancreas_rows.csv.gz
python code/ancillary/depmap24q2_sensitivity.py generated/pancreas_rows.csv.gz inputs/Model.csv generated/depmap
```

Expected manuscript dimensions are 47 pancreas-annotated models, 45 PAAD models and 18,443 gene columns. The official matrix is external. The three manuscript baseline group distributions are available in `data/derived/depmap_baseline_summary.tsv` and `depmap24q2_baseline.tsv`.

## Nalm6 CD19 comparison

Place `GSE130819_short_term_mageck.count.csv.gz`, `GSE130819_long_term_mageck.count.csv.gz` and `GSE130819_series_matrix.txt.gz` in `inputs/gse130819/`. Keep source count-column names intact.

```bash
python code/ancillary/nalm6_boundary_recheck.py inputs/gse130819 generated/nalm6
```

The released summary retains the day-1 descriptive pool. The script also records the separate day-36 archived pair and its control-label discrepancy.

## Fröse in vivo and Frey source checks

Under a local input root, place the Fröse processed pool workbook at `01_data/frose_primary/AllPools_assembled.xlsx` and the Frey processed file at `01_data/cart_human_gse180834/GSE180834_gene_summary_allreps.txt.gz`. These are relative layout conventions for the existing comparison reader, not additional data distributed with the repository.

```bash
python code/ancillary/frose_pool_summary.py inputs/immune generated/frose
```

The main p1112 comparison is 9 positive hits among 923 eligible genes overall and 0 positive hits among the 20 covered GPI genes. The focused Fröse and U87 rows in the main figure are reproduced through the packaged source summaries. Their full source-file preparation routes are not included as independently tested commands in this release.

## PRISM and pathway annotation

The PRISM fixed Pearson coefficients and sample counts are given in `prism_pearson_ci.tsv`. `analysis_extras/check_summaries.py` checks the Fisher-z intervals from the rounded coefficients. Regenerating both Pearson correlations would require the original PRISM-to-expression and PRISM-to-dependency cell-line joins, including the exact 12-gene baseline dependency input. The aggregate results support figure and table checking. The release does not claim end-to-end reconstruction of those joins.

For KEGG, the package includes the small derived coverage/hit indicator table and verifies 26 of 30 aligned and 7 overlaps. The source pathway entry remains available from KEGG. The independent-study 0/27 coverage check refers to the source files identified in [SOURCES.md](SOURCES.md).
