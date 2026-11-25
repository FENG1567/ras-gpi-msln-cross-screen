# Reproduce the AsPC-1 / MSLN core screen results

`reproduce_core.py` recalculates the AsPC-1 MRTX1133 effect from the raw guide-count matrix, reads the published HuP-T3 and KP4 MSLN escape scores, applies the declared eligibility and hit thresholds, and runs the CPM-decile nearest-neighbour rank comparison. It uses no historical UM53 score file or project-local intermediate table.

## Requirements

- Python 3.10 or later
- `openpyxl` 3.1 or later

All other imports are from the Python standard library. No network connection is needed during analysis.

## Source files

The source files are deliberately not bundled. Obtain them from the original public sources and check that their SHA-256 values match before running:

| Input | Public source | Expected SHA-256 |
| --- | --- | --- |
| `GSE249537_counts.txt.gz` | GEO series [GSE249537](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249537), supplementary files | `4e4d965edb6ece3a75c60337fdb324b0ac38fa1290d8c8162b2e000de6428517` |
| `can22_2245_tableS1.xlsx` | Hagel et al., *Cancer Research* (2023), Supplementary Table S1, [article record](https://doi.org/10.1158/0008-5472.CAN-22-2245) | `caa2c446e4504fef42742d9dcbc9933eee3f2937eb479ee1a503d1b53cb13843` |

The code validates the source workbook's gene counts and locked positive-hit totals. `qc_summary.json` records only input basenames and hashes, not local absolute paths. The input files remain subject to their source terms; this repository does not grant a new license to third-party inputs.

## Run

Place the downloaded files in any local directories, then pass their paths explicitly:

```bash
python code/reproduce_core.py \
  --drug-counts inputs/GSE249537_counts.txt.gz \
  --car-table-s1 inputs/can22_2245_tableS1.xlsx \
  --output-dir results
```

PowerShell:

```powershell
python .\code\reproduce_core.py `
  --drug-counts .\inputs\GSE249537_counts.txt.gz `
  --car-table-s1 .\inputs\can22_2245_tableS1.xlsx `
  --output-dir .\results
```

## Analysis definitions

1. Normalize each of the four AsPC-1 profiles separately to counts per million (CPM), sum guide counts by gene symbol, and define baseline expression as the mean CPM in the two DMSO profiles. Retain genes with mean DMSO CPM ≥5.
2. The primary drug score is `log2[(mean gene-summed DMSO CPM + 1)/(mean gene-summed MRTX1133 CPM + 1)]`; positive values indicate relative sensitization. The robustness score is the median sgRNA-level `log2[(mean per-guide DMSO CPM + 0.5)/(mean per-guide MRTX1133 CPM + 0.5)]`.
3. Read the first three columns of each HuP-T3/KP4 worksheet as gene symbol, escape log2 fold-change, and the source averaged −log10(p) score. Duplicate gene symbols use the last source row, matching the locked source reader.
4. The baseline dual-hit rule is drug score ≥0.5 and, in either CAR screen, escape log2 fold-change ≥0.5 plus source averaged −log10(p) ≥2. The script also reports the full 4 × 2 × 2 threshold grid for both drug-score definitions.
5. For the descriptive rank comparison, calculate midrank percentiles within the genes passing the AsPC-1 CPM filter and present in both CAR sources. For each covered GPI gene, form a pool of its 25 closest non-GPI genes in the same AsPC-1 DMSO CPM decile, using Euclidean distance across HuP-T3 and KP4 CAR rank percentiles. Draw one unique non-GPI comparison per GPI gene for 20,000 iterations (seed 261125); report the empirical plus-one tail.

## Outputs

- `eligible_gene_effects.tsv`: three-source gene estimates, guide summaries, and midrank percentiles.
- `primary_dual_gene_results.tsv`: all primary dual hits and GPI panel genes, with eligibility and per-guide concordance.
- `threshold_sensitivity.tsv`: all declared thresholds for both drug-score methods.
- `car_source_enrichment.tsv`: source-level MSLN hit counts and one-sided hypergeometric enrichment calculations.
- `conditional_drug_rank_matching.tsv`: the two locked matched-rank designs, first matching on AsPC-1 DMSO CPM deciles alone, then matching on CPM deciles plus HuP-T3 and KP4 escape-rank quintiles.
- `matched_rank_summary.json` and `matched_pool_support.tsv`: a separate nearest-neighbour sensitivity comparison using the 25 closest non-GPI genes per GPI gene within each CPM decile.
- `aspc1_all_2v2_allocations.tsv` and `aspc1_column_loo_qc.tsv`: six treatment/control column allocations and six column-level sensitivity diagnostics. These describe the observed profiles and do not establish independent biological replication.
- `qc_summary.json`: source hashes, input checks, key counts, definitions, and interpretation limits.

## Interpretation limits

This code reanalyses existing screens; it does not create an independent experiment. The AsPC-1 comparison has two treatment and two control profiles, and public metadata do not establish their biological independence. Guide-level rows are nested within genes. The nearest-neighbour tail is a descriptive conditional summary; it does not correct post-hoc pathway selection or gene dependence. CAR p-score fields are source values, not new pooled tests.

