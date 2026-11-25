# GPI-anchor biosynthesis in pancreatic cancer screening

This repository accompanies *GPI-anchor biosynthesis in pancreatic cancer CRISPR screens: overlap between relative depletion under MRTX1133 and mesothelin CAR T-cell escape*. It contains the authors' code and a limited set of derived results associated with the manuscript.

## Quick check and figure generation

Run from the repository root with Python 3.10 or later.

```bash
python -m pip install -r requirements.txt
python analysis_extras/check_summaries.py
python code/render_figures.py --output-dir generated/figures
python analysis_extras/ranking_bias.py --output generated/ranking_bias.csv
```

The summary check verifies the manuscript's Fisher tails, p1112 denominators, GPI threshold grid and both matched-rank anchors. The figure command redraws Figures 1–4 from the released figure-level tables. It creates vector PDF and SVG files and 300-dpi TIFF and PNG files. The figure names follow the manuscript order. The renderer uses an available sans-serif font by default. Pass `--arial-bold /path/to/arialbd.ttf` with your own licensed Arial Bold font to match the manuscript typography. The ranking-bias simulation regenerates the supplement's same-readout illustration of approximately 83% without external data.

## Primary analysis from public inputs

Obtain the GSE249537 raw count file and Hagel et al. Supplementary Table S1 from their original sources. See [SOURCES.md](SOURCES.md) for input identification, download locations and the checksums of the copies used. Keep downloaded data in the ignored `inputs/` directory. For example:

```bash
python code/reproduce_core.py \
  --drug-counts inputs/GSE249537_ASPC1-UM53-TKOv3.count.txt.gz \
  --car-table-s1 inputs/can22_2245_tableS1.xlsx \
  --output-dir generated/core
```

See [the core-analysis instructions](code/README.md) for the design, outputs and checkpoints. The raw primary workflow recomputes AsPC-1 eligibility, the MSLN escape readouts, dual hits, the full threshold grid, and both the quintile and nearest-neighbour matched-rank analyses. The manuscript's primary anchors are 18,001 eligible genes, 16,872 aligned genes, 25/27 GPI coverage and seven primary GPI dual hits; the primary matched rank comparison is 72.8 vs 46.8 percentile and the nearest-neighbour sensitivity is 72.8 vs 46.6.

The separate AsPC-1/UM53 coverage analysis and held-out ranking comparison use the same two public input files. Run `python analysis_extras/alternative_eligibility.py --drug-counts inputs/GSE249537_ASPC1-UM53-TKOv3.count.txt.gz --car-table-s1 inputs/can22_2245_tableS1.xlsx`. The distinct source-native gene-symbol alignment reproduces the 15,782-gene universe, 22/27 GPI coverage and the 67-to-66 held-out count with p≈0.448. This background is separate from the primary AsPC-1-only analysis.

## Package contents

- `code/` and `analysis_extras/` contain original analysis, checking and figure-generation code.
- `data/derived/` contains the checked-in author summaries supporting the figures, QC, alternative thresholds, pathway coverage, baseline associations and the ranking check. The few externally sourced effect values used as Figure 1 coordinates have the attribution and terms described in [SOURCES.md](SOURCES.md).
- `README.md`, `SOURCES.md`, `CITATION.cff`, `requirements.txt`, `.gitignore`, `LICENSE` and `LICENSE-DATA.md` describe how to use and attribute the package.

Keep third-party raw files and generated outputs in the corresponding ignored directories. The offline figure workflow uses released summaries for ancillary estimates. [Supplemental workflow notes](SUPPLEMENTAL_WORKFLOWS.md) distinguish the included external-data utilities from checks using released summaries.

## License

The authors' code is [MIT licensed](LICENSE). The authors' original derived summaries are [CC BY 4.0 licensed](LICENSE-DATA.md). Source-data rights remain with the applicable providers. `CITATION.cff` identifies the authors and prepared package version.
