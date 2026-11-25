# Data sources and redistribution

The two primary input files are retrieved from their public repositories and are kept in the local, ignored `inputs/` directory. `code/README.md` specifies the expected files and the numerical checkpoints. The SHA-256 hashes below identify the copies used for this release. A repackaged download may have different byte-level hash even if its data are unchanged. Always verify the result anchors when using a differently packaged copy.

| Source | Input or published endpoint | Release reference | Role |
|---|---|---|---|
| Kumarasamy et al., GSE249537 | `GSE249537_ASPC1-UM53-TKOv3.count.txt.gz`, containing the four AsPC-1 profiles GSM7949449–GSM7949452 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249537), [count file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249537/suppl/GSE249537_ASPC1-UM53-TKOv3.count.txt.gz), [DOI](https://doi.org/10.1158/0008-5472.CAN-23-2504) | Raw input for the primary screen. UM53 contributes only to the separate coverage sensitivity in the manuscript. |
| Hagel et al. | Spreadsheet Supplementary Table S1, local input name `can22_2245_tableS1.xlsx` | [PMC article and supplementary files](https://pmc.ncbi.nlm.nih.gov/articles/PMC9929516/), [DOI](https://doi.org/10.1158/0008-5472.CAN-22-2245) | HuP-T3 and KP4_MSLN gene-level LFC and source screening scores. |
| Frey et al. | GSE180834, `GSE180834_gene_summary_allreps.txt.gz` | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180834), [DOI](https://doi.org/10.1038/s41467-022-29412-3) | KPC-OVA and OT-I comparison. |
| Fröse et al. | Focused hCD19 comparison and in vivo p1112 pool from the processed source files | [Public data](https://zenodo.org/records/21650788), [DOI](https://doi.org/10.1126/sciadv.aeh8973) | hCD19 comparison and focused-library coverage. |
| Singh et al. | GSE130819 short-term and long-term MAGeCK count matrices and series matrix | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130819), [DOI](https://doi.org/10.1158/2159-8290.CD-19-0813) | Nalm6 day 1 and the separately documented archival day 36 pair. |
| Larson et al. | GSE179147 processed U87 EGFR CAR screening files | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179147), [DOI](https://doi.org/10.1038/s41586-022-04585-5) | Descriptive comparison under the manuscript threshold rule. |
| Broad DepMap | 24Q2 `CRISPRGeneEffect.csv` and `Model.csv` | [Versioned release](https://doi.org/10.25452/figshare.plus.25880521.v1) | Baseline Chronos score summaries for pancreas and PAAD models. |
| Dilly et al. and CCLE | PRISM AUC in Supplementary Table S3 and CCLE GSE36133; dependency sensitivity used a 12-gene subset | [Dilly DOI](https://doi.org/10.1158/2159-8290.CD-24-0177), [CCLE GSE36133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36133) | Published exploratory baseline associations and Fisher-z intervals. |
| KEGG | Human pathway hsa00563, pathway database dated 2026-11-20, retrieved 2026-11-23 | [Entry](https://www.kegg.jp/pathway/hsa00563), [REST entry](https://rest.kegg.jp/get/hsa00563) | Alternate GPI set definition. The coverage summary is included. The original KEGG flat file is external. |
| Long et al. | Pa14C and PANC-1 supplementary Data Files 2 and 3 | [PMC article and files](https://pmc.ncbi.nlm.nih.gov/articles/PMC13053059/), [DOI](https://doi.org/10.1158/0008-5472.CAN-25-1835) | Coverage check. Each inspected table covered 0/27 panel genes. |
| Hallin et al. | Supplementary Tables 6 and 7, guide-library design | [DOI](https://doi.org/10.1038/s41591-022-02007-7) | Coverage check. The library tables covered 0/27 panel genes. |

## Primary input checksums

| Local copy | SHA-256 |
|---|---|
| `GSE249537_counts.txt.gz` | `4e4d965edb6ece3a75c60337fdb324b0ac38fa1290d8c8162b2e000de6428517` |
| `can22_2245_tableS1.xlsx` | `caa2c446e4504fef42742d9dbc9933eee3f2937eb479ee1a503d1b53cb13843` |

The full source workbooks, count matrices and licensed databases are obtained by the reader under their own source terms. The bundle contains original analysis code and selective figure-level or aggregate results, including the small set of cited source-derived candidate measurements needed to redraw Figure 1. Source attribution is retained in those summaries.
