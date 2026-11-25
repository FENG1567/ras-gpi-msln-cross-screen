# Verification record

Checks for the 2026-11-24 submission candidate used the public source-file copies identified in [SOURCES.md](SOURCES.md).

| Workflow | Reproduction check |
|---|---|
| Primary AsPC-1 and two MSLN screens | 18,001 eligible AsPC-1 genes, 16,872 aligned genes, 25/27 GPI coverage, and seven main-threshold GPI dual hits. Both MSLN positives and one-sided enrichment tails reproduced. Seven candidates' drug and CAR effect values matched the checked study outputs. |
| Primary threshold and column sensitivity | Both 16-setting grids matched the reference counts and hit lists in all 32 rows. Six label allocations and six column diagnostics agreed with the reference output. |
| Matched ranks | Both conditional summaries matched all 10 reference rows. Drug percentile 72.794 against 46.844 after quintile matching and 46.641 in the nearest-neighbour sensitivity. |
| Alternative AsPC-1/UM53 background | 15,782 genes, 22/27 GPI coverage, no null draw reached the observed count in each grid cell. The forward held-out count was 67 to 66 at K=200, plus-one p≈0.448. |
| Same-readout null simulation | Mean relative reduction of 83.44% for independent input scores and penalty weight one. This is a bias illustration rather than an independent test. |
| Nalm6 and Fröse/Frey supporting checks | Ran the portable readers on the public processed source files. Reproduced the p1112 scope of 923 eligible genes, nine overall positive hits and zero positive GPI hits among 20 covered GPI genes. The Nalm6 day-1 and day-36 results were reproduced separately. |
| DepMap 24Q2 baseline | Re-ran the row projection and summary on the checksum-verified pancreas-row extract and official Model.csv. The manuscript's group medians and 18,443/47/45 dimensions reproduced. This test used the retained rows of the release; the full upstream matrix remains an external download. |
| Figures 1–4 | The renderer produced four vector PDFs and SVGs with corresponding 300-dpi PNG previews from the checked-in figure-level tables. Font names may differ by system. |

The offline summary checker also verifies the alternative KEGG coverage and the approximate PRISM correlation intervals computed from reported coefficients. Remaining coverage and external-file requirements are listed in [SUPPLEMENTAL_WORKFLOWS.md](SUPPLEMENTAL_WORKFLOWS.md). Local working files and inputs are outside this repository.
