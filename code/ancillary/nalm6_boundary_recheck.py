#!/usr/bin/env python3
"""Bounded descriptive recheck of the existing GSE130819 Nalm6 boundary."""
import csv
import argparse
import gzip
import json
import math
import os
import statistics
from collections import defaultdict
from bisect import bisect_left, bisect_right

P = None
D = None
GPI = set("PIGA PIGB PIGC PIGF PIGG PIGH PIGK PIGL PIGM PIGN PIGO PIGP PIGQ PIGS PIGT PIGU PIGV PIGW PIGY GPAA1 DPM1 DPM2 DPM3 MPDU1 PGAP1 PGAP2 PGAP3".split())
FOCAL = ["PIGC", "PIGF", "PIGM", "PIGP", "PIGS", "PIGT"]
FILES = {
    "short": "GSE130819_short_term_mageck.count.csv.gz",
    "long": "GSE130819_long_term_mageck.count.csv.gz",
}


def read_counts(filename):
    with gzip.open(os.path.join(D, filename), "rt", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        samples = header[2:]
        gene_counts = defaultdict(lambda: [0.0] * len(samples))
        n_rows = 0
        for row in reader:
            if len(row) < len(samples) + 2:
                continue
            n_rows += 1
            gene = row[1].strip()
            if not gene:
                continue
            for j, value in enumerate(row[2:2 + len(samples)]):
                try:
                    gene_counts[gene][j] += float(value)
                except ValueError:
                    pass
    return samples, dict(gene_counts), n_rows


def group_counts(samples, counts, members):
    """Sum count columns listed in members; retains gene-level raw counts."""
    idx = [samples.index(s) for s in members]
    return {g: sum(vals[i] for i in idx) for g, vals in counts.items()}


def effect_table(case, control):
    total_case = sum(case.values()) or 1.0
    total_control = sum(control.values()) or 1.0
    genes = sorted(set(case) | set(control))
    lfc = {
        g: math.log2(((case.get(g, 0.0) / total_case * 1e6) + 1.0) /
                     ((control.get(g, 0.0) / total_control * 1e6) + 1.0))
        for g in genes
    }
    ranked = sorted(lfc, key=lambda g: (-lfc[g], g))
    ranks = {g: i + 1 for i, g in enumerate(ranked)}
    return lfc, ranks, total_case, total_control


def summarize(name, case_label, control_label, case, control, note):
    lfc, ranks, case_total, control_total = effect_table(case, control)
    eligible = sorted(lfc)
    ascending = sorted(lfc.values())
    def effect_midrank(gene):
        value = lfc[gene]
        return (bisect_left(ascending, value) + bisect_right(ascending, value)) / (2 * len(ascending))
    covered = sorted(g for g in eligible if g.upper() in GPI)
    background = sorted(g for g in eligible if g.upper() not in GPI)
    overall_hits = sorted(g for g in eligible if lfc[g] >= 0.5)
    gpi_hits = sorted(g for g in covered if lfc[g] >= 0.5)
    bg_hits = sorted(g for g in background if lfc[g] >= 0.5)
    gpi_ranks = [ranks[g] for g in covered]
    bg_ranks = [ranks[g] for g in background]
    gpi_effects = [lfc[g] for g in covered]
    bg_effects = [lfc[g] for g in background]
    gpi_rows = [
        {"gene": g, "lfc_log2": lfc[g], "rank_descending": ranks[g],
         "rank_percentile": ranks[g] / len(eligible), "effect_midrank_percentile": effect_midrank(g),
         "hit_ge_0_5": lfc[g] >= 0.5}
        for g in sorted(covered, key=lambda x: (ranks[x], x))
    ]
    focal_rows = [
        {"gene": g, "eligible": g in lfc,
         "lfc_log2": lfc.get(g), "rank_descending": ranks.get(g),
         "rank_percentile": (ranks[g] / len(eligible)) if g in ranks else None,
         "hit_ge_0_5": (lfc[g] >= 0.5) if g in lfc else None}
        for g in FOCAL
    ]
    cd19_row = {"gene": "CD19", "eligible": "CD19" in lfc,
                "lfc_log2": lfc.get("CD19"), "rank_descending": ranks.get("CD19"),
                "rank_percentile": (ranks["CD19"] / len(eligible)) if "CD19" in ranks else None,
                "hit_ge_0_5": (lfc["CD19"] >= 0.5) if "CD19" in lfc else None}
    return {
        "contrast": name, "case_group": case_label, "control_group": control_label,
        "definition": note, "direction": "positive LFC means higher normalized abundance in case than control; candidate escape-direction enrichment only",
        "normalization": "sum MAGeCK raw sgRNA counts by gene, divide each group by its total counts across all measured genes, CPM scale, then log2((case CPM + 1)/(control CPM + 1)); pooled groups sum replicate raw counts before normalization",
        "rank_method": "descending LFC; exact ties broken by gene symbol alphabetically",
        "background_definition": "all measured gene identifiers outside the fixed 27-gene GPI set; no additional exclusions",
        "n_gene_rows_eligible": len(eligible), "n_gpi_covered": len(covered), "n_gpi_total": len(GPI),
        "overall_lfc_ge_0_5": len(overall_hits), "gpi_lfc_ge_0_5": len(gpi_hits),
        "background_lfc_ge_0_5": len(bg_hits), "gpi_hit_genes": gpi_hits,
        "median_gpi_lfc": statistics.median(gpi_effects) if gpi_effects else None,
        "median_gpi_rank_descending": statistics.median(gpi_ranks) if gpi_ranks else None,
        "median_gpi_rank_percentile": (statistics.median(gpi_ranks) / len(eligible)) if gpi_ranks else None,
        "median_gpi_effect_midrank_percentile": statistics.median(effect_midrank(g) for g in covered),
        "background_n": len(background), "background_median_lfc": statistics.median(bg_effects) if bg_effects else None,
        "background_median_rank_descending": statistics.median(bg_ranks) if bg_ranks else None,
        "background_median_rank_percentile": (statistics.median(bg_ranks) / len(eligible)) if bg_ranks else None,
        "gene_count_total_case": case_total, "gene_count_total_control": control_total,
        "all_gpi": gpi_rows, "focal_six": focal_rows, "CD19_positive_control": cd19_row,
        "descriptive_only": True,
    }


def read_geo_provenance():
    path = os.path.join(D, "GSE130819_series_matrix.txt.gz")
    lines = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("!Series_title", "!Series_geo_accession", "!Series_status",
                                "!Series_pubmed_id", "!Series_summary", "!Series_overall_design",
                                "!Series_supplementary_file", "!Series_relation", "!Sample_title",
                                "!Sample_geo_accession", "!Sample_characteristics_ch1",
                                "!Sample_description", "!Sample_data_processing")):
                lines.append(line.rstrip("\n"))
    return lines


def main():
    global P, D
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", help="Directory with processed GSE130819 count files and series matrix")
    parser.add_argument("output_dir", help="Directory for the JSON summary and notes")
    args = parser.parse_args()
    D = os.path.abspath(args.input_dir)
    P = os.path.abspath(args.output_dir)
    os.makedirs(P, exist_ok=True)
    provenance = read_geo_provenance()
    short_samples, short_counts, short_rows = read_counts(FILES["short"])
    long_samples, long_counts, long_rows = read_counts(FILES["long"])
    short_groups = {
        "d0_rep1": ["S2_d0_short_term_rep1"], "d0_rep2": ["S3_d0_short_term_rep2"],
        "car_rep1": ["S4_d1_P_CAR19_rep1"], "car_rep2": ["S5_d1_P_CAR19_rep2"],
        "ntd_rep1": ["S6_d1_P_NTD_rep1"], "ntd_rep2": ["S7_d1_P_NTD_rep2"],
    }
    sg = {k: group_counts(short_samples, short_counts, v) for k, v in short_groups.items()}
    contrasts = []
    contrasts.append(summarize("short_first_pair_recomputed", "S4_d1_P_CAR19_rep1", "S6_d1_P_NTD_rep1", sg["car_rep1"], sg["ntd_rep1"], "First-pair short comparison: CAR19 replicate 1 versus NTD replicate 1."))
    contrasts.append(summarize("short_CAR19_vs_NTD_rep1", "S4_d1_P_CAR19_rep1", "S6_d1_P_NTD_rep1", sg["car_rep1"], sg["ntd_rep1"], "Matched replicate 1 at 24 h: CART19 co-culture versus control T-cell (NTD) co-culture."))
    contrasts.append(summarize("short_CAR19_vs_NTD_rep2", "S5_d1_P_CAR19_rep2", "S7_d1_P_NTD_rep2", sg["car_rep2"], sg["ntd_rep2"], "Matched replicate 2 at 24 h: CART19 co-culture versus control T-cell (NTD) co-culture."))
    car_pooled = group_counts(short_samples, short_counts, ["S4_d1_P_CAR19_rep1", "S5_d1_P_CAR19_rep2"])
    ntd_pooled = group_counts(short_samples, short_counts, ["S6_d1_P_NTD_rep1", "S7_d1_P_NTD_rep2"])
    d0_pooled = group_counts(short_samples, short_counts, ["S2_d0_short_term_rep1", "S3_d0_short_term_rep2"])
    contrasts.append(summarize("short_CAR19_vs_NTD_pooled", "S4+S5 d1_P_CAR19 (n=2)", "S6+S7 d1_P_NTD (n=2)", car_pooled, ntd_pooled, "Pool raw counts across both d1 replicates per group, then normalize. Primary replicate-inclusive short contrast."))
    for rep in (1, 2):
        contrasts.append(summarize(f"short_CAR19_vs_d0_rep{rep}", f"S{3+rep}_d1_P_CAR19_rep{rep}", f"S{1+rep}_d0_short_term_rep{rep}", sg[f"car_rep{rep}"], sg[f"d0_rep{rep}"], "Matched replicate-number d1 CART19 versus d0 no-T-cell baseline; time/culture change and CAR pressure are combined."))
        contrasts.append(summarize(f"short_NTD_vs_d0_rep{rep}", f"S{5+rep}_d1_P_NTD_rep{rep}", f"S{1+rep}_d0_short_term_rep{rep}", sg[f"ntd_rep{rep}"], sg[f"d0_rep{rep}"], "Matched replicate-number d1 control-T co-culture versus d0 no-T-cell baseline; time/culture change is included."))
    contrasts.append(summarize("short_CAR19_vs_d0_pooled", "S4+S5 d1_P_CAR19 (n=2)", "S2+S3 d0 (n=2)", car_pooled, d0_pooled, "Pool raw counts across both d1 CART19 and both d0 baseline replicates, then normalize; descriptive baseline contrast."))
    contrasts.append(summarize("short_NTD_vs_d0_pooled", "S6+S7 d1_P_NTD (n=2)", "S2+S3 d0 (n=2)", ntd_pooled, d0_pooled, "Pool raw counts across both d1 NTD and both d0 baseline replicates, then normalize; descriptive baseline contrast."))
    long_case = group_counts(long_samples, long_counts, ["S9_d36_P_CAR19"])
    long_ctrl = group_counts(long_samples, long_counts, ["S8_d36_alone"])
    contrasts.append(summarize("long_d36_single_pair", "S9_d36_P_CAR19 (n=1)", "S8_d36_alone raw label (n=1)", long_case, long_ctrl, "Retains the original single pair. GEO sample metadata describes S8 as 36-day control-T-cell co-culture; the processed count header calls it alone."))


    result = {
        "analysis": "Bounded descriptive Nalm6 boundary recheck; no other boundaries analyzed.",
        "dataset": "GSE130819; human Nalm6 B-ALL; Brunello genome-wide CRISPR sgRNA count tables.",
        "raw_inputs": {"short": FILES["short"], "long": FILES["long"], "geo_series_matrix": "GSE130819_series_matrix.txt.gz"},
        "provenance_lines": provenance,
        "sample_headers": {"short": short_samples, "long": long_samples},
        "raw_sgrna_rows": {"short": short_rows, "long": long_rows},
        "replicate_design": {
            "short": {"d0_baseline": ["S2_d0_short_term_rep1", "S3_d0_short_term_rep2"],
                      "d1_CART19": ["S4_d1_P_CAR19_rep1", "S5_d1_P_CAR19_rep2"],
                      "d1_NTD_control_T": ["S6_d1_P_NTD_rep1", "S7_d1_P_NTD_rep2"]},
            "long": {"d36_CART19": ["S9_d36_P_CAR19"], "d36_control_raw_header_alone": ["S8_d36_alone"]},
            "geo_design": "Series matrix describes 24 h and 36 d co-cultures with CART19 or non-transduced/control T cells, each condition in duplicate; processed long table contains one d36 CART19 and one d36 control-labeled column plus a d0 column."
        },
        "gpi_set": sorted(GPI), "focal_six": FOCAL,
        "threshold": "descriptive positive-direction LFC >= 0.5; no inferential test attached",
        "contrasts": contrasts,
        "interpretation_limits": [
            "LFC estimates are descriptive normalized count contrasts from sgRNA count tables, not replicate-aware MAGeCK MLE estimates or inferential tests.",
            "The >=0.5 threshold is a screening description, not evidence of statistical significance or absence of a biological effect.",
            "The two d1 CART19-versus-NTD replicate pairs are compared separately and by pooled raw counts; pooled values are not an estimate of between-replicate uncertainty.",
            "d0-referenced comparisons mix 24-hour culture/time effects with T-cell or CAR exposure and are included only to expose baseline behavior."
        ]
    }
    with open(os.path.join(P, "nalm6_boundary_results.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2, ensure_ascii=False)

    lines = ["# GSE130819 Nalm6 boundary recheck", "",
             "Scope: existing GSE130819 data only; no other boundary analyses. Counts were aggregated by gene and normalized as stated in the JSON. LFC >= 0.5 is descriptive; no inferential test is attached.", "",
             "## Sample provenance", "",
             "The GEO series matrix identifies GSE130819 as a Brunello genome-wide CRISPR screen in human Nalm6 B-ALL under CART19 pressure (PMID 32001516; also lists PMID 35690611), with short 24-hour and long 36-day co-cultures. The series-level design says control/non-transduced T-cell comparator conditions were run in duplicate. Short counts contain two d0, two d1 CART19, and two d1 NTD columns. Long counts contain S9 d36 CART19 and S8 d36_alone raw header labels only; GEO sample metadata describes S8 as 36-day control-T-cell co-culture. The original note's 'alone' wording therefore follows the processed header, not the GEO treatment annotation.", "",
             "## Results", "",
             "| Contrast | Eligible genes | GPI covered | Total LFC >=0.5 | GPI LFC >=0.5 | Median GPI rank | Median GPI LFC | Median background rank | Median background LFC | CD19 LFC / rank |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for c in contrasts:
        cd19 = c["CD19_positive_control"]
        cd19_value = (f"{cd19['lfc_log2']:+.3f} / {cd19['rank_descending']}"
                      if cd19["eligible"] else "not measured")
        lines.append(f"| {c['contrast']} | {c['n_gene_rows_eligible']} | {c['n_gpi_covered']}/{c['n_gpi_total']} | {c['overall_lfc_ge_0_5']} | {c['gpi_lfc_ge_0_5']} | {c['median_gpi_rank_descending']:.1f} ({c['median_gpi_rank_percentile']:.3f}) | {c['median_gpi_lfc']:+.3f} | {c['background_median_rank_descending']:.1f} ({c['background_median_rank_percentile']:.3f}) | {c['background_median_lfc']:+.3f} | {cd19_value} |")
    lines.extend(["", "The first-pair short comparison is reported alongside the pooled replicate summary.", "",
                  "The original short result used S4 (CAR19 rep1) versus S6 (NTD rep1); the updated primary short summary pools raw counts from both CAR19 replicates (S4+S5) against both NTD replicates (S6+S7), while also retaining each matched replicate contrast. d0 contrasts are shown separately as culture/time baselines and should not be read as CAR-specific effects.", "",
                  "## Requested six PIG genes", "",
                  "| Gene | CAR19 vs NTD rep1 LFC / rank | CAR19 vs NTD rep2 LFC / rank | Pooled CAR19 vs NTD LFC / rank | Long pair LFC / rank |", "|---|---:|---:|---:|---:|"])
    key_contrasts = {c["contrast"]: c for c in contrasts}
    for gene in FOCAL:
        vals = []
        for contrast_name in ("short_CAR19_vs_NTD_rep1", "short_CAR19_vs_NTD_rep2", "short_CAR19_vs_NTD_pooled", "long_d36_single_pair"):
            item = next(x for x in key_contrasts[contrast_name]["focal_six"] if x["gene"] == gene)
            vals.append(f"{item['lfc_log2']:+.3f} / {item['rank_descending']}")
        lines.append("| " + gene + " | " + " | ".join(vals) + " |")
    lines.extend(["",
                  "## Requested GPI details", "",
                  "The JSON contains all 27 GPI genes' LFC, descending rank, rank percentile, and threshold status for each contrast, plus the background median LFC/rank. It also contains PIGC, PIGF, PIGM, PIGP, PIGS, and PIGT estimates and ranks for each contrast. All values are descriptive; no claim of no effect follows from zero thresholded hits.", "",
                  "## Reproducibility", "",
                  "Source: nalm6_boundary_recheck.py. Inputs are the public GSE130819 processed count tables and GEO series matrix; see SOURCES.md for access details.", ""])
    with open(os.path.join(P, "nalm6_boundary_recheck.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(json.dumps({"output": os.path.join(P, "nalm6_boundary_results.json"),
                      "notes": os.path.join(P, "nalm6_boundary_recheck.md"),
                      "short_rows": short_rows, "long_rows": long_rows,
                      "short_samples": short_samples, "long_samples": long_samples,
                      "contrasts": [{"name": c["contrast"], "genes": c["n_gene_rows_eligible"], "gpi_covered": c["n_gpi_covered"], "gpi_hits": c["gpi_lfc_ge_0_5"], "all_hits": c["overall_lfc_ge_0_5"], "median_gpi_rank": c["median_gpi_rank_descending"], "median_gpi_lfc": c["median_gpi_lfc"]} for c in contrasts]}, indent=2))


if __name__ == "__main__":
    main()
