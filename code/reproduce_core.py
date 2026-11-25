#!/usr/bin/env python3
"""Reproduce the AsPC-1 MRTX1133 / MSLN CAR-T core screen analyses.

The script accepts the public source files as inputs and writes derived tables
and a machine-readable QC summary. It intentionally does not use the historical
UM53 score file or other private/derived intermediate files.

Dependencies: Python >=3.10 and openpyxl >=3.1.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import random
import statistics
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from openpyxl import load_workbook


ASPC1_COLUMNS = (
    "ASPC1-DMSO-p11",
    "ASPC1-MRTX-6-p11",
    "ASPC1-DMSO-p11-A",
    "ASPC1-MRTX-6-p11-A",
)
GPI_GENES = frozenset(
    "PIGA PIGB PIGC PIGF PIGG PIGH PIGK PIGL PIGM PIGN PIGO PIGP PIGQ "
    "PIGS PIGT PIGU PIGV PIGW PIGY GPAA1 DPM1 DPM2 DPM3 MPDU1 PGAP1 "
    "PGAP2 PGAP3".split()
)
FROZEN_DUAL_GPI = frozenset("PIGC PIGF PIGM PIGP PIGS PIGT".split())
BASELINE = (0.5, 0.5, 2.0)
THRESHOLDS = tuple(itertools.product((0.3, 0.5, 0.75, 1.0), (0.3, 0.5), (1.3, 2.0)))
MATCH_METRICS = (
    "gene_summed_drug_rank_percentile",
    "per_guide_drug_rank_percentile",
    "HuP_T3_rank_percentile",
    "KP4_rank_percentile",
)
MATCH_DRAWS = 20_000
MATCH_SEED = 261125
MATCH_NEIGHBOURS = 25


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def quantile(values: Iterable[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot calculate a quantile of an empty list")
    location = (len(ordered) - 1) * probability
    lower = int(location)
    upper = math.ceil(location)
    return ordered[lower] + (location - lower) * (ordered[upper] - ordered[lower])


def midrank_percentiles(values: dict[str, float]) -> dict[str, float]:
    """Midrank percentiles on [0, 100], retaining the locked tie convention."""
    if not values:
        return {}
    ordered = sorted(values.values())
    n = len(ordered)
    return {
        gene: 100.0 * (bisect_left(ordered, value) + bisect_right(ordered, value)) / (2 * n)
        for gene, value in values.items()
    }


def read_aspc1_counts(path: Path) -> tuple[dict, dict]:
    """Read the four locked AsPC-1 columns and calculate gene/guide effects."""
    counts_by_gene: dict[str, list[float]] = defaultdict(lambda: [0.0] * len(ASPC1_COLUMNS))
    guides: list[tuple[str, str, list[float]]] = []
    totals = [0.0] * len(ASPC1_COLUMNS)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Count file is empty: {path}") from exc
        missing = [name for name in ASPC1_COLUMNS if name not in header]
        if missing:
            raise ValueError(f"Missing required AsPC-1 count columns: {missing}")
        indexes = [header.index(name) for name in ASPC1_COLUMNS]
        for line_number, row in enumerate(reader, start=2):
            if len(row) <= max(indexes) or len(row) < 2:
                raise ValueError(f"Malformed count row at line {line_number}")
            guide_id = row[0].strip()
            gene = row[1].strip().upper()
            if not gene:
                continue
            try:
                values = [float(row[index]) for index in indexes]
            except ValueError as exc:
                raise ValueError(f"Non-numeric count at line {line_number}") from exc
            if any(not math.isfinite(value) or value < 0 for value in values):
                raise ValueError(f"Invalid count at line {line_number}: counts must be finite and non-negative")
            for index, value in enumerate(values):
                totals[index] += value
                counts_by_gene[gene][index] += value
            guides.append((guide_id, gene, values))

    if not guides or any(total <= 0 for total in totals):
        raise ValueError("No usable guide counts or a zero-size library was found")

    cpm_by_gene = {
        gene: [1_000_000.0 * value / total for value, total in zip(values, totals)]
        for gene, values in counts_by_gene.items()
    }
    per_gene_column_cpm = cpm_by_gene
    dmso_cpm = {gene: statistics.mean((values[0], values[2])) for gene, values in cpm_by_gene.items()}
    gene_sum_effect = {
        gene: math.log2((statistics.mean((values[0], values[2])) + 1.0) /
                        (statistics.mean((values[1], values[3])) + 1.0))
        for gene, values in cpm_by_gene.items()
    }

    guide_effects: dict[str, list[float]] = defaultdict(list)
    for _, gene, values in guides:
        per_column_cpm = [1_000_000.0 * value / total for value, total in zip(values, totals)]
        effect = math.log2((statistics.mean((per_column_cpm[0], per_column_cpm[2])) + 0.5) /
                           (statistics.mean((per_column_cpm[1], per_column_cpm[3])) + 0.5))
        guide_effects[gene].append(effect)
    guide_summary = {
        gene: {
            "n_guides": len(effects),
            "median": statistics.median(effects),
            "q1": quantile(effects, 0.25),
            "q3": quantile(effects, 0.75),
            "positive_guides": sum(value > 0 for value in effects),
        }
        for gene, effects in guide_effects.items()
    }
    qc = {
        "sample_columns": list(ASPC1_COLUMNS),
        "library_totals": dict(zip(ASPC1_COLUMNS, totals)),
        "guide_rows_with_gene_symbol": len(guides),
        "unique_gene_symbols": len(counts_by_gene),
        "dmso_cpm_definition": "arithmetic mean of CPM in ASPC1-DMSO-p11 and ASPC1-DMSO-p11-A after separate library-size normalization",
        "primary_drug_effect_definition": "log2[(mean gene-summed DMSO CPM + 1)/(mean gene-summed MRTX1133 CPM + 1)]; positive values indicate relative sensitization",
        "guide_robustness_definition": "median across sgRNAs of log2[(mean per-guide DMSO CPM + 0.5)/(mean per-guide MRTX1133 CPM + 0.5)]",
    }
    return {
        "dmso_cpm": dmso_cpm,
        "gene_sum_effect": gene_sum_effect,
        "guide_summary": guide_summary,
        "per_gene_column_cpm": per_gene_column_cpm,
    }, qc


def read_car_workbook(path: Path) -> tuple[dict, dict]:
    """Read the locked first three columns of the HuP-T3 and KP4 worksheets."""
    workbook = load_workbook(path, read_only=True, data_only=True)
    screens: dict[str, dict[str, tuple[float, float]]] = {}
    qc = {}
    try:
        for screen, token in (("HuP-T3", "hupt3"), ("KP4", "kp4")):
            names = [name for name in workbook.sheetnames
                     if token in name.lower().replace("-", "").replace(" ", "")]
            if len(names) != 1:
                raise ValueError(f"Could not uniquely identify {screen} worksheet in {workbook.sheetnames}")
            rows = workbook[names[0]].values
            header = next(rows, None)
            if header is None or len(header) < 3:
                raise ValueError(f"Missing gene/effect/source-score columns in {names[0]}")
            scores: dict[str, tuple[float, float]] = {}
            duplicates = 0
            for row in rows:
                if not row or len(row) < 3 or not row[0]:
                    continue
                gene = str(row[0]).strip().upper()
                try:
                    effect, neglog10p = float(row[1]), float(row[2])
                except (TypeError, ValueError):
                    continue
                if not math.isfinite(effect) or not math.isfinite(neglog10p):
                    continue
                if gene in scores:
                    duplicates += 1
                # Preserve the existing locked reader's source-order overwrite rule.
                scores[gene] = (effect, neglog10p)
            source_hits = {gene for gene, (effect, pscore) in scores.items()
                           if effect >= 0.5 and pscore >= 2.0}
            gpi_hits = source_hits & GPI_GENES
            expected = {
                "HuP-T3": (20_080, 173, 18),
                "KP4": (20_080, 26, 16),
            }[screen]
            observed = (len(scores), len(source_hits), len(gpi_hits))
            if observed != expected:
                raise ValueError(f"Locked-source check failed for {screen}: {observed}; expected {expected}")
            screens[screen] = scores
            qc[screen] = {
                "worksheet": names[0],
                "header": [str(value) if value is not None else "" for value in header],
                "unique_genes": len(scores),
                "duplicate_gene_rows_overwritten_in_source_order": duplicates,
                "source_cutoff_gene_hits": len(source_hits),
                "gpi_panel_coverage": len(GPI_GENES & set(scores)),
                "gpi_source_cutoff_hits": len(gpi_hits),
                "source_cutoff": "escape log2 fold-change >= 0.5 and source averaged -log10(p) >= 2",
            }
    finally:
        workbook.close()
    return screens, qc


def hypergeom_upper_tail(population: int, successes: int, draws: int, observed: int) -> float:
    """One-sided hypergeometric tail, evaluated in log space without SciPy."""
    low = max(observed, 0, draws - (population - successes))
    high = min(successes, draws)
    if low > high:
        return 0.0

    def log_choose(n: int, k: int) -> float:
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    log_denominator = log_choose(population, draws)
    log_terms = [log_choose(successes, k) + log_choose(population - successes, draws - k) - log_denominator
                 for k in range(low, high + 1)]
    peak = max(log_terms)
    return math.exp(peak) * sum(math.exp(value - peak) for value in log_terms)


def nearest_match_analysis(eligible_rows: list[dict]) -> tuple[dict, list[dict]]:
    """Run the locked CPM-decile, 25-neighbour descriptive rank comparison."""
    bases = {row["gene"]: float(row["dmsocpm_mean_of_normalized_columns"]) for row in eligible_rows}
    values = {row["gene"]: tuple(float(row[key]) for key in MATCH_METRICS) for row in eligible_rows}
    if len(bases) != len(eligible_rows):
        raise ValueError("Duplicate genes in the three-source AsPC-1/CAR table")
    base_values = sorted(bases.values())
    cutpoints = [quantile(base_values, decile / 10) for decile in range(1, 10)]
    gpi = sorted(row["gene"] for row in eligible_rows if row["in_gpi"] == 1)
    bins: dict[int, list[str]] = defaultdict(list)
    for row in eligible_rows:
        if row["in_gpi"] == 0:
            bins[bisect_left(cutpoints, bases[row["gene"]])].append(row["gene"])
    pools: dict[str, list[str]] = {}
    pool_rows = []
    for gene in gpi:
        decile = bisect_left(cutpoints, bases[gene])
        candidates = bins[decile]
        if len(candidates) < max(MATCH_NEIGHBOURS, len(gpi)):
            raise ValueError(f"Insufficient non-GPI matches in DMSO CPM stratum for {gene}")
        _, _, hup_rank, kp4_rank = values[gene]

        def distance(candidate: str) -> tuple[float, str]:
            candidate_hup, candidate_kp4 = values[candidate][2:4]
            return ((candidate_hup - hup_rank) ** 2 + (candidate_kp4 - kp4_rank) ** 2, candidate)

        pools[gene] = sorted(candidates, key=distance)[:MATCH_NEIGHBOURS]
        furthest = distance(pools[gene][-1])[0] ** 0.5
        pool_rows.append({
            "gene": gene,
            "cpm_decile": decile + 1,
            "same_decile_non_gpi_candidates": len(candidates),
            "neighbours_per_gene": MATCH_NEIGHBOURS,
            "max_two_car_rank_distance_among_neighbours": furthest,
        })

    random_generator = random.Random(MATCH_SEED)
    simulated = [[] for _ in MATCH_METRICS]
    for _ in range(MATCH_DRAWS):
        sampled: set[str] = set()
        sums = [0.0] * len(MATCH_METRICS)
        for gene in gpi:
            comparison = random_generator.choice([candidate for candidate in pools[gene] if candidate not in sampled])
            sampled.add(comparison)
            for index, value in enumerate(values[comparison]):
                sums[index] += value
        for index, value in enumerate(sums):
            simulated[index].append(value / len(gpi))

    results = []
    for index, metric in enumerate(MATCH_METRICS):
        observed = statistics.mean(values[gene][index] for gene in gpi)
        null_values = sorted(simulated[index])
        reached = sum(value >= observed for value in null_values)
        results.append({
            "metric": metric,
            "observed_mean_percentile": observed,
            "matched_null_median": statistics.median(null_values),
            "matched_null_q025": quantile(null_values, 0.025),
            "matched_null_q975": quantile(null_values, 0.975),
            "draws_at_or_above_observed": reached,
            "draws": MATCH_DRAWS,
            "seed": MATCH_SEED,
            "plus_one_tail_descriptive": (reached + 1) / (MATCH_DRAWS + 1),
        })
    return {
        "gpi_genes_matched": len(gpi),
        "neighbors_per_gpi_gene": MATCH_NEIGHBOURS,
        "matching": "within AsPC-1 DMSO CPM decile; choose 25 nearest non-GPI genes by Euclidean distance across HuP-T3 and KP4 CAR rank percentiles; sample without replacement across GPI genes",
        "draws": MATCH_DRAWS,
        "seed": MATCH_SEED,
        "interpretation": "Descriptive conditional comparison; not a population-level test and does not correct pathway selection or gene dependence.",
        "metrics": results,
    }, pool_rows


def quantile_bins(values: dict[str, float], n_bins: int) -> dict[str, int]:
    cutpoints = [quantile(values.values(), index / n_bins) for index in range(1, n_bins)]
    return {gene: bisect_left(cutpoints, value) for gene, value in values.items()}


def conditional_rank_matching(eligible_rows: list[dict]) -> list[dict]:
    """Reproduce the CPM-only and CPM+both-CAR-rank conditional matches."""
    genes = [row["gene"] for row in eligible_rows]
    covered = sorted(row["gene"] for row in eligible_rows if row["in_gpi"] == 1)
    if not covered:
        raise ValueError("The three-source universe has no covered GPI genes")
    non_gpi = sorted(set(genes) - set(covered))
    cpm = {row["gene"]: float(row["dmsocpm_mean_of_normalized_columns"]) for row in eligible_rows}
    cpm_bins = quantile_bins(cpm, 10)
    hup_rank = {row["gene"]: float(row["HuP_T3_rank_percentile"]) for row in eligible_rows}
    kp4_rank = {row["gene"]: float(row["KP4_rank_percentile"]) for row in eligible_rows}
    metrics = {
        "gene_summed_D_drug_rank": {row["gene"]: float(row["gene_summed_drug_rank_percentile"]) for row in eligible_rows},
        "per_guide_median_drug_rank": {row["gene"]: float(row["per_guide_drug_rank_percentile"]) for row in eligible_rows},
        "HuP_T3_escape_rank": hup_rank,
        "KP4_escape_rank": kp4_rank,
        "higher_CAR_escape_rank": {gene: max(hup_rank[gene], kp4_rank[gene]) for gene in genes},
    }

    def run(strata: dict[str, object], label: str, seed: int) -> list[dict]:
        counts: dict[object, int] = defaultdict(int)
        pools: dict[object, list[str]] = defaultdict(list)
        for gene in covered:
            counts[strata[gene]] += 1
        for gene in non_gpi:
            pools[strata[gene]].append(gene)
        unsupported = {
            stratum: {"required": count, "available": len(pools[stratum])}
            for stratum, count in counts.items()
            if len(pools[stratum]) < count
        }
        if unsupported:
            raise ValueError(f"Conditional matched comparison lacks common support: {unsupported}")
        observed = {name: statistics.mean(values[gene] for gene in covered) for name, values in metrics.items()}
        simulations = {name: [] for name in metrics}
        random_generator = random.Random(seed)
        for _ in range(MATCH_DRAWS):
            sample = [gene for stratum, count in counts.items()
                      for gene in random_generator.sample(pools[stratum], count)]
            for name, values in metrics.items():
                simulations[name].append(statistics.mean(values[gene] for gene in sample))
        rows = []
        for name, values in simulations.items():
            null_values = sorted(values)
            reached = sum(value >= observed[name] for value in null_values)
            rows.append({
                "statistic": name,
                "observed_gpi_mean_percentile": observed[name],
                "null_median": statistics.median(null_values),
                "null_q025": quantile(null_values, 0.025),
                "null_q975": quantile(null_values, 0.975),
                "draws_at_least_observed": reached,
                "plus_one_tail_descriptive": (reached + 1) / (MATCH_DRAWS + 1),
                "gpi_n": len(covered),
                "draws": MATCH_DRAWS,
                "matching": label,
                "strata_n": len(counts),
                "supported": True,
                "seed": seed,
            })
        return rows

    cpm_only = {gene: cpm_bins[gene] for gene in genes}
    hup_bins = quantile_bins(hup_rank, 5)
    kp4_bins = quantile_bins(kp4_rank, 5)
    cpm_and_car = {
        gene: (
            cpm_bins[gene],
            hup_bins[gene],
            kp4_bins[gene],
        )
        for gene in genes
    }
    return run(cpm_only, "AsPC1_DMSO_CPM_decile", 261123) + run(
        cpm_and_car,
        "AsPC1_DMSO_CPM_decile_plus_HuPT3_and_KP4_escape_rank_quintiles",
        261124,
    )


def aspc1_column_diagnostics(per_gene_column_cpm: dict[str, list[float]], genes: list[str], gpi: list[str]) -> tuple[list[dict], list[dict]]:
    """Report the six 2v2 allocations and column-level sensitivity diagnostics."""
    column_labels = list(ASPC1_COLUMNS)
    all_allocations = []
    for treated in itertools.combinations(range(4), 2):
        controls = [index for index in range(4) if index not in treated]
        effects = {
            gene: math.log2((statistics.mean(per_gene_column_cpm[gene][index] for index in controls) + 1.0) /
                            (statistics.mean(per_gene_column_cpm[gene][index] for index in treated) + 1.0))
            for gene in genes
        }
        ranks = midrank_percentiles(effects)
        all_allocations.append({
            "treated_columns": ",".join(column_labels[index] for index in treated),
            "is_observed_allocation": int(set(treated) == {1, 3}),
            "mean_gpi_drug_rank_percentile": statistics.mean(ranks[gene] for gene in gpi),
            "median_gpi_effect": statistics.median(effects[gene] for gene in gpi),
            "PIGF_effect": effects.get("PIGF", ""),
            "PIGF_rank_percentile": ranks.get("PIGF", ""),
            "ITGB1_effect": effects.get("ITGB1", ""),
            "ITGB1_rank_percentile": ranks.get("ITGB1", ""),
        })

    observed_control, observed_treatment = [0, 2], [1, 3]
    loo_rows = []
    for omit in range(4):
        controls = [index for index in observed_control if index != omit]
        treatments = [index for index in observed_treatment if index != omit]
        effects = {
            gene: math.log2((statistics.mean(per_gene_column_cpm[gene][index] for index in controls) + 1.0) /
                            (statistics.mean(per_gene_column_cpm[gene][index] for index in treatments) + 1.0))
            for gene in genes
        }
        ranks = midrank_percentiles(effects)
        loo_rows.append({
            "diagnostic": "leave_one_column_out",
            "omitted_column": column_labels[omit],
            "n_control_columns": len(controls),
            "n_treatment_columns": len(treatments),
            "mean_gpi_drug_rank_percentile": statistics.mean(ranks[gene] for gene in gpi),
            "median_gpi_effect": statistics.median(effects[gene] for gene in gpi),
            "PIGF_effect": effects.get("PIGF", ""),
            "PIGF_rank_percentile": ranks.get("PIGF", ""),
            "ITGB1_effect": effects.get("ITGB1", ""),
            "ITGB1_rank_percentile": ranks.get("ITGB1", ""),
        })
    for block, control_index, treatment_index in (("p11_pair", 0, 1), ("p11_A_pair", 2, 3)):
        effects = {
            gene: math.log2((per_gene_column_cpm[gene][control_index] + 1.0) /
                            (per_gene_column_cpm[gene][treatment_index] + 1.0))
            for gene in genes
        }
        ranks = midrank_percentiles(effects)
        loo_rows.append({
            "diagnostic": "single_paired_column_block",
            "omitted_column": block,
            "n_control_columns": 1,
            "n_treatment_columns": 1,
            "mean_gpi_drug_rank_percentile": statistics.mean(ranks[gene] for gene in gpi),
            "median_gpi_effect": statistics.median(effects[gene] for gene in gpi),
            "PIGF_effect": effects.get("PIGF", ""),
            "PIGF_rank_percentile": ranks.get("PIGF", ""),
            "ITGB1_effect": effects.get("ITGB1", ""),
            "ITGB1_rank_percentile": ranks.get("ITGB1", ""),
        })
    return all_allocations, loo_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--drug-counts", required=True, type=Path, help="GSE249537_counts.txt.gz")
    parser.add_argument("--car-table-s1", required=True, type=Path, help="Hagel et al. Table S1 workbook")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    drug_path = args.drug_counts.resolve()
    car_path = args.car_table_s1.resolve()
    out_dir = args.output_dir.resolve()
    if not drug_path.is_file() or not car_path.is_file():
        raise FileNotFoundError("Both public source files must exist; see the repository README for access details")
    out_dir.mkdir(parents=True, exist_ok=True)

    drug, drug_qc = read_aspc1_counts(drug_path)
    car, car_qc = read_car_workbook(car_path)
    eligible_drug = {gene for gene, cpm in drug["dmso_cpm"].items() if cpm >= 5.0}
    eligible_both = sorted(gene for gene in eligible_drug if gene in car["HuP-T3"] and gene in car["KP4"])
    eligible_any = sorted(gene for gene in eligible_drug if gene in car["HuP-T3"] or gene in car["KP4"])
    gene_sum_rank = midrank_percentiles({gene: drug["gene_sum_effect"][gene] for gene in eligible_both})
    guide_rank = midrank_percentiles({gene: drug["guide_summary"][gene]["median"] for gene in eligible_both})
    hup_rank = midrank_percentiles({gene: car["HuP-T3"][gene][0] for gene in eligible_both})
    kp4_rank = midrank_percentiles({gene: car["KP4"][gene][0] for gene in eligible_both})

    eligible_rows = []
    for gene in eligible_both:
        hup_effect, hup_neglog10p = car["HuP-T3"][gene]
        kp4_effect, kp4_neglog10p = car["KP4"][gene]
        guide = drug["guide_summary"][gene]
        eligible_rows.append({
            "gene": gene,
            "in_gpi": int(gene in GPI_GENES),
            "dmsocpm_mean_of_normalized_columns": drug["dmso_cpm"][gene],
            "gene_summed_drug_effect": drug["gene_sum_effect"][gene],
            "gene_summed_drug_rank_percentile": gene_sum_rank[gene],
            "per_guide_n": guide["n_guides"],
            "per_guide_median_effect": guide["median"],
            "per_guide_q1": guide["q1"],
            "per_guide_q3": guide["q3"],
            "per_guide_drug_rank_percentile": guide_rank[gene],
            "per_guide_positive_n": guide["positive_guides"],
            "HuP_T3_escape_lfc": hup_effect,
            "HuP_T3_source_neglog10p": hup_neglog10p,
            "HuP_T3_rank_percentile": hup_rank[gene],
            "KP4_escape_lfc": kp4_effect,
            "KP4_source_neglog10p": kp4_neglog10p,
            "KP4_rank_percentile": kp4_rank[gene],
        })
    write_tsv(out_dir / "eligible_gene_effects.tsv", eligible_rows)

    thresholds = []
    primary_hits: set[str] = set()
    baseline_all_methods: dict[str, set[str]] = {}
    score_methods = (
        ("gene_summed_primary", drug["gene_sum_effect"]),
        ("per_guide_median_robustness", {gene: item["median"] for gene, item in drug["guide_summary"].items()}),
    )
    for score_name, scores in score_methods:
        for drug_min, escape_min, neglog10p_min in THRESHOLDS:
            hits = set()
            for gene in eligible_any:
                if gene not in scores or scores[gene] < drug_min:
                    continue
                hup = car["HuP-T3"].get(gene)
                kp4 = car["KP4"].get(gene)
                hup_hit = bool(hup and hup[0] >= escape_min and hup[1] >= neglog10p_min)
                kp4_hit = bool(kp4 and kp4[0] >= escape_min and kp4[1] >= neglog10p_min)
                if hup_hit or kp4_hit:
                    hits.add(gene)
            baseline = (drug_min, escape_min, neglog10p_min) == BASELINE
            gpi_hits = sorted(hits & GPI_GENES)
            thresholds.append({
                "drug_score_method": score_name,
                "drug_min": drug_min,
                "escape_lfc_min": escape_min,
                "source_neglog10p_min": neglog10p_min,
                "all_dual_hit_n": len(hits),
                "all_dual_hit_genes": ",".join(sorted(hits)),
                "gpi_dual_hit_n": len(gpi_hits),
                "gpi_dual_hit_genes": ",".join(gpi_hits),
                "declared_baseline": int(baseline),
            })
            if baseline:
                baseline_all_methods[score_name] = hits
                if score_name == "gene_summed_primary":
                    primary_hits = hits
    write_tsv(out_dir / "threshold_sensitivity.tsv", thresholds)

    dual_rows = []
    for gene in sorted(primary_hits | GPI_GENES):
        hup = car["HuP-T3"].get(gene)
        kp4 = car["KP4"].get(gene)
        dual_rows.append({
            "gene": gene,
            "in_gpi": int(gene in GPI_GENES),
            "in_three_source_intersection": int(gene in eligible_both),
            "dmsocpm_ge5": int(gene in eligible_drug),
            "dmsocpm_mean": drug["dmso_cpm"].get(gene, ""),
            "gene_summed_drug_effect": drug["gene_sum_effect"].get(gene, ""),
            "per_guide_median_drug_effect": drug["guide_summary"].get(gene, {}).get("median", ""),
            "HuP_T3_escape_lfc": hup[0] if hup else "",
            "HuP_T3_source_neglog10p": hup[1] if hup else "",
            "KP4_escape_lfc": kp4[0] if kp4 else "",
            "KP4_source_neglog10p": kp4[1] if kp4 else "",
            "primary_baseline_dual_hit": int(gene in primary_hits),
            "gene_sum_and_guide_baseline_concordant": int(
                (gene in primary_hits) == (gene in baseline_all_methods.get("per_guide_median_robustness", set()))
            ),
        })
    write_tsv(out_dir / "primary_dual_gene_results.tsv", dual_rows)

    car_source_summary = []
    for screen, scores in car.items():
        hits = {gene for gene, (effect, neglog10p) in scores.items() if effect >= 0.5 and neglog10p >= 2.0}
        gpi_hits = hits & GPI_GENES
        p_value = hypergeom_upper_tail(len(scores), len(hits), len(GPI_GENES), len(gpi_hits))
        car_source_summary.append({
            "screen": screen,
            "scored_genes": len(scores),
            "source_cutoff_hits": len(hits),
            "gpi_panel_n": len(GPI_GENES),
            "gpi_covered_n": len(GPI_GENES & set(scores)),
            "gpi_source_cutoff_hits": len(gpi_hits),
            "one_sided_fisher_hypergeometric_p": p_value,
            "gpi_hit_genes": ",".join(sorted(gpi_hits)),
        })
    write_tsv(out_dir / "car_source_enrichment.tsv", car_source_summary)

    matching, pools = nearest_match_analysis(eligible_rows)
    (out_dir / "matched_rank_summary.json").write_text(json.dumps(matching, indent=2), encoding="utf-8")
    write_tsv(out_dir / "matched_pool_support.tsv", pools)
    conditional_rows = conditional_rank_matching(eligible_rows)
    write_tsv(out_dir / "conditional_drug_rank_matching.tsv", conditional_rows)
    gpi_covered = sorted(gene for gene in eligible_both if gene in GPI_GENES)
    allocation_rows, column_loo_rows = aspc1_column_diagnostics(
        drug["per_gene_column_cpm"], eligible_both, gpi_covered
    )
    write_tsv(out_dir / "aspc1_all_2v2_allocations.tsv", allocation_rows)
    write_tsv(out_dir / "aspc1_column_loo_qc.tsv", column_loo_rows)

    # Descriptive Spearman correlation of the two drug-score rankings.
    primary_rank = midrank_percentiles({gene: drug["gene_sum_effect"][gene] for gene in eligible_both})
    guide_rank_all = midrank_percentiles({gene: drug["guide_summary"][gene]["median"] for gene in eligible_both})
    x = [primary_rank[gene] for gene in eligible_both]
    y = [guide_rank_all[gene] for gene in eligible_both]
    x_mean, y_mean = statistics.mean(x), statistics.mean(y)
    denominator = math.sqrt(sum((value - x_mean) ** 2 for value in x) * sum((value - y_mean) ** 2 for value in y))
    spearman = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y)) / denominator

    qc = {
        "script": Path(__file__).name,
        "python_version": sys.version.split()[0],
        "inputs": {
            "drug_counts": {"filename": drug_path.name, "sha256": sha256(drug_path)},
            "car_table_s1": {"filename": car_path.name, "sha256": sha256(car_path)},
        },
        "aspc1": {
            **drug_qc,
            "eligible_at_mean_dmso_cpm_ge5": len(eligible_drug),
            "aligned_to_at_least_one_car_source": len(eligible_any),
            "aligned_to_both_car_sources": len(eligible_both),
            "gpi_panel_covered_in_three_source_intersection": len(GPI_GENES & set(eligible_both)),
            "gpi_panel_n": len(GPI_GENES),
            "primary_baseline_all_dual_hit_n": len(primary_hits),
            "primary_baseline_all_dual_hit_genes": sorted(primary_hits),
            "primary_baseline_gpi_hit_n": len(primary_hits & GPI_GENES),
            "primary_baseline_gpi_hit_genes": sorted(primary_hits & GPI_GENES),
            "frozen_six_gpi_hit_genes_retained_after_cpm_filter": sorted(FROZEN_DUAL_GPI & eligible_drug),
            "gene_sum_vs_per_guide_spearman_rank_correlation_both_car_universe": spearman,
        },
        "car_sources": car_qc,
        "threshold_grid": {
            "threshold_combinations_per_drug_score": len(THRESHOLDS),
            "drug_thresholds": [0.3, 0.5, 0.75, 1.0],
            "escape_lfc_thresholds": [0.3, 0.5],
            "source_neglog10p_thresholds": [1.3, 2.0],
            "car_sources_combined_by_or": True,
        },
        "matching": matching,
        "conditional_rank_matching": {
            "designs": ["AsPC1_DMSO_CPM_decile", "AsPC1_DMSO_CPM_decile_plus_HuPT3_and_KP4_escape_rank_quintiles"],
            "draws_per_design": MATCH_DRAWS,
            "seed_by_design": {"AsPC1_DMSO_CPM_decile": 261123,
                               "AsPC1_DMSO_CPM_decile_plus_HuPT3_and_KP4_escape_rank_quintiles": 261124},
            "rows": conditional_rows,
            "interpretation": "Descriptive conditional summaries; do not provide population-level inference or adjust for post-hoc pathway selection or gene dependence.",
        },
        "limits": [
            "This is a reanalysis of existing public screen data, not an independent experiment.",
            "The AsPC-1 comparison has two treatment and two control profiles; their biological independence is not established by the source metadata.",
            "Guide-level estimates are nested within genes and do not create biological replicate counts.",
            "The nearest-neighbour rank comparison is descriptive and does not correct post-hoc pathway selection or gene dependence.",
            "The source CAR -log10(p) values are carried from the publication and are not recomputed here.",
        ],
    }
    (out_dir / "qc_summary.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(json.dumps({
        "aspc1_eligible": len(eligible_drug),
        "three_source_intersection": len(eligible_both),
        "gpi_covered": len(GPI_GENES & set(eligible_both)),
        "gpi_primary_baseline_hits": sorted(primary_hits & GPI_GENES),
        "matching_metrics": matching["metrics"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
