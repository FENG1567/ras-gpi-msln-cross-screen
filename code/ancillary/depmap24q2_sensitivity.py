"""Summarize the full-gene DepMap 24Q2 matrix in pancreas-annotated models."""

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np


GPI = set(
    "PIGA PIGB PIGC PIGF PIGG PIGH PIGK PIGL PIGM PIGN PIGO PIGP PIGQ "
    "PIGS PIGT PIGU PIGV PIGW PIGY GPAA1 DPM1 DPM2 DPM3 MPDU1 PGAP1 PGAP2 PGAP3".split()
)
DUAL = {"PIGC", "PIGF", "PIGL", "PIGM", "PIGP", "PIGS", "PIGT"}


def summarize(values):
    values = values[np.isfinite(values)]
    if not len(values):
        raise ValueError("No measured gene scores in a required comparison.")
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return {"n_genes": len(values), "median": median, "q1": q1, "q3": q3}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path, help="Lossless pancreas model-row projection")
    parser.add_argument("model", type=Path, help="Official 24Q2 Model.csv")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.model.open(newline="") as handle:
        models = {
            row["ModelID"]: row
            for row in csv.DictReader(handle)
            if "Pancrea" in row["OncotreePrimaryDisease"]
        }
    with gzip.open(args.matrix, "rt", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header[0] != "":
            raise ValueError("Unexpected matrix identifier column.")
        names = np.array([re.sub(r"\s+\(\d+\)$", "", column) for column in header[1:]])
        model_ids = []
        data = []
        for row in reader:
            if len(row) != len(header) or row[0] not in models:
                raise ValueError("Unexpected model identifier or row width.")
            model_ids.append(row[0])
            data.append([float(value) if value else np.nan for value in row[1:]])
    if len(set(model_ids)) != len(model_ids) or len(set(names)) != len(names):
        raise ValueError("Duplicate matrix axis labels.")
    matrix = np.asarray(data)
    gpi_mask = np.isin(names, list(GPI))
    dual_mask = np.isin(names, list(DUAL))
    ribosome_mask = np.array([name.startswith(("RPL", "RPS")) for name in names])
    if gpi_mask.sum() != len(GPI) or dual_mask.sum() != len(DUAL):
        raise ValueError("Incomplete locked GPI or overlap panel in the full matrix.")

    selected_rows = {
        "pancreas_annotated": np.ones(len(model_ids), dtype=bool),
        "paad_annotated": np.array([models[model]["OncotreeCode"] == "PAAD" for model in model_ids]),
    }
    group_summary = {}
    medians = {}
    measured_counts = {}
    for label, row_mask in selected_rows.items():
        cohort = matrix[row_mask]
        finite = np.isfinite(cohort)
        measured = finite.sum(axis=0)
        gene_medians = np.full(len(names), np.nan)
        gene_medians[measured > 0] = np.nanmedian(cohort[:, measured > 0], axis=0)
        medians[label] = gene_medians
        measured_counts[label] = measured
        group_summary[label] = {
            "n_models": int(row_mask.sum()),
            "oncotree_codes": dict(
                Counter(models[model]["OncotreeCode"] for model, keep in zip(model_ids, row_mask) if keep)
            ),
            "all_scored_genes": summarize(gene_medians),
            "non_gpi_genes": summarize(gene_medians[~gpi_mask]),
            "gpi": summarize(gene_medians[gpi_mask]),
            "gpi_dual_hits": summarize(gene_medians[dual_mask]),
            "ribosomal_prefix_rpl_rps": summarize(gene_medians[ribosome_mask]),
        }
    payload = {
        "release": "DepMap Public 24Q2",
        "doi": "10.25452/figshare.plus.25880521.v1",
        "score": "Chronos gene effect across model-level screens",
        "n_gene_columns": len(names),
        "n_pancreas_metadata": len(models),
        "models_with_chronos_omitted_by_paad_filter": int(
            len(model_ids) - selected_rows["paad_annotated"].sum()
        ),
        "cohorts": group_summary,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "depmap24q2_sensitivity.json").write_text(json.dumps(payload, indent=2) + "\n")
    with (args.output / "depmap24q2_gene_medians.tsv").open("w", newline="") as handle:
        fields = ["gene", "gpi", "dual_hit", "ribosomal_prefix"]
        for label in selected_rows:
            fields += [label + "_median", label + "_n_models"]
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        for index, name in enumerate(names):
            row = [name, int(gpi_mask[index]), int(dual_mask[index]), int(ribosome_mask[index])]
            for label in selected_rows:
                row += [medians[label][index], measured_counts[label][index]]
            writer.writerow(row)
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()

