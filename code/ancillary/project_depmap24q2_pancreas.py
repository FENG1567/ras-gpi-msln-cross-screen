#!/usr/bin/env python3
"""Stream the DepMap Public 24Q2 gene-effect matrix to pancreas models only.

The selected matrix rows are written without numerical transformation. Run
`depmap24q2_sensitivity.py` on this projection to reproduce the summarized
pancreas and PAAD-only gene-effect distributions.
"""

import argparse
import csv
import gzip
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path, help="Official CRISPRGeneEffect.csv")
    parser.add_argument("model", type=Path, help="Official Model.csv from the same 24Q2 release")
    parser.add_argument("output", type=Path, help="Output .csv.gz path for pancreas-annotated model rows")
    args = parser.parse_args()

    with args.model.open(newline="", encoding="utf-8-sig") as handle:
        metadata = {
            row["ModelID"]: row
            for row in csv.DictReader(handle)
            if "Pancrea" in (row.get("OncotreePrimaryDisease") or "")
        }
    if not metadata:
        raise ValueError("No pancreas-annotated models found in Model.csv")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    retained = set()
    with args.matrix.open(newline="", encoding="utf-8-sig") as source, gzip.open(
        args.output, "wt", newline="", encoding="utf-8"
    ) as target:
        reader = csv.reader(source)
        writer = csv.writer(target, lineterminator="\n")
        header = next(reader)
        if not header or header[0] != "":
            raise ValueError("Unexpected gene-effect matrix header")
        if len(set(header[1:])) != len(header[1:]):
            raise ValueError("Duplicate gene columns in the matrix header")
        writer.writerow(header)
        for row in reader:
            if len(row) != len(header):
                raise ValueError("Unexpected row width in the gene-effect matrix")
            model_id = row[0]
            if not model_id or model_id in seen:
                raise ValueError(f"Missing or duplicate matrix ModelID: {model_id!r}")
            seen.add(model_id)
            if model_id in metadata:
                writer.writerow(row)
                retained.add(model_id)

    if not retained:
        raise ValueError("No pancreas-annotated matrix rows matched Model.csv")
    print(
        f"matrix_model_rows={len(seen)} metadata_pancreas_models={len(metadata)} "
        f"retained_rows={len(retained)} gene_columns={len(header) - 1}"
    )
    if len(retained) != len(metadata):
        missing = sorted(set(metadata) - retained)
        print(f"metadata_models_without_gene_effect_rows={len(missing)}")


if __name__ == "__main__":
    main()
