#!/usr/bin/env python3
"""Check pooled-library column and denominator provenance without modifying upstream files."""

import csv
import gzip
import json
import math
import sys
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
from pathlib import Path
from statistics import mean, median

from openpyxl import load_workbook

root = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
out.mkdir(parents=True, exist_ok=True)
panel = set("PIGA PIGB PIGC PIGF PIGG PIGH PIGK PIGL PIGM PIGN PIGO PIGP PIGQ PIGS PIGT PIGU PIGV PIGW PIGY GPAA1 DPM1 DPM2 DPM3 MPDU1 PGAP1 PGAP2 PGAP3".split())
book = load_workbook(root / "01_data/frose_primary/AllPools_assembled.xlsx", read_only=True, data_only=True)

def usable(number):
    return isinstance(number, (int, float)) and math.isfinite(number)

records = defaultdict(dict)
flags = Counter()
inconsistent = []
negative_positive_lfc_diff = Counter()
rows = book.active.values
header = next(rows)
for row in rows:
    record = dict(zip(header, row))
    gene = str(record["gene"] or "").strip().upper()
    pool = str(record["Pool"] or "").strip()
    tag = record["OnePerGene"]
    if tag is not None:
        flags[str(tag)] += 1
    neg_lfc = record["neg|lfc"]
    pos_lfc = record["pos|lfc"]
    fdr = record["pos|fdr"]
    if not gene or not pool or not usable(neg_lfc) or not usable(pos_lfc) or not usable(fdr) or not 0 <= fdr <= 1:
        continue
    if abs(neg_lfc - pos_lfc) > 1e-8:
        negative_positive_lfc_diff[pool] += 1
    record = (neg_lfc, pos_lfc, fdr)
    previous = records[pool].get(gene)
    if previous is not None and previous != record:
        inconsistent.append((pool, gene, previous, record))
    else:
        records[pool][gene] = record
book.close()

pool_summary = []
for pool, table in sorted(records.items()):
    gpi = sorted(table.keys() & panel)
    row = {"pool": pool, "genes": len(table), "gpi": len(gpi), "gpi_members": gpi,
           "all_pos_hits": sum(value[1] >= 0.5 and value[2] <= 0.1 for value in table.values()),
           "all_neg_lfc_hits": sum(value[0] >= 0.5 and value[2] <= 0.1 for value in table.values()),
           "gpi_pos_hits": [gene for gene in gpi if table[gene][1] >= 0.5 and table[gene][2] <= 0.1],
           "gpi_neg_lfc_hits": [gene for gene in gpi if table[gene][0] >= 0.5 and table[gene][2] <= 0.1],
           "median_gpi_pos_lfc": median(table[gene][1] for gene in gpi) if gpi else None,
           "cd19": table.get("CD19"), "gpi_min_and_max_positive_lfc": (
             [min(table[gene][1] for gene in gpi), max(table[gene][1] for gene in gpi)] if gpi else None),
    }
    pool_summary.append(row)
    print(json.dumps(row), flush=True)

both = defaultdict(list)
for pool, table in records.items():
    for gene, value in table.items():
        both[gene].append((pool, value))
shared = {gene: rows for gene, rows in both.items() if len(rows) > 1}
with gzip.open(root / "01_data/cart_human_gse180834/GSE180834_gene_summary_allreps.txt.gz", "rt", newline="") as handle:
    tcr = list(csv.DictReader(handle, delimiter="\t"))
summary = {
    "pool_summary": pool_summary, "unique_across_pools": len(both),
    "pool_gene_pairs": sum(len(table) for table in records.values()),
    "multi_pool_gene_count": len(shared), "gpi_multi_pool_genes": sorted(shared.keys() & panel),
    "gpi_across_pools": sorted(both.keys() & panel), "divergent_lfc_rows": dict(negative_positive_lfc_diff),
    "nonidentical_within_pool_gene_records": len(inconsistent),
    "nonidentical_within_pool_examples": inconsistent[:12],
    "onepergene_nonempty": sum(flags.values()), "onepergene_unique_values": len(flags),
    "onepergene_example_counts": flags.most_common(8),
    "tcr_rows": len(tcr), "tcr_nonempty_unique_gene_ids": len({r["id"].upper().strip() for r in tcr if r["id"].strip()}),
}
(out / "frose_pool_qc.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps({key: val for key, val in summary.items() if key != "pool_summary"}), flush=True)
gpi_pools = [row["pool"] for row in pool_summary if row["gpi"]]
assert gpi_pools == ["p1112"], gpi_pools
valid = records[gpi_pools[0]]
effects = {gene: row[1] for gene, row in valid.items()}
ordered = sorted(effects.values())
hits = {gene for gene, record in valid.items() if record[1] >= 0.5 and record[2] <= 0.1}
targets = set("PIGC PIGF PIGM PIGP PIGS PIGT".split())
panel_rows = [
    {
        "screen": "Fröse in vivo (p1112)", "gene": gene, "effect": effects[gene],
        "rank_pct": 50 * (bisect_left(ordered, effects[gene]) + bisect_right(ordered, effects[gene])) / len(ordered),
        "meets_screen_hit_cutoff": int(gene in hits), "frozen_dual_hit": int(gene in targets),
    }
    for gene in sorted(valid.keys() & panel)
]
with (out / "frose_p1112_gene_effects.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=panel_rows[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(panel_rows)
print("Pool p1112 median and mean rank percentiles", median(row["rank_pct"] for row in panel_rows),
      mean(row["rank_pct"] for row in panel_rows), flush=True)
print("DONE", flush=True)

