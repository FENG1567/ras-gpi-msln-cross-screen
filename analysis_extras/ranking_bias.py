"""Reproduce the same-readout ranking bias illustration in the supplement."""

import argparse
import csv
import math
from pathlib import Path
import random

N_GENES = 2000
TOP_K = 200
REPETITIONS = 200
SEED = 7


def mean(values):
    return sum(values) / len(values)


def standard_deviation(values):
    center = mean(values)
    return math.sqrt(sum((value - center) ** 2 for value in values) / (len(values) - 1))


def z_scores(values):
    center = mean(values)
    scale = standard_deviation(values) or 1.0
    return [(value - center) / scale for value in values]


def run_one(rng, penalty, correlation):
    drug = [rng.gauss(0, 1) for _ in range(N_GENES)]
    escape = [correlation * value + math.sqrt(1 - correlation**2) * rng.gauss(0, 1) for value in drug]
    drug_rank_score = z_scores(drug)
    escape_rank_score = z_scores(escape)
    ordered_escape = sorted(escape)
    percentile_index = (N_GENES - 1) * .75
    lower = math.floor(percentile_index)
    threshold = ordered_escape[lower] + (ordered_escape[lower + 1] - ordered_escape[lower]) * (percentile_index - lower)
    hits = [value >= threshold for value in escape]
    plain = sorted(range(N_GENES), key=lambda i: -drug_rank_score[i])[:TOP_K]
    penalised_score = [drug_rank_score[i] - penalty * max(0.0, escape_rank_score[i]) for i in range(N_GENES)]
    penalised = sorted(range(N_GENES), key=lambda i: -penalised_score[i])[:TOP_K]
    n_plain = sum(hits[i] for i in plain)
    n_penalised = sum(hits[i] for i in penalised)
    return 100 * (n_plain - n_penalised) / n_plain if n_plain else 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('generated/ranking_bias.csv'))
    args = parser.parse_args()
    rng = random.Random(SEED)
    rows = []
    for correlation in (-0.5, -0.25, 0.0, 0.25, 0.5):
        for penalty in (0.5, 1.0, 2.0):
            contrasts = [run_one(rng, penalty, correlation) for _ in range(REPETITIONS)]
            rows.append({'rho': correlation, 'lambda': penalty, 'delta_mean_percent': round(mean(contrasts), 2),
                         'delta_sd_percent': round(standard_deviation(contrasts), 2)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    example = next(r for r in rows if r['rho'] == 0 and r['lambda'] == 1)
    print(f"Same-readout null simulation: {example['delta_mean_percent']:.2f}% relative reduction.")


if __name__ == '__main__':
    main()
