"""Check released result tables against manuscript anchors."""

import argparse
import csv
from math import atanh, comb, isclose, sqrt, tanh
from pathlib import Path


def read_table(path, delimiter='\t'):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def fisher_tail(n_universe, screen_hits, panel_size, k):
    denominator = comb(n_universe, panel_size)
    return sum(comb(screen_hits, i) * comb(n_universe-screen_hits, panel_size-i)
               for i in range(k, min(screen_hits, panel_size)+1)) / denominator


def check(dataset):
    screen = {r['label']: r for r in read_table(dataset/'immune_screen_summary.tsv')}
    for name, reference in [('HuP-T3', 1.210095322e-31), ('KP4_MSLN', 2.075860734e-42)]:
        row = screen[name]
        tail = fisher_tail(*(int(row[f]) for f in ('n', 'background_hits', 'gpi', 'hits')))
        assert isclose(tail, reference, rel_tol=1e-7), (name, tail)
    row = screen['Fröse p1112']
    assert [int(row[f]) for f in ('n', 'background_hits', 'gpi', 'hits')] == [923, 9, 20, 0]
    grid = read_table(dataset/'primary_threshold_grid.tsv')
    assert len(grid) == 16 and {int(r['hits']) for r in grid} == {3,5,7,8,11,13}
    baseline = next(r for r in grid if tuple(map(float, (r['dmin'],r['emin'],r['pmin']))) == (.5,.5,2))
    assert int(baseline['hits']) == 7
    matched = {r['name']: r for r in read_table(dataset/'matched_rank_summary.tsv')}
    for label, median in [('CPM + both CAR ranks', 46.8), ('CPM + CAR neighbours', 46.6)]:
        row = matched[label]
        assert round(float(row['observed'])*100,1) == 72.8
        assert round(float(row['median'])*100,1) == median
        assert [int(row[f]) for f in ('reached','draws','n')] == [0,20000,25]
    kegg = [r for r in read_table(dataset/'KEGG_hsa00563_coverage.tsv')]
    assert (len(kegg),sum(int(r['aligned']) for r in kegg),sum(int(r['overlap_hit']) for r in kegg)) == (30,26,7)
    for r in read_table(dataset/'prism_pearson_ci.tsv'):
        if not r['association'].startswith('GPI '):
            continue
        z = atanh(float(r['r']))
        margin = 1.96/sqrt(int(r['n'])-3)
        assert abs(tanh(z-margin)-float(r['ci95_low'])) < .001
        assert abs(tanh(z+margin)-float(r['ci95_high'])) < .001
    print('Released immune counts, Fisher tails, GPI thresholds, matching, pathway coverage and PRISM intervals verified.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path(__file__).resolve().parents[1]/'data'/'derived')
    args = parser.parse_args()
    check(args.data)
