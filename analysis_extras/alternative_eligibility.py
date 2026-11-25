"""Reproduce the separate AsPC-1/UM53 coverage and held-out ranking checks."""

import argparse
from bisect import bisect_left
from collections import defaultdict
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import random
import statistics
from xml.etree import ElementTree
import zipfile

UM53_CONTROLS = ('UM53-TKO-DMSO-P8', 'UM53-TKO-DMSO-T2-P8', 'UM53-TKO-DMSO-T3-P8')
ASPC_CONTROLS = ('ASPC1-DMSO-p11','ASPC1-DMSO-p11-A')
ASPC_TREATED = ('ASPC1-MRTX-6-p11','ASPC1-MRTX-6-p11-A')
PANEL = set('PIGA PIGB PIGC PIGF PIGG PIGH PIGK PIGL PIGM PIGN PIGO PIGP PIGQ PIGS PIGT PIGU PIGV PIGW PIGY GPAA1 DPM1 DPM2 DPM3 MPDU1 PGAP1 PGAP2 PGAP3'.split())
GRID = ((0.3, 0.5, 0.75, 1.0), (0.3, 0.5), (1.3, 2.0))


def quantile(values, probability):
    ordered = sorted(values)
    pos = probability * (len(ordered) - 1)
    left = math.floor(pos)
    return ordered[left] + (pos-left) * (ordered[min(left+1,len(ordered)-1)]-ordered[left])


def read_drug(counts):
    columns = ASPC_CONTROLS + ASPC_TREATED + UM53_CONTROLS
    totals = [0] * len(columns)
    gene_totals = defaultdict(lambda: [0]*len(columns))
    with gzip.open(counts, 'rt', encoding='utf-8-sig', newline='') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader)
        indexes = [header.index(column) for column in columns]
        for row in reader:
            try:
                nums = [float(row[index]) for index in indexes]
            except (IndexError, ValueError):
                continue
            gene = row[1]
            for j, n in enumerate(nums):
                totals[j] += n
                gene_totals[gene][j] += n
    result = {}
    for gene, values in gene_totals.items():
        cpm = [n/total*1e6 for n,total in zip(values,totals)]
        control, drug, um_control = (statistics.mean(x) for x in (cpm[:2],cpm[2:4],cpm[4:]))
        if control >= 5 and um_control >= 5:
            result[gene] = (control, -math.log2((drug+1)/(control+1)))
    return result


def read_car(path):
    ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
    with zipfile.ZipFile(path) as book:
        root = ElementTree.fromstring(book.read('xl/sharedStrings.xml'))
        strings = [''.join(t.text or '' for t in entry.iter(ns+'t')) for entry in root.findall(ns+'si')]
        output = []
        for n in (1,2):
            sheet = ElementTree.fromstring(book.read(f'xl/worksheets/sheet{n}.xml'))
            score = {}
            for row in list(sheet.iter(ns+'row'))[1:]:
                values = []
                for cell in row.findall(ns+'c')[:3]:
                    value = cell.findtext(ns+'v')
                    values.append(strings[int(value)] if cell.get('t') == 's' and value is not None else value)
                if len(values) != 3 or not values[0]:
                    continue
                try:
                    score[values[0].strip()] = float(values[1]),float(values[2])
                except (ValueError, TypeError):
                    continue
            output.append(score)
    return output


def legacy_grid(rows):
    cpm = [r['cpm'] for r in rows]
    boundaries = [quantile(cpm, 0.1*j) for j in range(1,10)]
    strata = defaultdict(list)
    panel = []
    for r in rows:
        r['stratum'] = bisect_left(boundaries, r['cpm'])
        strata[r['stratum']].append(r)
        if r['gpi']:
            panel.append(r)
    rng = random.Random(20261121)
    output = []
    for d, e, p in itertools.product(*GRID):
        def meets(r):
            return r['drug'] >= d and ((r['hup'] >= e and r['hscore'] >= p) or (r['kp4'] >= e and r['kscore'] >= p))
        hit_names = [r['gene'] for r in panel if meets(r)]
        hit_pools = {key:[meets(r) for r in pool] for key,pool in strata.items()}
        n_reached = sum(sum(rng.choice(hit_pools[r['stratum']]) for r in panel) >= len(hit_names) for _ in range(2000))
        output.append({'drug_min':d,'escape_min':e,'source_score_min':p,'observed':len(hit_names),
                       'genes':hit_names,'draws_reached':n_reached,'draws':2000,'plus_one_tail':(n_reached+1)/2001})
    return len(panel), output


def standardize(values):
    center = statistics.mean(values)
    scale = statistics.stdev(values)
    return [(v-center)/scale for v in values]


def rank_rule(rows):
    drug = standardize([r['drug'] for r in rows])
    base = sorted(range(len(drug)), key=lambda i: -drug[i])
    output = []
    for build, assess, candidates in [('hup','kp4',(100,200,500)),('kp4','hup',(200,))]:
        build_z = standardize([r[build] for r in rows])
        held_out = [r[assess] for r in rows]
        threshold = quantile(held_out,.75)
        penalized = sorted(range(len(drug)), key=lambda i: -(drug[i]-max(0.,build_z[i])))
        for k in candidates:
            i_base, i_pen = base[:k], penalized[:k]

            def counts(values):
                return (sum(values[i]>=threshold for i in i_base), sum(values[i]>=threshold for i in i_pen))

            first, second = counts(held_out)
            delta = (first-second)/first
            p = None
            if k == 200 and build == 'hup':
                rng, shuffled = random.Random(11), held_out[:]
                reached = 0
                for _ in range(200):
                    rng.shuffle(shuffled)
                    plain, penalized_count = counts(shuffled)
                    null_reduction = (plain-penalized_count)/plain if plain else 0
                    reached += null_reduction >= delta
                p = (reached+1)/201
            output.append({'build':build,'evaluate':assess,'K':k,'R1_count':first,'R2_count':second,
                           'relative_reduction':delta,'R1_drug_mean':statistics.mean(drug[i] for i in i_base),
                           'R2_drug_mean':statistics.mean(drug[i] for i in i_pen),'plus_one_p':p})
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--drug-counts', required=True, type=Path)
    parser.add_argument('--car-table-s1', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=Path('generated/alternative_eligibility.json'))
    args = parser.parse_args()
    drug = read_drug(args.drug_counts)
    hup, kp4 = read_car(args.car_table_s1)
    # The alternative calculation retains the original case-sensitive gene-symbol alignment.
    rows = [{'gene':g, 'cpm':drug[g][0], 'drug':drug[g][1], 'gpi':g in PANEL,
             'hup':hup[g][0], 'hscore':hup[g][1], 'kp4':kp4[g][0], 'kscore':kp4[g][1]}
            for g in sorted(drug.keys() & hup.keys() & kp4.keys())]
    n_panel, grid = legacy_grid(rows)
    result = {'aligned_genes':len(rows),'covered_GPI':n_panel,'threshold_grid':grid,'ranking':rank_rule(rows)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(f"Alternative background {len(rows):,} genes, panel {n_panel}/27. Main held-out result {result['ranking'][1]}.")


if __name__ == '__main__':
    main()
