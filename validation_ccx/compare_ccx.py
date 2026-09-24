#!/usr/bin/env python3
"""Compare a CalculiX results table against the stored Abaqus one, row by row.

    python3 compare_ccx.py <abaqus.csv> <calculix.csv> [col ...]

Reports the relative difference per run_id and column. Plain csv/math, so it
runs without pandas.
"""
import csv
import sys

DEFAULT_COLS = ['E_eff', 'nu_eff', 'G_eff', 'porosity',
                'D_rve', 'D_classical', 'D_ratio', 'l', 'E_bending']

# Below this, a value is floating-point zero rather than a small number, and a
# relative difference between two of them is noise reported as a percentage.
# The homogeneous cube's porosity lands at ~1e-8 in one solver and ~1e-9 in the
# other; calling that "110% disagreement" would be the single most misleading
# line in the table.
NEGLIGIBLE = 1e-6


def load(path):
    with open(path) as f:
        return dict((r['run_id'], r) for r in csv.DictReader(f))


def num(row, col):
    v = row.get(col, '')
    if v in ('', 'ERROR', 'MISSING', None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def main(argv):
    ref, new = load(argv[1]), load(argv[2])
    cols = argv[3:] or DEFAULT_COLS
    shared = [k for k in ref if k in new]
    if not shared:
        print('no run_id in common')
        return 1

    print('{:<14} {:<13} {:>18} {:>18} {:>12}'.format(
        'run_id', 'column', 'abaqus', 'calculix', 'rel.diff'))
    print('-' * 80)
    worst = {}
    for rid in shared:
        for c in cols:
            a, b = num(ref[rid], c), num(new[rid], c)
            if a is None or b is None:
                continue
            scale = max(abs(a), abs(b))
            if scale < NEGLIGIBLE:
                print('{:<14} {:<13} {:>18.10g} {:>18.10g} {:>12}'.format(
                    rid, c, a, b, 'both ~0'))
                continue
            rel = (b - a) / scale
            print('{:<14} {:<13} {:>18.10g} {:>18.10g} {:>11.3g}%'.format(
                rid, c, a, b, 100 * rel))
            if abs(rel) > abs(worst.get(c, (0, ''))[0]):
                worst[c] = (rel, rid)
        print('-' * 80)

    print('\nworst relative difference per column')
    for c in cols:
        if c in worst:
            rel, rid = worst[c]
            print('  {:<14} {:>10.3g}%   ({})'.format(c, 100 * rel, rid))

    # Columns the CalculiX side reports and the Abaqus side has no counterpart
    # for are skipped by the loop above, because there is nothing to difference.
    # equilibrium_gap is the important one: it is the evidence that the solve
    # converged, and on a large cell it is the ONLY such evidence, so print it
    # rather than let it fall silently out of the report.
    for c in ('equilibrium_gap', 'solver'):
        vals = [(rid, new[rid].get(c, '')) for rid in shared
                if new[rid].get(c, '') not in ('', None)]
        if vals and not any(ref[rid].get(c, '') for rid, _ in vals):
            print('\n{} (CalculiX only)'.format(c))
            for rid, v in vals:
                try:
                    print('  {:<14} {:>12.3e}'.format(rid, float(v)))
                except (TypeError, ValueError):
                    print('  {:<14} {:>12}'.format(rid, v))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
