#!/usr/bin/env python3
"""Report the manufactured-solution sweep collected in mms.csv.

    python3 elements_ccx/tests/mms_report.py out_mms/mms.csv

Prints one block per (layout, kg, arm) with the errors per mesh and the
observed rate over the two finest meshes, then the K/G summary at the
finest mesh of each layout.
"""
from __future__ import annotations

import csv
import math
import sys
from typing import Dict, List, Tuple

Row = Tuple[int, str, float, str, float, float, float, float, float]


def load(path: str) -> List[Row]:
    """Read the sweep CSV.

    Args:
        path: Path of mms.csv.

    Returns:
        The rows as tuples.
    """
    rows: List[Row] = []
    with open(path) as handle:
        for rec in csv.DictReader(handle):
            rows.append((int(rec['n']), rec['layout'], float(rec['kg']),
                         rec['arm'], float(rec['h']), float(rec['l2']),
                         float(rec['h1']), float(rec['l2_incl']),
                         float(rec['h1_incl'])))
    return rows


def rate(a: float, b: float, na: int, nb: int) -> float:
    """Observed convergence rate between two meshes.

    Args:
        a: Error on the coarser mesh.
        b: Error on the finer mesh.
        na: Coarse mesh size.
        nb: Fine mesh size.

    Returns:
        log(a/b) / log(nb/na).
    """
    return float('nan') if b <= 0 else math.log(a / b) / math.log(nb / na)


def main() -> int:
    """Print the tables.

    Returns:
        Process exit status.
    """
    rows = load(sys.argv[1])
    groups: Dict[Tuple[str, float, str], List[Row]] = {}
    for r in rows:
        groups.setdefault((r[1], r[2], r[3]), []).append(r)
    finest: Dict[str, int] = {}
    for layout in ('single', 'two'):
        ns = [r[0] for r in rows if r[1] == layout]
        if ns:
            finest[layout] = max(ns)

    for (layout, kg, arm), rs in sorted(groups.items()):
        rs.sort()
        print('\n%s  K/G=%g  %s' % (layout, kg, arm))
        for n, _, _, _, h, l2, h1, l2i, h1i in rs:
            print('  n=%-3d  h %.4f  L2 %.4e  H1 %.4e' % (n, h, l2, h1))
        if len(rs) >= 2:
            rate_l2 = rate(rs[-2][5], rs[-1][5], rs[-2][0], rs[-1][0])
            rate_h1 = rate(rs[-2][6], rs[-1][6], rs[-2][0], rs[-1][0])
            print('  last-pair rates: L2 %.2f  H1 %.2f' % (rate_l2, rate_h1))

    for layout in ('single', 'two'):
        if layout not in finest:
            continue
        print('\nK/G summary at n=%d (%s phase)'
              % (finest[layout], layout))
        print('  %-10s %-14s %10s %10s %10s' % ('arm', 'K/G', 'L2', 'H1', 'incl H1'))
        for kg in sorted({r[2] for r in rows if r[1] == layout}):
            for arm in sorted({r[3] for r in rows
                               if r[1] == layout and r[2] == kg}):
                sel = [r for r in rows if r[1] == layout and r[2] == kg
                       and r[3] == arm and r[0] == finest[layout]]
                if not sel:
                    continue
                _, _, _, _, _, l2, h1, l2i, h1i = sel[0]
                print('  %-10s %-14g %10.3e %10.3e %10.3e'
                      % (arm, kg, l2, h1, h1i))
    return 0


if __name__ == '__main__':
    sys.exit(main())
