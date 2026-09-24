#!/usr/bin/env python3
"""Mesh convergence for F-barES-FEM-T4 on a cell whose soft slab is resolved.

Per mesh size the driver generates the cell with make_slabconv.py, converts
the undrained deck to F-barES-FEM-T4 with fbares.py, solves three arms with
the given CalculiX binary (plain C3D4 drained and undrained, F-bar c = 1
undrained), and reports R = C1111(und)/C1111(drn) with report_slabconv.py.
The drained denominator is always plain C3D4, so the arms differ only in the
element under test.

    python3 meshconv.py --ccx /path/to/ccx_fbar
    python3 meshconv.py --ccx /path/to/ccx_fbar --ns 10 20 30 40 --kg 500
    python3 meshconv.py --dry-run

The same cell is rebuilt at every n.  The drained twin is the same mesh with
one elastic card rewritten, so R carries no packing noise: the only thing
that changes between points is h.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from typing import List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
FBARES = os.path.normpath(os.path.join(HERE, os.pardir, 'fbares.py'))
ARMS = ('und/m_ccx', 'drn/m_ccx', 'und_fbar1/m_ccx')


def pmake(script: str, args: List[str]) -> List[str]:
    """Python command for one helper script.

    Args:
        script: Script name inside this directory.
        args: Script arguments.

    Returns:
        The full command.
    """
    return [sys.executable, os.path.join(HERE, script)] + args


def solve(ccx: str, jobdir: str, threads: int) -> int:
    """Run one CalculiX arm, capturing the log.

    Args:
        ccx: CalculiX binary.
        jobdir: Directory holding ``m_ccx.inp``.
        threads: Value for OMP_NUM_THREADS.

    Returns:
        The solver exit status.
    """
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(os.path.join(jobdir, 'solve.log'), 'w') as log:
        return subprocess.call([ccx, 'm_ccx'], cwd=jobdir, env=env,
                               stdout=log, stderr=subprocess.STDOUT)


def main() -> int:
    """Run the sweep and the report.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ccx', default=os.environ.get('CCX_FBAR', 'ccx_fbar'),
                    help='patched CalculiX binary')
    ap.add_argument('--ns', type=int, nargs='+', default=[10, 20, 30, 40])
    ap.add_argument('--kg', type=float, default=500.0)
    ap.add_argument('--jitter', type=float, default=0.3)
    ap.add_argument('--bridge', default='one')
    ap.add_argument('--load', choices=['x', 'y', 'z'], default='x')
    ap.add_argument('--conf', default='sym')
    ap.add_argument('--root', default=None,
                    help='sweep directory (default out_meshconv/kg<kg>_<...>)')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--dry-run', action='store_true',
                    help='generate decks and print the solve commands')
    a = ap.parse_args()

    root = a.root or 'out_meshconv/kg%g_%s_%s%s' % (a.kg, a.bridge, a.load,
                                                    a.conf)
    ns = [str(n) for n in a.ns]

    for n in ns:
        for state in ('und', 'drn'):
            d = os.path.join(root, 'n' + n, state)
            os.makedirs(d, exist_ok=True)
            subprocess.run(pmake('make_slabconv.py',
                                 [os.path.join(d, 'm'), n, state, '%g' % a.kg,
                                  '%g' % a.jitter, a.bridge, a.load, a.conf]),
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT)
        w = os.path.join(root, 'n' + n, 'und_fbar1')
        os.makedirs(w, exist_ok=True)
        subprocess.run([sys.executable, FBARES,
                        os.path.join(root, 'n' + n, 'und', 'm_ccx.inp'),
                        os.path.join(w, 'm_ccx.inp'),
                        '--elset', 'Sphere_Only', '--cycles', '1'],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.STDOUT)

        for arm in ARMS:
            jobdir = os.path.join(root, 'n' + n, os.path.dirname(arm))
            if a.dry_run:
                print('  (cd %s && OMP_NUM_THREADS=%d %s m_ccx)'
                      % (jobdir, a.threads, a.ccx))
                continue
            rc = solve(a.ccx, jobdir, a.threads)
            if rc != 0:
                print('  WARNING n=%s %s: ccx exited %d (see solve.log)'
                      % (n, arm, rc))

    if a.dry_run:
        print('\ndry run: %d decks written under %s' % (2 * len(ns) + len(ns),
                                                        root))
        return 0

    rc = subprocess.call(pmake('report_slabconv.py',
                               [root, '%g' % a.kg, a.load] + ns), cwd=HERE)
    return rc


if __name__ == '__main__':
    sys.exit(main())
