#!/usr/bin/env python3
"""Run the F-barES-FEM-T4 test suite.

    python3 run_tests.py                  # element checks
    python3 run_tests.py --quick          # operators, walk and deck only
    python3 run_tests.py --ccx /path/to/ccx   # adds the CalculiX checks

`verify_element.py` carries the checks that need no CalculiX binary.
`verify_nltan.py` compiles the Fortran tangent driver with gfortran and runs
always.  With `--ccx`, the suite adds `verify_nlccx.py` (the NLGEOM force
through CalculiX) and `verify_pairing.py --metrics --ccx` (the two pairings
cross-checked against the compiled element).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from typing import List, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))


def run(name: str, args: List[str]) -> int:
    """Run one suite member and print its banner.

    Args:
        name: Display name of the member.
        args: Command and arguments.

    Returns:
        The member's exit status.
    """
    print('\n==================== %s ====================' % name)
    rc = subprocess.call([sys.executable] + args, cwd=HERE)
    print('-------------------- %s: %s' % (name, 'OK' if rc == 0
                                          else 'FAILED (exit %d)' % rc))
    return rc


def main() -> int:
    """Run the selected suite members.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ccx', default=None,
                    help='path to the patched CalculiX binary')
    ap.add_argument('--quick', action='store_true',
                    help='pass --quick to verify_element.py')
    ap.add_argument('--skip-stability', action='store_true')
    ap.add_argument('--skip-tangent', action='store_true')
    ap.add_argument('--skip-fortran', action='store_true',
                    help='skip verify_nltan.py (needs gfortran)')
    a = ap.parse_args()

    element = ['verify_element.py']
    if a.quick:
        element.append('--quick')
    if a.skip_stability:
        element.append('--skip-stability')
    if a.skip_tangent:
        element.append('--skip-tangent')

    results: List[Tuple[str, int]] = []
    results.append(('element checks', run('element checks', element)))
    if not a.skip_fortran:
        results.append(('fortran tangent', run('fortran tangent',
                                               ['verify_nltan.py'])))
    if a.ccx:
        results.append(('NLGEOM through CalculiX',
                        run('NLGEOM through CalculiX',
                            ['verify_nlccx.py', '--ccx', a.ccx])))
        results.append(('pairing cross-check',
                        run('pairing cross-check',
                            ['verify_pairing.py', '--metrics',
                             '--ccx', a.ccx])))

    print('\n==================== summary ====================')
    failed = 0
    for name, rc in results:
        print('  %-28s %s' % (name, 'OK' if rc == 0 else 'FAILED'))
        failed += rc != 0
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
