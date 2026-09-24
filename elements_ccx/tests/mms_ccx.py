#!/usr/bin/env python3
"""Manufactured-solution convergence through the real CalculiX element.

Generates the divergence-free manufactured problem of
`verify_convergence.py` as a CalculiX deck, optionally converts the F-bar
element set with `fbares.py`, runs `ccx_fbar`, parses the nodal
displacements from the .dat and measures the same relative L2 and H1
errors as the prototype.  The two levels therefore quote the same norms.

Layouts and arms:

    single   one material; arms c3d4, fbar0, fbar1
    two      cube inclusion [0.25, 0.75]^3; arms c3d4,
             fbar1_incl (inclusion only, the production configuration)
             and fbar1_all

    python3 elements_ccx/tests/mms_ccx.py --out out_mms --n 8 \
        --layout single --kg 5000 --arms c3d4 fbar1

    elements_ccx/tests/mms_ccx.sh        # the full sweep
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as proto  # noqa: E402
import verify_convergence as vc  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', '..'))
FBARES = os.path.join(REPO, 'elements_ccx', 'fbares.py')


def e_nu(K: float, G: float) -> Tuple[float, float]:
    """Convert bulk and shear modulus to the pair CalculiX wants.

    Args:
        K: Bulk modulus.
        G: Shear modulus.

    Returns:
        The pair (Young's modulus, Poisson's ratio).
    """
    e = 9.0 * K * G / (3.0 * K + G)
    nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return e, nu


def nset(ids: Sequence[int], name: str, per_line: int = 12) -> List[str]:
    """Emit an *NSET block.

    Args:
        ids: One-based node ids.
        name: Set name.
        per_line: Entries per data line.

    Returns:
        The deck lines.
    """
    out = ['*NSET,NSET=' + name]
    s = [str(int(i)) for i in ids]
    out += [', '.join(s[c:c + per_line]) for c in range(0, len(s), per_line)]
    return out


def deck(
    nodes: np.ndarray, tets: np.ndarray, mat: np.ndarray, vol: np.ndarray,
    props: Dict[int, Tuple[float, float]], layout: str,
) -> List[str]:
    """Write the model and step of the manufactured problem.

    Args:
        nodes: Node coordinates.
        tets: Element connectivity (0-based).
        mat: Material index per element.
        vol: Element volumes.
        props: Material moduli keyed by material index.
        layout: 'single' or 'two'.

    Returns:
        The deck as a list of lines.
    """
    lines = ['*HEADING',
             ' MMS divergence-free convergence problem',
             '*NODE']
    for i, p in enumerate(nodes):
        lines.append('%d, %.12e, %.12e, %.12e' % (i + 1, p[0], p[1], p[2]))

    if layout == 'single':
        lines.append('*ELEMENT,TYPE=C3D4,ELSET=BODY')
        for e, t in enumerate(tets):
            lines.append('%d, %s' % (e + 1, ', '.join(str(int(x) + 1) for x in t)))
        lines += nset(range(1, len(nodes) + 1), 'ALLNODES')
        lines += ['*SOLID SECTION,ELSET=BODY,MATERIAL=Mat_Body',
                  '*MATERIAL,NAME=Mat_Body', '*ELASTIC']
        e, nu = e_nu(*props[1])
        lines.append('%.12e, %.12e' % (e, nu))
    else:
        for m, name in ((0, 'MATRIX'), (1, 'INCLUSION')):
            ids = [e for e in range(len(tets)) if mat[e] == m]
            lines.append('*ELEMENT,TYPE=C3D4,ELSET=' + name)
            for e in ids:
                lines.append('%d, %s' % (e + 1, ', '.join(str(int(x) + 1)
                                                          for x in tets[e])))
        lines += nset(range(1, len(nodes) + 1), 'ALLNODES')
        for m, name in ((0, 'MATRIX'), (1, 'INCLUSION')):
            lines += ['*SOLID SECTION,ELSET=%s,MATERIAL=Mat_%s' % (name, name),
                      '*MATERIAL,NAME=Mat_%s' % name, '*ELASTIC']
            e, nu = e_nu(*props[m])
            lines.append('%.12e, %.12e' % (e, nu))

    fixed = vc.boundary_mask(nodes)
    ux = vc.u_star(nodes)
    f = vc.body_loads(nodes, tets, mat, vol, props)
    lines.append('*STEP')
    lines.append('*STATIC,SOLVER=PARDISO')
    lines.append('*BOUNDARY')
    for i in np.flatnonzero(fixed):
        for c in range(3):
            lines.append('%d, %d, %d, %.12e' % (i + 1, c + 1, c + 1, ux[i, c]))
    lines.append('*CLOAD')
    for i in np.flatnonzero(~fixed):
        for c in range(3):
            lines.append('%d, %d, %.12e' % (i + 1, c + 1, f[3 * i + c]))
    lines += ['*NODE PRINT,NSET=ALLNODES', 'U', '*END STEP']
    return lines


def read_dat_u(path: str) -> Optional[np.ndarray]:
    """Read the last displacements block of a CalculiX .dat file.

    Args:
        path: Path of the .dat file.

    Returns:
        The values ordered by node id as an (n_node, 3) array, or None when
        no displacements block was written.
    """
    values: Dict[int, List[float]] = {}
    active = False
    started = False
    with open(path, 'r', errors='replace') as handle:
        for line in handle:
            s = line.strip()
            low = s.lower()
            if low.startswith('displacements'):
                active = True
                started = False
                values = {}
                continue
            if not active:
                continue
            if not s:
                if started:
                    active = False
                continue
            parts = s.split()
            try:
                nid = int(parts[0])
                vals = [float(parts[1]), float(parts[2]), float(parts[3])]
            except (ValueError, IndexError):
                if started:
                    active = False
                continue
            values[nid] = vals
            started = True
    if not values:
        return None
    n = max(values)
    out = np.zeros((n, 3))
    for nid, vals in values.items():
        out[nid - 1] = vals
    return out


def run_arm(
    outdir: str, n: int, layout: str, kg: float, arm: str, ccx: str,
    cycles: Optional[int], elsets: Optional[List[str]],
    mesh: str = 'box', symmetric: bool = False,
) -> Optional[Tuple[float, float, float, float, float]]:
    """Generate, solve and post-process one arm.

    Args:
        outdir: Directory that receives the deck and results.
        n: Elements per edge.
        layout: 'single' or 'two'.
        kg: Inclusion K/G.
        arm: Arm label.
        ccx: CalculiX binary.
        cycles: Cycles for the conversion, None for plain C3D4.
        elsets: Element sets to convert.
        mesh: 'box' or 'delaunay' (single phase only).
        symmetric: Generate for the Galerkin pairing and run with
            CCX_FBAR_SYM=1.

    Returns:
        The tuple (h, rel L2, rel H1, rel L2 inclusion), or None on failure.
    """
    os.makedirs(outdir, exist_ok=True)
    if layout == 'single':
        geom = vc.geom_single
        props = {1: vc.iso(kg)}
    else:
        geom = vc.geom_cube
        props = {0: (2.1666667 * vc.G_INCLUSION, vc.G_INCLUSION),
                 1: vc.iso(kg)}
    if mesh == 'delaunay':
        nodes, tets, mat = vc.build_mesh('delaunay', n)
    else:
        nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0, geom=geom)
    g, vol = proto.grads(nodes, tets)
    src = os.path.join(outdir, 'm.inp')
    with open(src, 'w') as handle:
        handle.write('\n'.join(deck(nodes, tets, mat, vol, props, layout)) + '\n')
    deckname = 'm'
    if cycles is not None:
        deckname = 'mf'
        dst = os.path.join(outdir, 'mf.inp')
        cmd = [sys.executable, FBARES, src, dst, '--cycles', str(cycles)]
        for es in elsets or []:
            cmd += ['--elset', es]
        if symmetric:
            cmd += ['--symmetric']
        subprocess.run(cmd, check=True)
    cmd = [ccx, deckname]
    env = dict(os.environ)
    if cycles is not None:
        env['CCX_FBAR_C'] = str(cycles)
    with open(os.path.join(outdir, 'solve.log'), 'w') as log:
        subprocess.run(cmd, cwd=outdir, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    dat = os.path.join(outdir, deckname + '.dat')
    if not os.path.exists(dat):
        print('   %s: no .dat written' % arm)
        return None
    u = read_dat_u(dat)
    if u is None or len(u) != len(nodes):
        print('   %s: displacements not parsed (%s)'
              % (arm, None if u is None else len(u)))
        return None
    l2, h1 = vc.errors(nodes, tets, mat, g, vol, u.reshape(-1))
    l2i, h1i = vc.errors(nodes, tets, mat, g, vol, u.reshape(-1),
                         phase=None if layout == 'single' else 1)
    return (1.0 / n, l2, h1, l2i, h1i)


def main() -> int:
    """Command-line entry point.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='out_mms')
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--layout', choices=['single', 'two'], required=True)
    ap.add_argument('--kg', type=float, required=True)
    ap.add_argument('--arms', nargs='+', required=True)
    ap.add_argument('--mesh', choices=['box', 'delaunay'], default='box')
    ap.add_argument('--symmetric', action='store_true',
                    help='Galerkin pairing via CCX_FBAR_SYM=1')
    ap.add_argument('--ccx', default=os.environ.get('CCX_MMS', 'ccx_fbar'))
    a = ap.parse_args()
    if a.mesh == 'delaunay' and a.layout == 'two':
        raise SystemExit('mms_ccx: the two-phase cube geometry needs the '
                         'structured mesh')
    tag = (a.mesh if a.mesh != 'box' else '') + ('_sym' if a.symmetric else '')

    arms = {
        'c3d4': (None, None),
        'fbar0': (0, ['BODY'] if a.layout == 'single' else ['INCLUSION']),
        'fbar1': (1, ['BODY'] if a.layout == 'single' else ['INCLUSION']),
        'fbar1_incl': (1, ['INCLUSION']),
        'fbar1_all': (1, ['MATRIX', 'INCLUSION']),
    }
    for arm in a.arms:
        cycles, elsets = arms[arm]
        d = os.path.join(a.out, 'n%d_%s_kg%g_%s%s'
                         % (a.n, a.layout, a.kg, arm, tag))
        res = run_arm(d, a.n, a.layout, a.kg, arm, a.ccx, cycles, elsets,
                      mesh=a.mesh, symmetric=a.symmetric)
        if res is not None:
            h, l2, h1, l2i, h1i = res
            print('  n=%-3d %-10s %-10s h %.4f  L2 %.4e  H1 %.4e  incl L2 %.4e'
                  % (a.n, arm, tag or a.mesh, h, l2, h1, l2i))
            csv_path = os.path.join(a.out, 'mms.csv')
            with open(csv_path, 'a') as handle:
                handle.write('%d,%s,%g,%s,%s,%.8e,%.8e,%.8e,%.8e,%.8e\n'
                             % (a.n, a.layout, a.kg, arm, tag or a.mesh, h,
                                l2, h1, l2i, h1i))
    return 0


if __name__ == '__main__':
    sys.exit(main())
