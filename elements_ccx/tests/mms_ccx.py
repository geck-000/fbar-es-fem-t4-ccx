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

    python3 elements_ccx/tests/mms_ccx.py --sweep

--sweep runs n = 8 12 16 24 over K/G = 10 100 1000 5000 (single phase, plus
fbar0 at 5000) and the same K/G list for the two-phase layout, then prints
the rates and the K/G summary.
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


CSV_HEADER = 'n,layout,kg,arm,h,l2,h1,l2_incl,h1_incl\n'


def arm_spec(layout: str, arm: str
             ) -> Tuple[Optional[int], Optional[List[str]]]:
    """Return (cycles, elsets) for one arm name.

    Args:
        layout: 'single' or 'two'.
        arm: Arm label.

    Returns:
        The cycle count and the element sets to convert.
    """
    return {
        'c3d4': (None, None),
        'fbar0': (0, ['BODY'] if layout == 'single' else ['INCLUSION']),
        'fbar1': (1, ['BODY'] if layout == 'single' else ['INCLUSION']),
        'fbar1_incl': (1, ['INCLUSION']),
        'fbar1_all': (1, ['MATRIX', 'INCLUSION']),
    }[arm]


def run_arms(out: str, n: int, layout: str, kg: float, arms: List[str],
             ccx: str, mesh: str = 'box', symmetric: bool = False) -> None:
    """Run a list of arms on one mesh and append them to ``out/mms.csv``.

    Args:
        out: Output directory.
        n: Elements per edge.
        layout: 'single' or 'two'.
        kg: Inclusion K/G.
        arms: Arm labels.
        ccx: CalculiX binary.
        mesh: 'box' or 'delaunay'.
        symmetric: Generate and run the Galerkin pairing.
    """
    tag = (mesh if mesh != 'box' else '') + ('_sym' if symmetric else '')
    os.makedirs(out, exist_ok=True)
    csv_path = os.path.join(out, 'mms.csv')
    if not os.path.exists(csv_path):
        with open(csv_path, 'w') as handle:
            handle.write(CSV_HEADER)
    for arm in arms:
        cycles, elsets = arm_spec(layout, arm)
        d = os.path.join(out, 'n%d_%s_kg%g_%s%s' % (n, layout, kg, arm, tag))
        res = run_arm(d, n, layout, kg, arm, ccx, cycles, elsets,
                      mesh=mesh, symmetric=symmetric)
        if res is None:
            continue
        h, l2, h1, l2i, h1i = res
        print('  n=%-3d %-10s %-10s h %.4f  L2 %.4e  H1 %.4e  incl L2 %.4e'
              % (n, arm, tag or mesh, h, l2, h1, l2i, h1i))
        with open(csv_path, 'a') as handle:
            handle.write('%d,%s,%g,%s,%s,%.8e,%.8e,%.8e,%.8e,%.8e\n'
                         % (n, layout, kg, arm, tag or mesh, h, l2, h1, l2i,
                            h1i))


def load_rows(path: str) -> List[Tuple[int, str, float, str, float, float,
                                       float, float, float]]:
    """Read the sweep CSV.

    Args:
        path: Path of mms.csv.

    Returns:
        The rows as tuples.
    """
    import csv
    rows = []
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
        log(a/b) / log(nb/na), or nan when b is zero.
    """
    import math
    return float('nan') if b <= 0 else math.log(a / b) / math.log(nb / na)


def report(path: str) -> None:
    """Print errors per mesh, last-pair rates and the K/G summary.

    Args:
        path: Path of mms.csv.
    """
    rows = load_rows(path)
    groups: Dict[Tuple[str, float, str], list] = {}
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
        print('\nK/G summary at n=%d (%s phase)' % (finest[layout], layout))
        print('  %-10s %-14s %10s %10s %10s'
              % ('arm', 'K/G', 'L2', 'H1', 'incl H1'))
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


def main() -> int:
    """Command-line entry point.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='out_mms')
    ap.add_argument('--sweep', action='store_true',
                    help='run the default sweep and print the report')
    ap.add_argument('--ns', type=int, nargs='+', default=None)
    ap.add_argument('--kgs', type=float, nargs='+', default=None)
    ap.add_argument('--two-kgs', type=float, nargs='+', default=None)
    ap.add_argument('--n', type=int, default=None)
    ap.add_argument('--layout', choices=['single', 'two'], default=None)
    ap.add_argument('--kg', type=float, default=None)
    ap.add_argument('--arms', nargs='+', default=None)
    ap.add_argument('--mesh', choices=['box', 'delaunay'], default='box')
    ap.add_argument('--symmetric', action='store_true',
                    help='Galerkin pairing via CCX_FBAR_SYM=1')
    ap.add_argument('--ccx', default=os.environ.get('CCX_MMS', 'ccx_fbar'))
    a = ap.parse_args()

    if a.sweep:
        ns = a.ns or [8, 12, 16, 24]
        kgs = a.kgs or [10, 100, 1000, 5000]
        two_kgs = a.two_kgs or [10, 100, 1000, 5000]
        for n in ns:
            for kg in kgs:
                run_arms(a.out, n, 'single', kg, ['c3d4', 'fbar1'], a.ccx,
                         mesh=a.mesh, symmetric=a.symmetric)
            run_arms(a.out, n, 'single', 5000, ['fbar0'], a.ccx,
                     mesh=a.mesh, symmetric=a.symmetric)
            for kg in two_kgs:
                run_arms(a.out, n, 'two', kg,
                         ['c3d4', 'fbar1_incl', 'fbar1_all'], a.ccx,
                         mesh=a.mesh, symmetric=a.symmetric)
        report(os.path.join(a.out, 'mms.csv'))
        return 0

    if a.mesh == 'delaunay' and a.layout == 'two':
        raise SystemExit('mms_ccx: the two-phase cube geometry needs the '
                         'structured mesh')
    missing = [k for k in ('n', 'layout', 'kg', 'arms') if getattr(a, k) is None]
    if missing:
        ap.error('missing arguments: ' + ', '.join('--' + m for m in missing))
    run_arms(a.out, a.n, a.layout, a.kg, a.arms, a.ccx, mesh=a.mesh,
             symmetric=a.symmetric)
    return 0


if __name__ == '__main__':
    sys.exit(main())
