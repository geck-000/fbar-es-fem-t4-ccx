#!/usr/bin/env python3
"""Coercivity of the volumetric pairing: Petrov--Galerkin against Galerkin.

The volumetric block of the tangent is

    K_vol(PG) = G^T W H,   G = E D,   H = E A^c D,

the Petrov--Galerkin pairing of the method (test from the unsmoothed edge
strain, trial from the smoothed one).  It is not symmetric, and its
symmetric part can lose positive semi-definiteness.  The total symmetric
part is Sym(K) = K_dev + (K/G) Sym(K_vol), and the method stays coercive
only while the smallest eigenvalue of that sum is positive.  With K_vol
assembled at K = 1, the crossing is the generalized eigenvalue problem

    Sym(K_vol) x = mu K_dev x,      K/G* = -1 / mu_min,

on the Dirichlet-reduced system.  The Galerkin pairing G = H gives the Gram
matrix H^T W H, positive semi-definite for every K/G, so its threshold is
infinite.

The script prints the threshold per mesh and cycle count, the Galerkin
check, the number of negative eigenvalues of the true tangent on both sides
of the threshold, the fluctuation of both pairings on the two-phase sphere
cell, and a compact manufactured-solution table across the transition.

    python3 elements_ccx/tests/verify_pairing.py --ns 6 8 12

`--metrics --ccx <binary>` adds the CalculiX cross-check of the two pairings
on the two-phase cube-inclusion MMS cell, through the real U3 element.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Optional, Tuple

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as proto  # noqa: E402
import verify_convergence as vc  # noqa: E402

G = vc.G_INCLUSION


def blocks(
    n: int, cycles: int, symmetric: bool
) -> Tuple[sp.csr_matrix, sp.csr_matrix]:
    """Deviatoric and volumetric blocks of a single-phase cube.

    Args:
        n: Elements per edge.
        cycles: Cycle count of the volumetric chain.
        symmetric: Use the Galerkin pairing (test = trial = smoothed).

    Returns:
        The pair (K_dev, K_vol) on the Dirichlet-free degrees of freedom.
    """
    nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0, geom=vc.geom_single)
    g, vol = proto.grads(nodes, tets)
    faces, patch = proto.topology(tets, mat)
    fixed = np.repeat(vc.boundary_mask(nodes), 3)
    free = ~fixed
    scheme = 'fbar_%d' % cycles
    if symmetric:
        os.environ['FBAR_SYM'] = '1'
    K_dev = proto.assemble('c3d4', nodes, tets, mat, g, vol, faces, patch,
                           {1: (0.0, 1.0)})[free][:, free].tocsc()
    K_vol = proto.assemble(scheme, nodes, tets, mat, g, vol, faces, patch,
                           {1: (1.0, 0.0)})[free][:, free].tocsc()
    os.environ.pop('FBAR_SYM', None)
    return K_dev, K_vol


def threshold(
    K_dev: sp.csc_matrix, K_vol: sp.csc_matrix, dense: bool = True
) -> Tuple[float, float, int]:
    """Coercivity threshold of a volumetric block.

    Args:
        K_dev: Deviatoric stiffness (SPD on the free space).
        K_vol: Volumetric block at K = 1.
        dense: Solve the generalized problem densely (small meshes only).

    Returns:
        The triple (mu_min, K/G threshold, count of negative eigenvalues of
        Sym(K_vol) with respect to K_dev).
    """
    A = 0.5 * (K_vol + K_vol.T)
    if dense:
        mu = sla.eigvalsh(A.toarray(), K_dev.toarray())
    else:  # pragma: no cover - only used for the largest meshes
        import scipy.sparse.linalg as spl
        mu = spl.eigsh(A.tocsc(), k=A.shape[0] - 1, M=K_dev.tocsc(),
                       which='SA', return_eigenvectors=False)
        mu = np.sort(mu)
    mu_min = float(mu.min())
    neg = int((mu < -1e-10 * max(abs(mu.min()), abs(mu.max()))).sum())
    limit = float('inf') if mu_min >= 0 else float(-1.0 / mu_min)
    return mu_min, limit, neg


def negcount(
    K_dev: sp.csc_matrix, K_vol: sp.csc_matrix, kg: float
) -> Tuple[int, float]:
    """Negative eigenvalues of the true tangent at one K/G.

    Args:
        K_dev: Deviatoric stiffness.
        K_vol: Volumetric block at K = 1.
        kg: The ratio K/G.

    Returns:
        The pair (number of negative eigenvalues, lam_min / lam_max).
    """
    K = (K_dev + kg * K_vol).toarray()
    S = 0.5 * (K + K.T)
    lam = np.linalg.eigvalsh(S)
    scale = max(abs(lam).max(), 1e-300)
    return int((lam < -1e-10 * scale).sum()), float(lam.min() / scale)


def two_phase_blocks(
    n: int, cycles: int, geom: str, jitter: float = 0.3
) -> Tuple[sp.csc_matrix, sp.csc_matrix]:
    """Base stiffness and inclusion volumetric block of a two-phase cell.

    The base holds the whole matrix phase at its own moduli and the
    inclusion deviatoric part; the volumetric block carries the inclusion
    penalty at K = G_i, so the threshold returned by `threshold` is in units
    of the inclusion K/G.

    Args:
        n: Elements per edge.
        cycles: Cycle count.
        geom: 'sphere' or 'bridged' geometry of `smoothing_proto`.
        jitter: Interior jitter.

    Returns:
        The pair (K_base, K_vol_inclusion) on the Dirichlet-free dofs.
    """
    gfun = proto.GEOM[geom]
    nodes, tets, mat = proto.mesh_box(n, 0.375, 0.625, jitter, geom=gfun)
    g, vol = proto.grads(nodes, tets)
    faces, patch = proto.topology(tets, mat)
    fixed = np.repeat(vc.boundary_mask(nodes), 3)
    free = ~fixed
    stiff = (9.37e9 / (3.0 * (1.0 - 2.0 * 0.33)),
             9.37e9 / (2.0 * (1.0 + 0.33)))
    K_base = proto.assemble('fbar_%d' % cycles, nodes, tets, mat, g, vol,
                            faces, patch, {0: stiff, 1: (0.0, G)})[free][:, free]
    K_vol = proto.assemble('fbar_%d' % cycles, nodes, tets, mat, g, vol,
                           faces, patch,
                           {0: (0.0, 0.0), 1: (G, 0.0)})[free][:, free]
    return K_base.tocsc(), K_vol.tocsc()


def pairing_metrics(
    cell: str, kg: float, n: int = 8, ccx: Optional[str] = None
) -> None:
    """Compare the two pairings on a two-phase periodic cell.

    Args:
        cell: 'sphere' or 'bridged'.
        kg: Inclusion K/G.
        n: Elements per edge.
        ccx: CalculiX binary for the cube-inclusion MMS cross-check, or
            None to print the prototype table only.  The periodic cell
            itself is a prototype construction (periodic boundary
            conditions and a homogenisation read-out); the ccx
            cross-check therefore uses `mms_ccx`'s two-phase layout,
            which carries the inclusion arm as the production
            configuration.
    """
    stiff = (9.37e9 / (3.0 * (1.0 - 2.0 * 0.33)),
             9.37e9 / (2.0 * (1.0 + 0.33)))
    props = {0: stiff, 1: (kg * G, G)}
    print('two-phase %s cell at K/G = %g (periodic, n = %d)'
          % (cell, kg, n))
    print('  %-10s %12s %12s %12s' %
          ('pairing', 'C1111', 'fluc', 'p-sch'))
    for sym in (False, True):
        if sym:
            os.environ['FBAR_SYM'] = '1'
        res = proto.run('fbar_1', n, (0.375, 0.625), props, jitter=0.3,
                        geomname=cell)
        os.environ.pop('FBAR_SYM', None)
        print('  %-10s %12.4e %12.4e %12.4e'
              % ('Galerkin' if sym else 'PG', res[0], res[1], res[6]))
    if ccx:
        import mms_ccx
        print('  ccx two-phase MMS cube inclusion at K/G = %g (n = %d)'
              % (kg, n))
        print('  %-10s %-10s %12s %12s %12s'
              % ('pairing', 'arm', 'L2', 'H1', 'incl L2'))
        for sym in (False, True):
            for arm, elsets in (('fbar1_incl', ['INCLUSION']),
                                ('fbar1_all', ['MATRIX', 'INCLUSION'])):
                outdir = os.path.join('out_pairing', 'n%d_%s_kg%g_%s%s'
                                      % (n, cell, kg, arm,
                                         '_sym' if sym else ''))
                res = mms_ccx.run_arm(outdir, n, 'two', kg, arm, ccx, 1,
                                      elsets, symmetric=sym)
                if res is not None:
                    _, l2, h1, l2i, _ = res
                    print('  %-10s %-10s %12.4e %12.4e %12.4e'
                          % ('Galerkin' if sym else 'PG', arm, l2, h1, l2i))


def main() -> int:
    """Run the pairing study.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--ns', type=int, nargs='+', default=[6, 8, 12])
    ap.add_argument('--cell', choices=['cube', 'sphere', 'bridged'],
                    default='cube')
    ap.add_argument('--metrics', action='store_true',
                    help='pairing comparison on a two-phase cell')
    ap.add_argument('--ccx', default=None,
                    help='CalculiX binary for the two-phase MMS cross-check')
    a = ap.parse_args()

    if a.metrics:
        for cell in ('sphere', 'bridged'):
            for kg in (5000.0, 1.0e5):
                pairing_metrics(cell, kg, ccx=a.ccx)
        return 0

    if a.cell != 'cube':
        print('two-phase %s cell: inclusion K/G coercivity threshold of the '
              'Petrov-Galerkin pairing' % a.cell)
        print('  %-4s %-4s %12s %14s %8s' %
              ('n', 'c', 'mu_min', 'K/G*', 'neg mu'))
        for n in a.ns:
            for c in (1, 2):
                K_base, K_vol = two_phase_blocks(n, c, a.cell)
                K_base = 0.5 * (K_base + K_base.T)
                mu_min, limit, neg = threshold(K_base, K_vol, dense=n <= 10)
                print('  %-4d %-4d %12.4e %14s %8d'
                      % (n, c, mu_min,
                         ('inf' if limit == float('inf') else '%.4g' % limit),
                         neg))
        return 0

    print('coercivity threshold of the Petrov-Galerkin pairing')
    print('  %-4s %-4s %12s %14s %8s' %
          ('n', 'c', 'mu_min', 'K/G*', 'neg mu'))
    for n in a.ns:
        for c in (1, 2):
            K_dev, K_vol = blocks(n, c, symmetric=False)
            mu_min, limit, neg = threshold(K_dev, K_vol, dense=n <= 12)
            print('  %-4d %-4d %12.4e %14s %8d'
                  % (n, c, mu_min,
                     ('inf' if limit == float('inf') else '%.4g' % limit),
                     neg))
        K_dev, K_vol = blocks(n, 1, symmetric=True)
        mu_min, limit, neg = threshold(K_dev, K_vol, dense=n <= 12)
        print('  %-4d %-4s %12.4e %14s %8d   (Galerkin)'
              % (n, 'sym', mu_min,
                 ('inf' if limit == float('inf') else '%.4g' % limit), neg))

    print('\nnegative eigenvalues of Sym(K) around the threshold (n = 8, c = 1)')
    K_dev, K_vol = blocks(8, 1, symmetric=False)
    for kg in (1e2, 1e3, 5e3, 1e4, 1e5, 1e6):
        neg, ratio = negcount(K_dev, K_vol, kg)
        print('  K/G=%-8g  n_neg = %-4d  lam_min/lam_max = %11.3e'
              % (kg, neg, ratio))

    print('\ntwo-phase sphere cell, fluctuation at K/G = 5000')
    stiff = (9.37e9 / (3 * (1 - 2 * 0.33)), 9.37e9 / (2 * (1 + 0.33)))
    props = {0: stiff, 1: (5000.0 * G, G)}
    for sym in (False, True):
        if sym:
            os.environ['FBAR_SYM'] = '1'
        res = proto.run('fbar_1', 8, (0.375, 0.625), props, jitter=0.3,
                        geomname='sphere')
        os.environ.pop('FBAR_SYM', None)
        print('  %-8s  fluc = %10.3e   c1111 = %10.3e'
              % ('Galerkin' if sym else 'PG', res[1], res[0]))

    print('\nmanufactured solution across the transition (single phase, n = 10)')
    print('  %-10s %-8s %12s' % ('pairing', 'K/G', 'H1 error'))
    for sym in (False, True):
        if sym:
            os.environ['FBAR_SYM'] = '1'
        for kg in (1e3, 1e4, 1e5, 1e6):
            rows = vc.study('fbar_1', vc.geom_single, {1: (kg * G, G)},
                            'fbar_1', [10], None)
            print('  %-10s %-8g %12.4e'
                  % ('Galerkin' if sym else 'PG', kg, rows[-1][3]))
        os.environ.pop('FBAR_SYM', None)
    return 0


if __name__ == '__main__':
    sys.exit(main())
