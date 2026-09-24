#!/usr/bin/env python3
"""Verify the F-barES-FEM-T4 element against closed-form answers.

One entry point for the checks that need no CalculiX binary:

  V  smoothing operators, patch test, c = 0 identity, constraint rank
  W  the u3vol.f element-weight walk against the operator S = E A^c
  N  finite-strain internal force (N1--N7)
  T  consistent tangent and quadratic Newton (T1--T3)
  S  spurious-mode census on the two-phase clamped cube
  D  deck generator walk, set equality both ways, deck format rules
  A  assembled-structure audit of a generated deck (informational)

usage:
    python3 verify_element.py [--quick] [--skip-stability] [--skip-tangent]
                              [--n N] [--deck ORIGINAL GENERATED [CYCLES]]
                              [--sample N]

Without --deck a small cell is generated with make_slabconv.py and fbares.py
and the deck checks run on it.
"""
from __future__ import annotations

import argparse
import os
import random
import re
import subprocess
import sys
import tempfile
import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import scipy.sparse.linalg as spl

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import smoothing_proto as P  # noqa: E402
import smoothing_proto as S  # noqa: E402

EI, NI = 9.37e9, 0.33
ICE = (EI / (3.0 * (1.0 - 2.0 * NI)), EI / (2.0 * (1.0 + NI)))
GB = 4.4e5
FAILS: List[str] = []


def chk(name: str, got: float, want: float = 0.0, tol: float = 0.0) -> None:
    """Print and record one numeric check.

    Args:
        name: Label printed in the table.
        got: Measured value.
        want: Expected value.
        tol: Absolute tolerance on ``got - want``.
    """
    ok = abs(got - want) <= tol
    print('   %-58s %12.4e  %s' % (name, got, 'PASS' if ok else 'FAIL'))
    if not ok:
        FAILS.append(name)


def chk_ok(name: str, ok: bool, detail: str = '') -> None:
    """Print and record one boolean check.

    Args:
        name: Label printed in the table.
        ok: Whether the check passed.
        detail: Extra text printed on failure.
    """
    print('   %-58s %s%s' % (name, 'PASS' if ok else 'FAIL',
                             '' if ok else '  ' + detail))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------------------
# V  smoothing operators, patch test, c = 0 identity, constraint rank
# ---------------------------------------------------------------------------

def check_operators() -> None:
    """Run the operator and patch-test checks (V1--V4)."""
    print('V1  smoothing operators: unit row sums, constant field preserved')
    nodes, tets, mat = P.mesh_box(6, 0.375, 0.625, 0.3, geom=P.GEOM['sphere'])
    g, vol = P.grads(nodes, tets)
    for jin in ('elem', 'edge'):
        for c in (0, 1, 2, 3):
            sel, idx, ek, ed, Vh, E, S, R = P.fbar_es_operators(
                nodes, tets, mat, g, vol, c, 1, jinput=jin)
            for nm, Op in (('E', E), ('R', R), ('S', S)):
                rs = np.asarray(Op.sum(axis=1)).ravel()
                chk('%s c=%d  max|rowsum(%s) - 1|' % (jin, c, nm),
                    np.abs(rs - 1).max(), 0.0, 1e-12)
            one = np.ones(len(sel))
            chk('%s c=%d  max|S.1 - 1| (constant J preserved)' % (jin, c),
                np.abs(S @ one - 1).max(), 0.0, 1e-11)

    print('\nV2  patch test, homogeneous block, exact C1111 = K + 4G/3')
    exact = ICE[0] + 4.0 * ICE[1] / 3.0
    for jin in ('elem', 'edge'):
        os.environ['FBAR_JIN'] = jin
        for c in (0, 1, 2, 3):
            r = P.run('fbar_%d' % c, 4, (2.0, 3.0), {0: ICE, 1: ICE},
                      jitter=0.0)
            chk('%s fbar_%d  rel err in C1111' % (jin, c),
                abs(r[0] / exact - 1), 0.0, 1e-10)
    os.environ['FBAR_JIN'] = 'elem'

    print('\nV3  c = 0 must reproduce selective ES-FEM-T4 exactly (S = E)')
    sel, idx, ek, ed, Vh, E, S, R = P.fbar_es_operators(
        nodes, tets, mat, g, vol, 0, 1, jinput='elem')
    chk('max|S - E| at c=0', abs(S - E).max() if (S - E).nnz else 0.0,
        0.0, 1e-14)

    print('\nV4  rank of the volumetric operator')
    for c in (0, 1, 2):
        sel, idx, ek, ed, Vh, E, S, R = P.fbar_es_operators(
            nodes, tets, mat, g, vol, c, 1, jinput='elem')
        nnode = len({int(a) for e in sel for a in tets[e]})
        rE = np.linalg.matrix_rank(E.toarray())
        rS = np.linalg.matrix_rank(S.toarray())
        print('   c=%d  brine: %d elems, %d nodes, %d edges | rank E = %d, '
              'rank S = %d  (c >= 1 caps the trial space at the incidence rank)'
              % (c, len(sel), nnode, len(ek), rE, rS))


# ---------------------------------------------------------------------------
# W  the u3vol.f walk against S = E A^c
# ---------------------------------------------------------------------------

def walk(
    nodes: np.ndarray,
    tets: np.ndarray,
    mat: np.ndarray,
    g: object,
    vol: np.ndarray,
    imat: int,
    edge: Tuple[int, int],
    ncyc: int,
) -> Tuple[Optional[float], Optional[Dict[int, float]], Optional[Dict[int, float]]]:
    """Mirror u3vol.f: eq. (8) ring, then ncyc cycles of A = P Q."""
    at: Dict[int, List[int]] = {}
    for e in range(len(tets)):
        if mat[e] != imat:
            continue
        for a in tets[e]:
            at.setdefault(int(a), []).append(e)
    na, nb = edge
    cur: Dict[int, float] = {}
    for e in at.get(na, []):
        t = [int(x) for x in tets[e]]
        if nb in t:
            cur[e] = vol[e] / 6.0
    vh = sum(cur.values())
    if vh <= 0:
        return None, None, None
    cur = {e: w / vh for e, w in cur.items()}
    ring = dict(cur)                                  # tbar's element weights
    for _ in range(ncyc):
        nxt: Dict[int, float] = {}
        for je, wgt in cur.items():
            for k in (int(x) for x in tets[je]):
                # V_n of eq. (6), THIS MATERIAL ONLY
                vn = sum(vol[e] / 4.0 for e in at[k])
                if vn <= 0:
                    continue
                for nl in at[k]:
                    nxt[nl] = nxt.get(nl, 0.0) \
                        + wgt * 0.25 * vol[nl] / 4.0 / vn
        cur = nxt
    return vh, ring, cur


def check_walk() -> None:
    """Run the u3vol.f row-walk checks against the operator chain."""
    print('\nW   u3vol.f element-weight walk against S = E A^c')
    for ncyc in (1, 2):
        for jitter in (0.0, 0.3):
            nodes, tets, mat = P.mesh_box(6, 0.375, 0.625, jitter,
                                          geom=P.GEOM['sphere'])
            g, vol = P.grads(nodes, tets)
            m = 1
            sel, idx, ekeys, edges, Vh, E, S, R = P.fbar_es_operators(
                nodes, tets, mat, g, vol, ncyc, m)
            eref = {e: i for i, e in enumerate(sel)}
            worst_v, worst_w, worst_rs = 0.0, 0.0, 0.0
            rng = np.random.default_rng(3)
            pick = rng.choice(len(ekeys), size=min(150, len(ekeys)),
                              replace=False)
            for h in pick:
                k = ekeys[h]
                vh, ring, cur = walk(nodes, tets, mat, g, vol, m, k, ncyc)
                worst_v = max(worst_v, abs(vh / Vh[h] - 1.0))
                row = S.getrow(h).toarray().ravel()
                got = np.zeros_like(row)
                for e, w in cur.items():
                    got[eref[e]] = w
                den = np.abs(row).max()
                worst_w = max(worst_w, np.abs(got - row).max() / den)
                worst_rs = max(worst_rs, abs(sum(cur.values()) - 1.0))
            for nm, val, tol in (('V_h vs prototype', worst_v, 1e-13),
                                 ('chain row vs S = E A^c', worst_w, 1e-12),
                                 ('unit row sum (constant J)', worst_rs, 1e-12)):
                ok = val <= tol
                print('   c=%d jit=%.1f  %-30s %11.3e  %s'
                      % (ncyc, jitter, nm, val, 'PASS' if ok else 'FAIL'))
                if not ok:
                    FAILS.append('walk c=%d jit=%.1f %s' % (ncyc, jitter, nm))


# ---------------------------------------------------------------------------
# N  finite-strain internal force
# ---------------------------------------------------------------------------

def check_force() -> None:
    """Run the finite-strain force checks (N1--N7)."""
    print('\nN1  f(0) = 0')
    nodes, tets, mat = P.mesh_box(6, 0.375, 0.625, 0.3, geom=P.GEOM['sphere'])
    g, vol = P.grads(nodes, tets)
    props = {0: ICE, 1: (500.0 * GB, GB)}
    for c in (0, 1, 2, 3):
        nl = P.FbarNL(nodes, tets, mat, g, vol, props, c)
        f = nl.force(np.zeros(3 * len(nodes)))
        chk('c=%d  max|f(0)|' % c, np.abs(f).max(), 0.0, 1e-6)

    print('\nN2  rigid translation gives f = 0')
    u = np.zeros(3 * len(nodes))
    u[0::3] = 0.37
    u[1::3] = -0.21
    u[2::3] = 0.08
    for c in (0, 2):
        nl = P.FbarNL(nodes, tets, mat, g, vol, props, c)
        chk('c=%d  max|f(rigid)|' % c, np.abs(nl.force(u)).max(), 0.0, 1e-6)

    print('\nN3  homogeneous uniform stretch: affine field is exact equilibrium')
    # ONE material.  The smoothing is not allowed to cross a material
    # interface, so a patch test needs a single-material mesh; the periodic
    # statement (N3b) uses a badly distorted cell instead.
    hom = {0: ICE}
    print('  N3a  unjittered box, free surfaces, analytic Hencky resultant')
    n0, t0, m0 = P.mesh_box(6, 0.375, 0.625, 0.0, geom=P.GEOM['sphere'])
    g0, v0 = P.grads(n0, t0)
    m0 = np.zeros_like(m0)
    int0 = np.ones(len(n0), dtype=bool)
    for d in range(3):
        int0 &= (n0[:, d] > 1e-9) & (n0[:, d] < 1.0 - 1e-9)
    face0 = np.abs(n0[:, 0] - 1.0) < 1e-9
    for c in (0, 1, 2, 3):
        nlh = P.FbarNL(n0, t0, m0, g0, v0, hom, c)
        for lam in (1.0 + 1e-3, 1.2):
            uu = np.zeros(3 * len(n0))
            uu[0::3] = (lam - 1.0) * n0[:, 0]
            f = nlh.force(uu)
            scale = np.abs(f).max() or 1.0
            chk('c=%d lam=%.3f  max|f_interior| / max|f|' % (c, lam),
                np.abs(f.reshape(-1, 3)[int0]).max() / scale, 0.0, 1e-9)
            # F = diag(lam,1,1), so T_xx = (K + 4G/3) ln(lam) and the x = 1
            # face keeps unit current area: its resultant is T_xx itself.
            Txx = (ICE[0] + 4.0 * ICE[1] / 3.0) * np.log(lam)
            fx = f.reshape(-1, 3)[face0, 0].sum()
            chk('c=%d lam=%.3f  face resultant vs analytic Hencky' % (c, lam),
                abs(fx / Txx - 1.0), 0.0, 1e-9)

    print('  N3b  jittered periodic cell: affine field is an exact solution')
    for c in (0, 1, 2, 3):
        cn, fln, it, rn = P.run_nl('fbar_%d' % c, 6, (0.375, 0.625),
                                   {0: ICE, 1: ICE}, eps=1e-3, jitter=0.3,
                                   geomname='sphere')
        chk('c=%d  fluctuation of the homogeneous cell' % c, fln, 0.0, 1e-7)
        exact = ICE[0] + 4.0 * ICE[1] / 3.0
        chk('c=%d  C1111 vs K + 4G/3 (Hencky, eps=1e-3)' % c,
            abs(cn / exact - 1.0), 0.0, 2e-3)

    print('\nN4  finite strain -> small strain as eps -> 0 (two-phase, K/G = 500)')
    for c in (0, 1, 2):
        lin = P.run('fbar_%d' % c, 6, (0.375, 0.625), props, eps=1e-3,
                    jitter=0.3, geomname='sphere')[0]
        for eps in (1e-3, 1e-5):
            cn, fln, it, rn = P.run_nl('fbar_%d' % c, 6, (0.375, 0.625), props,
                                       eps=eps, jitter=0.3, geomname='sphere')
            rel = abs(cn / lin - 1.0)
            print('   c=%d eps=%.0e  C1111_nl=%.6e  vs lin %.6e  rel %8.2e  '
                  '(%d its, |r|=%.1e)' % (c, eps, cn, lin, rel, it, rn))
            if eps == 1e-5:
                chk('c=%d  |C1111_nl/C1111_lin - 1| at eps=1e-5' % c,
                    rel, 0.0, 5e-4)

    print('\nN5  frame indifference: a rigid rotation carries no force')
    th = 0.7
    ax = np.array([1.0, 2.0, 3.0])
    ax /= np.linalg.norm(ax)
    Kx = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    Q = np.eye(3) + np.sin(th) * Kx + (1 - np.cos(th)) * (Kx @ Kx)
    urot = ((Q - np.eye(3)) @ nodes.T).T.ravel()
    for c in (0, 1, 2, 3):
        nl = P.FbarNL(nodes, tets, mat, g, vol, props, c)
        fr = nl.force(urot)
        # scale by the force a comparable STRETCH produces
        ustr = np.zeros(3 * len(nodes))
        ustr[0::3] = 0.1 * nodes[:, 0]
        ref = np.abs(nl.force(ustr)).max()
        chk('c=%d  max|f(rigid rotation)| / max|f(stretch)|' % c,
            np.abs(fr).max() / ref, 0.0, 1e-11)

    print('\nN6  the linear operator IS df/du at the reference configuration')
    n4, t4, m4 = P.mesh_box(4, 0.375, 0.625, 0.3, geom=P.GEOM['sphere'])
    g4, v4 = P.grads(n4, t4)
    fc4, pt4 = P.topology(t4, m4)
    nd = 3 * len(n4)
    for c in (0, 1, 2):
        nlj = P.FbarNL(n4, t4, m4, g4, v4, props, c)
        Klin = P.assemble('fbar_%d' % c, n4, t4, m4, g4, v4, fc4, pt4, props,
                          0.0, False, None).toarray()
        hstep = 1e-7 / 4.0
        Jnum = np.zeros((nd, nd))
        for j in range(nd):
            e = np.zeros(nd)
            e[j] = hstep
            Jnum[:, j] = (nlj.force(e) - nlj.force(-e)) / (2.0 * hstep)
        den = np.abs(Klin).max()
        chk('c=%d  max|df/du - K_lin| / max|K_lin|' % c,
            np.abs(Jnum - Klin).max() / den, 0.0, 2e-6)
        asym = np.abs(Klin - Klin.T).max() / den
        print('       c=%d  asymmetry max|K - K^T|/max|K| = %.3e %s'
              % (c, asym, '(non-symmetric, as eq. 17 requires)' if c else
                 '(symmetric, as S = E at c=0 requires)'))
        if c == 0 and asym > 1e-12:
            FAILS.append('c=0 must be symmetric')
        if c > 0 and asym < 1e-6:
            FAILS.append('c=%d must be non-symmetric' % c)

    print('\nN7  rigid-body content of the assembled operator')
    for c in (0, 1, 2, 3):
        Klin = P.assemble('fbar_%d' % c, n4, t4, m4, g4, v4, fc4, pt4, props,
                          0.0, False, None).toarray()
        Ksym = 0.5 * (Klin + Klin.T)
        ev = np.linalg.eigvalsh(Ksym)
        scale = ev.max()
        nz = int((np.abs(ev) < 1e-9 * scale).sum())
        print('       c=%d  zero modes = %d (expect 6)   '
              'lambda_7/lambda_max = %.3e' % (c, nz, ev[6] / scale))
        if nz != 6:
            FAILS.append('c=%d has %d zero modes' % (c, nz))


# ---------------------------------------------------------------------------
# T  consistent tangent and quadratic Newton
# ---------------------------------------------------------------------------

def fd_tangent(nl: P.FbarNL, u: np.ndarray, h: float = 1e-6,
               ndof: Optional[int] = None) -> np.ndarray:
    """Central-difference tangent of one `FbarNL` force.

    Args:
        nl: The finite-strain element.
        u: The configuration.
        h: Difference step.
        ndof: Number of unknowns, ``3 * nl.nn`` when omitted.

    Returns:
        The dense (ndof, ndof) finite-difference Jacobian.
    """
    ndof = ndof or 3 * nl.nn
    J = np.empty((ndof, ndof))
    e = np.zeros(ndof)
    for j in range(ndof):
        e[j] = h
        J[:, j] = (nl.force(u + e) - nl.force(u - e)) / (2.0 * h)
        e[j] = 0.0
    return J


def check_tangent() -> None:
    """Run the consistent-tangent checks (T1--T3)."""
    n4, t4, m4 = P.mesh_box(4, 0.375, 0.625, 0.3, geom=P.GEOM['sphere'])
    g4, v4 = P.grads(n4, t4)
    fc4, pt4 = P.topology(t4, m4)
    nd = 3 * len(n4)
    props = {0: ICE, 1: (500.0 * GB, GB)}
    props100 = {0: ICE, 1: (100.0 * GB, GB)}

    print('\nT1  tangent(0) is the small-strain operator, to round-off')
    for c in (0, 1, 2, 3):
        nl = P.FbarNL(n4, t4, m4, g4, v4, props, c)
        Klin = P.assemble('fbar_%d' % c, n4, t4, m4, g4, v4, fc4, pt4,
                          props, 0.0, False, None).toarray()
        dK = nl.tangent(np.zeros(nd)).toarray() - Klin
        chk('c=%d  max|tangent(0) - K_lin| / max|K_lin|' % c,
            np.abs(dK).max() / np.abs(Klin).max(), 0.0, 1e-12)

    print('\nT2  central differences at finite strain, K/G = 0, 100, 5000')
    shear = np.zeros(nd)
    shear[2::3] = 0.25 * n4[:, 0]              # gamma_xz = 0.25
    shear[0::3] = 0.25 * n4[:, 2]
    stretch = np.zeros(nd)
    stretch[0::3] = 0.5 * n4[:, 0]             # 50 % stretch
    for tag, u in (('stretch 20%', 0.4 * stretch), ('stretch 50%', stretch),
                   ('two-phase shear', shear)):
        for kg in (0.0, 100.0, 5000.0):
            pr = {0: ICE, 1: (kg * GB, GB)}
            nl = P.FbarNL(n4, t4, m4, g4, v4, pr, 1)
            Kt = nl.tangent(u)
            Jnum = fd_tangent(nl, u)
            err = np.abs(Kt.toarray() - Jnum).max() / np.abs(Jnum).max()
            chk('c=1 %-15s K/G=%-6g  max rel tangent error'
                % (tag, kg), err, 0.0, 1e-7)

    print('\nT3  quadratic Newton with the consistent tangent')
    for macro, eps, ninc in (('xx', 0.2, 10), ('xx', 0.5, 20),
                             ('xz', 0.2, 10), ('xz', 0.5, 20)):
        hist: List[float] = []
        cn, fl, it, rn = P.run_nl('fbar_1', 4, (0.375, 0.625), props100,
                                  eps=eps, jitter=0.3, geomname='sphere',
                                  consistent=True, history=hist, ninc=ninc,
                                  macro=macro, maxit=120)
        start = max(i for i, r in enumerate(hist) if r >= 1.0)
        leg = hist[start + 1:]
        tail = '  '.join('%.1e' % r for r in leg[-5:])
        nquad = sum(1 for a, b in zip(leg[:-1], leg[1:])
                    if a > 0 and b <= a ** 1.5)
        label = 'stretch' if macro == 'xx' else 'shear  '
        print('   %s eps=%.2f  its=%d  C1111=%.6e  |r|=%.1e'
              % (label, eps, it, cn, rn))
        print('      last increment |r| tail: %s   superlinear steps %d'
              % (tail, nquad))
        chk('%s eps=%.2f  converged' % (label.strip(), eps), rn, 0.0, 1e-10)
        chk('%s eps=%.2f  superlinear steps in the last leg (>= 2)'
            % (label.strip(), eps), min(nquad, 2), 2.0, 0.0)
    print('   K/G uniformity of the tangent (stretch eps=0.2)')
    for kg in (0.0, 10.0, 100.0):
        pr = {0: ICE, 1: (kg * GB, GB)}
        hist = []
        cn, fl, it, rn = P.run_nl('fbar_1', 4, (0.375, 0.625), pr, eps=0.2,
                                  jitter=0.3, geomname='sphere',
                                  consistent=True, history=hist, ninc=10,
                                  maxit=120)
        chk('K/G=%-6g  eps=0.2  converged' % kg, rn, 0.0, 1e-10)
        print('      K/G=%-6g  its=%d  |r| tail: %s'
              % (kg, it, '  '.join('%.1e' % r for r in hist[-4:])))


# ---------------------------------------------------------------------------
# S  spurious-mode census
# ---------------------------------------------------------------------------

def census(scheme: str, n: int, K: float, G: float, jitter: float = 0.3,
           stab: float = 0.0, bubble: bool = False, nev: int = 30,
           two_phase: bool = False, ice: Optional[Tuple[float, float]] = None,
           radius: float = 0.30) -> Optional[List[Tuple[float, float]]]:
    """Lowest modes of a clamped two-phase cube and their volumetric content.

    A homogeneous cube cannot show the instability: every node patch reaches
    the clamped boundary, so the null mode has nowhere to live.  The soft
    phase has to be bounded by the stiff one.

    Args:
        scheme: Smoothing scheme name for ``smoothing_proto.assemble``.
        n: Cells per edge.
        K: Bulk modulus of the stiff phase.
        G: Shear modulus.
        jitter: Mesh jitter.
        stab: Stabilisation parameter.
        bubble: Use the subcell bubble gradients.
        nev: Number of eigenpairs.
        two_phase: Embed a brine sphere in ice.
        ice: Properties of the soft phase.
        radius: Sphere radius.

    Returns:
        A list of (eigenvalue, volumetric ratio) pairs, or None when the
        free system is too small.
    """
    if two_phase:
        nodes, tets, mat = S.mesh_box(n, 0.0, radius, jitter,
                                      geom=S.GEOM['sphere'])
    else:
        nodes, tets, mat = S.mesh_box(n, 0.0, 0.0, jitter, geom=S.GEOM['slab'])
        mat[:] = 0
    g, vol = S.grads(nodes, tets)
    faces, patch = S.topology(tets, mat)
    sbg = S.subcell_bubble_grads(nodes, tets, g) if bubble else None
    props = {0: (ice if two_phase else (K, G)), 1: (K, G)}
    Kg = S.assemble(scheme, nodes, tets, mat, g, vol, faces, patch,
                    props, stab, bubble, sbg)

    tol = 1e-9
    onb = ((nodes <= tol) | (nodes >= 1.0 - tol)).any(axis=1)
    nn = len(nodes)
    free = np.ones(Kg.shape[0], dtype=bool)
    free[:3 * nn][np.repeat(onb, 3)] = False
    Kf = Kg[free][:, free].tocsc()
    if Kf.shape[0] <= nev + 2:
        return None
    w, V = spl.eigsh(Kf, k=nev, sigma=0.0, which='LM')
    order = np.argsort(w)
    w, V = w[order], V[:, order]

    # volumetric content of each mode
    B = [S.bmat(g[e]) for e in range(len(tets))]
    idx = np.flatnonzero(free)
    out: List[Tuple[float, float]] = []
    for m in range(len(w)):
        u = np.zeros(Kg.shape[0])
        u[idx] = V[:, m]
        num = den = 0.0
        for e, t in enumerate(tets):
            if two_phase and mat[e] != 1:
                continue                      # only the soft phase can host it
            eps = B[e] @ u[S.dofs_of(t)]
            dv = eps[0] + eps[1] + eps[2]
            ee = (eps[0] ** 2 + eps[1] ** 2 + eps[2] ** 2
                  + 0.5 * (eps[3] ** 2 + eps[4] ** 2 + eps[5] ** 2))
            num += vol[e] * dv * dv
            den += vol[e] * ee
        out.append((w[m], num / den if den else 0.0))
    return out


def check_stability(n: int = 10) -> None:
    """Run the spurious-mode census (informational table).

    Args:
        n: Cells per edge of the clamped cube.
    """
    G = 4.4e5
    ice = (9.37e9 / (3 * (1 - 2 * 0.33)), 9.37e9 / (2 * 1.33))
    arms = [('c3d4', 0.0, False, 'c3d4'),
            ('ns_vol', 0.0, False, 'ns_vol (NS-FEM) s=0'),
            ('ns_vol', 0.07, False, 'ns_vol s=0.07'),
            ('fs_ns', 0.0, False, 'fs_ns (FS/NS-FEM)')]
    print('\nS   two-phase clamped cube: brine sphere r=0.30 in ice, n=%d, '
          'jitter 0.3' % n)
    print('    r is measured over the SOFT PHASE only.\n')
    for ratio in (50, 500, 5000):
        print('  K/G = %-6d (nu = %.5f)'
              % (ratio, (3 * ratio - 2) / (2 * (3 * ratio + 1))))
        print('     %-24s %12s %8s %8s %10s'
              % ('scheme', 'lambda_1/G', 'r_1', 'n_bad', 'max r'))
        for sc, stab, bub, tag in arms:
            c = census(sc, n, ratio * G, G, stab=stab, bubble=bub,
                       two_phase=True, ice=ice)
            if c is None:
                continue
            lam1, r1 = c[0]
            bad = [r for _l, r in c if r > 0.3]
            print('     %-24s %12.4e %8.3f %8d %10.3f'
                  % (tag, lam1 / G, r1, len(bad), max(r for _l, r in c)))
        print()


# ---------------------------------------------------------------------------
# D  deck generator walk; A  assembled-structure audit
# ---------------------------------------------------------------------------

def read_deck(path: str):
    """Elements by type, declarations, declaration order and format errors."""
    els, decl, order, badline = defaultdict(dict), {}, [], []
    elsetof = {}
    mode, cur, curtype, pend, nwant = None, None, None, None, 0
    for lno, raw in enumerate(open(path), 1):
        s = raw.strip()
        if not s or s.startswith('**'):
            continue
        if s.startswith('*'):
            u = s.upper().replace(' ', '')
            mode = None
            if u.startswith('*USERELEMENT'):
                t = next(p.split('=')[1] for p in s.split(',')
                         if p.strip().upper().startswith('TYPE='))
                n = int(next(p.split('=')[1] for p in s.split(',')
                             if p.strip().upper().startswith('NODES=')))
                decl[t.strip()] = n
                order.append(('decl', t.strip()))
            elif u.startswith('*ELEMENT'):
                curtype = next(p.split('=')[1] for p in s.split(',')
                               if p.strip().upper().startswith('TYPE=')).strip()
                cur = next((p.split('=')[1] for p in s.split(',')
                            if p.strip().upper().startswith('ELSET=')),
                           'ALL').strip()
                order.append(('use', curtype))
                elsetof[curtype] = cur
                mode = 'e'
                pend, nwant = None, decl.get(curtype, 4)
            continue
        if mode != 'e':
            continue
        f = [x.strip() for x in s.split(',')]
        if len(f) > 16:
            badline.append((lno, 'more than 16 fields'))
        if f and f[-1] == '':
            badline.append((lno, 'trailing comma'))
        f = [x for x in f if x]
        if pend is None:
            pend = [int(x) for x in f]
        else:
            pend += [int(x) for x in f]
        if len(pend) >= nwant + 1:
            els[curtype][pend[0]] = pend[1:nwant + 1]
            pend = None
    return els, decl, order, badline, elsetof


def walk_deck(tets, mats, imat, edge, ncyc) -> Set[int]:
    """u3vol.f's walk, as a set of elements."""
    at = defaultdict(list)
    for e, t in tets.items():
        if mats[e] != imat:
            continue
        for a in t:
            at[a].append(e)
    na, nb = edge
    cur = {e for e in at[na] if nb in tets[e]}
    for _ in range(ncyc):
        nxt = set()
        for e in cur:
            for k in tets[e]:
                nxt.update(at[k])
        cur = nxt
    return cur


def check_deck(orig: str, gen: str, ncyc: int = 1, nsample: int = 500) -> None:
    """Check that a deck written by fbares.py says what the elements walk.

    Args:
        orig: The original C3D4 deck.
        gen: The generated F-barES deck.
        ncyc: Cycle count used by the generator.
        nsample: Elements sampled in the O(stencil) connectivity walk; the
            format checks are exhaustive regardless.
    """
    print('\nD   deck generator: %s' % os.path.basename(gen))
    els, decl, order, badline, elsetof = read_deck(gen)

    chk_ok('deck: no line over 16 fields, no trailing comma', not badline,
           str(badline[:3]))
    seen = set()
    dup = [e for t in els for e in els[t] if e in seen or seen.add(e)]
    chk_ok('deck: element ids unique', not dup, str(dup[:5]))
    firstuse = {}
    undeclared = []
    for kind, t in order:
        if kind == 'use' and t not in firstuse:
            firstuse[t] = True
            if t.startswith('U') and t not in decl and t != 'U4':
                undeclared.append(t)
    chk_ok('deck: every user type declared before use', not undeclared,
           str(undeclared[:5]))
    bad = [t for t in decl if not t[2:].isalpha() and t[2:]]
    chk_ok('deck: type suffixes are letters only (elements.f digit rule)',
           not bad, str(bad[:5]))

    # material of every U4 tet, from the *SOLID SECTION of its elset
    tets, mats = {}, {}
    esof, matof = {}, {}
    mode, cur = None, None
    for raw in open(gen):
        s = raw.strip()
        if not s or s.startswith('**'):
            continue
        if s.startswith('*'):
            u = s.upper().replace(' ', '')
            mode = None
            if u.startswith('*ELEMENT') and 'TYPE=U4' in u:
                cur = next((p.split('=')[1] for p in s.split(',')
                            if p.strip().upper().startswith('ELSET=')),
                           'ALL').strip()
                mode = 'u4'
            elif u.startswith('*SOLIDSECTION'):
                es = next(p.split('=')[1] for p in s.split(',')
                          if p.strip().upper().startswith('ELSET=')).strip()
                mt = next(p.split('=')[1] for p in s.split(',')
                          if p.strip().upper().startswith('MATERIAL=')).strip()
                matof[es] = mt
            continue
        if mode == 'u4':
            f = [x.strip() for x in s.split(',') if x.strip()]
            if len(f) >= 5:
                tets[int(f[0])] = [int(x) for x in f[1:5]]
                esof[int(f[0])] = cur
    for e in tets:
        mats[e] = matof[esof[e]]
    chk_ok('deck: U4 tets carried over from the original', len(tets) > 0,
           '%d found' % len(tets))

    # every U2 / U3 element, against an independent walk
    nd2 = nd3 = 0
    bad2, bad3, badr = [], [], []
    LET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    rnd = random.Random(11)
    for t, table in els.items():
        if not (t.startswith('U2') or t.startswith('U3')):
            continue
        items = list(table.items())
        if nsample and len(items) > nsample:
            items = rnd.sample(items, nsample)
        for eid, conn in items:
            na, nb = conn[0], conn[1]
            # the element's OWN material, from its *SOLID SECTION
            imat = matof[elsetof[t]]
            cyc = 0 if t.startswith('U2') else ncyc
            got = walk_deck(tets, mats, imat, (na, nb), cyc)
            want = set()
            for e in got:
                want.update(tets[e])
            have = set(conn)
            if t.startswith('U2'):
                nd2 += 1
                if have != want:
                    bad2.append((eid, sorted(want - have), sorted(have - want)))
            else:
                nd3 += 1
                if have != want:
                    bad3.append((eid, sorted(want - have), sorted(have - want)))
                # RING FIRST, AND THE LABEL AGREES WITH IT
                nr = LET.find(t[2]) + 1
                ringset = set()
                for e in walk_deck(tets, mats, imat, (na, nb), 0):
                    ringset.update(tets[e])
                if nr < 1 or nr != len(ringset) or set(conn[:nr]) != ringset:
                    badr.append((eid, t, nr, len(ringset)))
    chk_ok('U2 connectivity == the ring u2edge.f walks (%d elements)' % nd2,
           not bad2, str(bad2[:2]))
    chk_ok('U3 connectivity == the E A^c support u3vol.f walks (%d elements)'
           % nd3, not bad3, str(bad3[:2]))
    chk_ok('U3 ring is first in the connectivity and matches lakon(3:3)',
           not badr, str(badr[:2]))

    # the edge nodes must be konl(1), konl(2) and must share a tet
    at = defaultdict(list)
    for e, t4 in tets.items():
        for a in t4:
            at[a].append(e)
    bade = []
    for t, table in els.items():
        if not (t.startswith('U2') or t.startswith('U3')):
            continue
        items = list(table.items())
        if nsample and len(items) > nsample:
            items = rnd.sample(items, nsample)
        for eid, conn in items:
            if not any(conn[1] in tets[e] for e in at[conn[0]]):
                bade.append(eid)
    chk_ok('konl(1),konl(2) are a real mesh edge', not bade, str(bade[:5]))


def check_structure(deck: str) -> None:
    """Print the assembled-structure audit of a converted deck.

    Args:
        deck: The generated F-barES deck.
    """
    print('\nA   assembled-structure audit: %s' % os.path.basename(deck))
    types = {}
    for line in open(deck):
        if line.startswith('*USER ELEMENT'):
            m = re.search(r'TYPE=(\w+).*NODES=(\d+)', line)
            if m and m.group(1).startswith('U3'):
                nm = m.group(1)
                types[nm] = (ord(nm[2]) - ord('A') + 1, int(m.group(2)))

    conn = {}
    cur = None
    toks = []

    def flush():
        if cur and toks:
            conn.setdefault(cur, []).extend(toks)

    for line in open(deck):
        if line[0] == '*':
            flush()
            toks = []
            m = re.match(r'\*ELEMENT,TYPE=(\w+)', line.strip())
            cur = m.group(1) if m and m.group(1) in types else None
            continue
        if cur:
            toks.extend(line.split(','))
    flush()

    M = 1 << 21
    keys = []
    nel = transient = plan_alloc = 0
    for t, tk in conn.items():
        r, s = types[t]
        a = np.array(tk, dtype=np.int64).reshape(-1, s + 1)
        c = a.shape[0]
        nel += c
        rr, ss = 3 * r, 3 * s
        transient += c * (rr * ss - rr * (rr - 1) // 2)
        plan_alloc += c * rr * ss
        node = a[:, 1:]
        ring = node[:, :r]
        p = np.broadcast_to(ring[:, :, None], (c, r, s)).reshape(-1)
        q = np.broadcast_to(node[:, None, :], (c, r, s)).reshape(-1)
        keys.append(np.unique(np.minimum(p, q) * M + np.maximum(p, q)))
        del a, node, ring, p, q

    allk = np.unique(np.concatenate(keys))
    lo, hi = allk // M, allk % M
    ndist, nself = int((lo != hi).sum()), int((lo == hi).sum())
    final = 9 * ndist + 3 * nself
    print('    U3 elements             : %d' % nel)
    print('    transient insert() calls: %.4e' % transient)
    print('    plan preallocation      : %.4e  (%.2fx the transient)'
          % (plan_alloc, plan_alloc / transient))
    print('    true final structure    : %.4e  (%.1fx smaller than transient)'
          % (final, transient / final))


def make_decks(tmpdir: str, n: int) -> Tuple[str, str]:
    """Generate a small original deck and its F-barES conversion.

    Args:
        tmpdir: Directory for the generated files.
        n: Mesh size for make_slabconv.py.

    Returns:
        The paths (original C3D4 deck, generated F-barES deck).
    """
    stem = os.path.join(tmpdir, 'm')
    subprocess.run([sys.executable, os.path.join(HERE, 'make_slabconv.py'),
                    stem, str(n), 'und', '500', '0.3', 'one', 'x', 'sym'],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT)
    orig = stem + '_ccx.inp'
    gen = os.path.join(tmpdir, 'mf_ccx.inp')
    subprocess.run([sys.executable, os.path.join(HERE, os.pardir, 'fbares.py'),
                    orig, gen, '--elset', 'Sphere_Only', '--cycles', '1'],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT)
    return orig, gen


def main() -> int:
    """Run the checks selected on the command line.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true',
                    help='operators, walk and deck checks only; skip N, T, S')
    ap.add_argument('--skip-stability', action='store_true')
    ap.add_argument('--skip-tangent', action='store_true')
    ap.add_argument('--n', type=int, default=10,
                    help='mesh size for the stability census')
    ap.add_argument('--deck', nargs='+', metavar='FILE',
                    help='original and generated decks, plus optional cycles')
    ap.add_argument('--sample', type=int, default=500,
                    help='elements sampled in the deck walk')
    a = ap.parse_args()

    check_operators()
    check_walk()
    if not a.quick:
        check_force()
        if not a.skip_tangent:
            check_tangent()
        if not a.skip_stability:
            check_stability(a.n)

    if a.deck:
        orig, gen = a.deck[0], a.deck[1]
        ncyc = int(a.deck[2]) if len(a.deck) > 2 else 1
        check_deck(orig, gen, ncyc, a.sample)
        check_structure(gen)
    else:
        try:
            with tempfile.TemporaryDirectory(prefix='fbar_verify_') as tmp:
                orig, gen = make_decks(tmp, 10)
                check_deck(orig, gen, 1, a.sample)
                check_structure(gen)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            print('\nD   deck checks skipped: %s' % exc)
            FAILS.append('deck checks skipped')

    print('\n%s' % ('all checks passed' if not FAILS
                    else 'FAILED: ' + ', '.join(FAILS)))
    return 1 if FAILS else 0


if __name__ == '__main__':
    sys.exit(main())
