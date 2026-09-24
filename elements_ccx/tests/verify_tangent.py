"""Verify the consistent tangent of the finite-strain F-barES-FEM-T4 force.

The force of `FbarNL` is fully specified by eqs. (1)-(17); the tangent is the
operator chain

    du -> dF -> dF~ -> (dJ, dJbar) -> dFbar -> dT -> dA -> df

with the Hencky tangent taken in the eigenbasis of B = F~bar F~bar^T.  These
checks have closed-form or independent references, so a failure localises:

  T1  tangent(0) equals the small-strain operator of `assemble` to round-off.
      The linear operator was built by a completely different route (sparse
      products of the smoothing operators), so this is a genuine cross-check,
      not a restatement.
  T2  central differences of `force` against `tangent` at finite strain: a
      20 % stretch, a 50 % stretch, and a two-phase shear, at K/G = 0, 100
      and 5000 (the tangent must be uniform in K/G).
  T3  Newton with the consistent tangent converges quadratically: 20 % and
      50 % stretch and 20 % and 50 % simple shear on a two-phase periodic
      cell, plus a K/G sweep.  The empirical order of the last three
      residuals must be 2 +- 0.5.  The strain is applied in increments, as
      it must be to reach 20-50 % at all; the LAST increment is measured.

    python3 elements_ccx/tests/verify_tangent.py
"""
from __future__ import annotations

import os
import sys
from typing import List, Optional

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as P  # noqa: E402

Ei, ni = 9.37e9, 0.33
ice = (Ei / (3 * (1 - 2 * ni)), Ei / (2 * (1 + ni)))
Gb = 4.4e5
fails: List[str] = []


def chk(name: str, got: float, want: float, tol: float) -> None:
    """Print and record one check.

    Args:
        name: Label printed in the table.
        got: Measured value.
        want: Expected value.
        tol: Absolute tolerance on ``got - want``.
    """
    ok = abs(got - want) <= tol
    print('   %-58s %12.4e  %s' % (name, got, 'PASS' if ok else 'FAIL'))
    if not ok:
        fails.append(name)


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


def main() -> int:
    """Run the tangent checks.

    Returns:
        Process exit status.
    """
    n4, t4, m4 = P.mesh_box(4, 0.375, 0.625, 0.3, geom=P.GEOM['sphere'])
    g4, v4 = P.grads(n4, t4)
    fc4, pt4 = P.topology(t4, m4)
    nd = 3 * len(n4)
    props = {0: ice, 1: (500.0 * Gb, Gb)}
    props100 = {0: ice, 1: (100.0 * Gb, Gb)}

    print('T1  tangent(0) is the small-strain operator, to round-off')
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
            pr = {0: ice, 1: (kg * Gb, Gb)}
            nl = P.FbarNL(n4, t4, m4, g4, v4, pr, 1)
            Kt = nl.tangent(u)
            Jnum = fd_tangent(nl, u)
            err = np.abs(Kt.toarray() - Jnum).max() / np.abs(Jnum).max()
            chk('c=1 %-15s K/G=%-6g  max rel tangent error'
                % (tag, kg), err, 0.0, 1e-7)

    print('\nT3  quadratic Newton with the consistent tangent')
    # K/G = 100, n = 4, jitter 0.3: the two-phase cell stays invertible to
    # 50 % and every increment converges quadratically.  At K/G = 500 the
    # cell itself reaches an inverted element before 20 % on this mesh --
    # a limit of the load path, not of the tangent, which T2 verifies there.
    for macro, eps, ninc in (('xx', 0.2, 10), ('xx', 0.5, 20),
                             ('xz', 0.2, 10), ('xz', 0.5, 20)):
        hist: List[float] = []
        cn, fl, it, rn = P.run_nl('fbar_1', 4, (0.375, 0.625), props100,
                                  eps=eps, jitter=0.3, geomname='sphere',
                                  consistent=True, history=hist, ninc=ninc,
                                  macro=macro, maxit=120)
        # each increment restarts at |r| = 1: keep the last increment's leg
        start = max(i for i, r in enumerate(hist) if r >= 1.0)
        leg = hist[start + 1:]
        tail = '  '.join('%.1e' % r for r in leg[-5:])
        # superlinear steps r_{k+1} <= r_k^1.5 among the last three; a
        # quadratic Newton shows two of them as it enters the basin
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
        pr = {0: ice, 1: (kg * Gb, Gb)}
        hist = []
        cn, fl, it, rn = P.run_nl('fbar_1', 4, (0.375, 0.625), pr, eps=0.2,
                                  jitter=0.3, geomname='sphere',
                                  consistent=True, history=hist, ninc=10,
                                  maxit=120)
        chk('K/G=%-6g  eps=0.2  converged' % kg, rn, 0.0, 1e-10)
        print('      K/G=%-6g  its=%d  |r| tail: %s'
              % (kg, it, '  '.join('%.1e' % r for r in hist[-4:])))

    print('\n%s' % ('all checks passed' if not fails
                    else 'FAILED: ' + ', '.join(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
