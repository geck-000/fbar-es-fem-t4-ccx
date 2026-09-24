"""Verify the Fortran finite-strain tangent u3nltan against the prototype.

`nl_tan_driver.f` calls `u3nltan` once per U3 edge element on a small
single-material mesh and writes the global (row, column, value) triplets
of the assembled matrix.  The check is twofold:

  * at u = 0 the assembled tangent must equal the small-strain operator
    `assemble('fbar_1')` -- the linear limit of the eq. (17) chain;
  * at a finite displacement field it must equal `FbarNL.tangent(u)`, the
    prototype operator chain verified against central differences in
    `verify_element.py`.

The driver and the element sources live in this repository; gfortran is
the only external requirement.

    python3 elements_ccx/tests/verify_nltan.py
"""
from __future__ import annotations

import os
import subprocess
import sys
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as P  # noqa: E402
import verify_convergence as vc  # noqa: E402
import mms_ccx  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', '..'))
DRIVER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'nl_tan_driver.f')
CCX_SRC = os.path.expanduser('~/ccx_sym/CalculiX/ccx_2.23/src')
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


def u3_elements(nodes: np.ndarray, tets: np.ndarray, g: np.ndarray,
                vol: np.ndarray
                ) -> List[Tuple[int, int, List[int]]]:
    """Build the U3 element list in the generator's ring-first order.

    Args:
        nodes: Node coordinates.
        tets: Element connectivity.
        g: Shape-function gradients.
        vol: Element volumes.

    Returns:
        One tuple (nring, nope, konl) per edge, konl zero-based.
    """
    mat = np.ones(len(tets), dtype=int)
    sel, idx, ekeys, edges, Vh, E, S, R = P.fbar_es_operators(
        nodes, tets, mat, g, vol, 1, 1)
    Sc = S.tocoo()
    at: Dict[int, List[int]] = {}
    for e, t in enumerate(tets):
        for nd in t:
            at.setdefault(int(nd), []).append(e)
    out = []
    for h, key in enumerate(ekeys):
        na, nb = key
        ring = [e for e in at[int(na)] if e in at[int(nb)]]
        ringnodes = set()
        for e in ring:
            ringnodes.update(int(x) for x in tets[e])
        chain = Sc.col[Sc.row == h]
        supp = set()
        for e in chain:
            supp.update(int(x) for x in tets[e])
        assert ringnodes <= supp, 'ring not contained in the chain support'
        inner = [na, nb] + sorted(ringnodes - {na, nb})
        outer = sorted(supp - ringnodes)
        konl = inner + outer
        out.append((len(ringnodes), len(supp), konl))
    return out


def assemble(nn: int, elems: List[Tuple[int, int, List[int]]],
             triplets: List[Tuple[int, int, float]]) -> np.ndarray:
    """Assemble the global matrix from the driver triplets.

    Args:
        nn: Number of nodes.
        elems: The element list (unused, kept for symmetry with the driver).
        triplets: (row, column, value), one-based.

    Returns:
        The dense (3 nn, 3 nn) matrix.
    """
    K = np.zeros((3 * nn, 3 * nn))
    for gr, gc, val in triplets:
        K[gr - 1, gc - 1] += val
    return K


def run_case(amplitude: float, out: str = '/tmp/nltan') -> Tuple[np.ndarray,
                                                                 np.ndarray]:
    """Run one driver case and return (K_fortran, u).

    Args:
        amplitude: Amplitude of the random displacement field.
        out: Working directory.

    Returns:
        The pair (assembled Fortran tangent, displacement vector).
    """
    n = 3
    nodes, tets, mat = P.mesh_box(n, 0.0, 1.0, 0.0, geom=vc.geom_single)
    props = {1: vc.iso(100.0)}
    g, vol = P.grads(nodes, tets)
    elems = u3_elements(nodes, tets, g, vol)
    rng = np.random.default_rng(1)
    u = rng.standard_normal(3 * len(nodes)) * amplitude
    U = u.reshape(-1, 3)
    lines = ['%d %d' % (len(nodes), len(tets))]
    for p in nodes:
        lines.append('%.17e %.17e %.17e' % tuple(p))
    for t in tets:
        lines.append(' '.join(str(int(x) + 1) for x in t))
    lines.append('2 1')
    lines.append('2 1')
    lines.append('%.17e %.17e' % mms_ccx.e_nu(*props[1]))
    for i in range(len(nodes)):
        lines.append('%.17e %.17e %.17e' % (U[i, 0], U[i, 1], U[i, 2]))
    lines.append(str(len(elems)))
    for nring, nope, konl in elems:
        lines.append('%d %d' % (nope, nring))
        lines.append(' '.join(str(x + 1) for x in konl))
    os.makedirs(out, exist_ok=True)
    inp = '\n'.join(lines) + '\n'
    exe = os.path.join(out, 'tldrv')
    subprocess.run(['gfortran', '-O2', '-o', exe, DRIVER,
                    f'{REPO}/elements_ccx/u3nltan.f',
                    f'{REPO}/elements_ccx/u3nl.f',
                    f'{REPO}/elements_ccx/u3vol.f',
                    f'{CCX_SRC}/shape4tet.f',
                    f'{REPO}/elements_ccx/fbar_lock.c'], check=True)
    res = subprocess.run([exe], input=inp, text=True, capture_output=True,
                         check=True)
    triplets = []
    for line in res.stdout.splitlines():
        parts = line.split()
        if len(parts) == 3:
            triplets.append((int(parts[0]), int(parts[1]),
                             float(parts[2])))
    return assemble(len(nodes), elems, triplets), u


def main() -> int:
    """Run both tangent comparisons.

    Returns:
        Process exit status.
    """
    n = 3
    nodes, tets, mat = P.mesh_box(n, 0.0, 1.0, 0.0, geom=vc.geom_single)
    g, vol = P.grads(nodes, tets)
    faces, patch = P.topology(tets, mat)
    props = {1: vc.iso(100.0)}

    print('u = 0: the tangent must be the small-strain operator')
    K0, _ = run_case(0.0, '/tmp/nltan0')
    Klin = P.assemble('fbar_1', nodes, tets, mat, g, vol, faces, patch,
                      props).toarray()
    chk('max|u3nltan(0) - K_lin| / max|K_lin|',
        np.abs(K0 - Klin).max() / np.abs(Klin).max(), 0.0, 1e-12)

    print('\nu = random, amplitude 5e-2: against FbarNL.tangent')
    K, u = run_case(5e-2, '/tmp/nltan1')
    nl = P.FbarNL(nodes, tets, mat, g, vol, props, 1)
    Kp = nl.tangent(u).toarray()
    chk('max|u3nltan - FbarNL.tangent| / max|K|',
        np.abs(K - Kp).max() / np.abs(Kp).max(), 0.0, 1e-10)

    print('\n%s' % ('all checks passed' if not fails
                    else 'FAILED: ' + ', '.join(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
