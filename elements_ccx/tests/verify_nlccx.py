"""Finite-strain U3 force through CalculiX (NLGEOM), against the prototype.

The boundary of a structured cube is prescribed with the divergence-free
manufactured field of `verify_convergence`, the interior is left free, and
the step is `*STEP, NLGEOM`.  Under NLGEOM the U3 element returns the
finite-strain force of eqs. (1)-(17) through `u3nlforce` (and U2 returns no
separate force, because the finite-strain force already carries both
halves).  At equilibrium the reactions at the prescribed degrees of freedom
are the internal force there, so they can be compared node by node with
`FbarNL.force` evaluated at the CalculiX solution -- an independent
implementation of the same equations.

    python3 elements_ccx/tests/verify_nlccx.py --ccx <binary>

The iteration matrix is still the small-strain stiffness (the consistent
tangent is the next stage), so the test runs at a small boundary amplitude
where the modified Newton converges; the check is on the force, not on the
convergence rate.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from typing import Dict, Optional, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as proto  # noqa: E402
import verify_convergence as vc  # noqa: E402
import mms_ccx  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', '..'))
FBARES = os.path.join(REPO, 'elements_ccx', 'fbares.py')
fails = []


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


def read_dat_block(path: str, key: str) -> Optional[np.ndarray]:
    """Read the last block of a CalculiX .dat file.

    Args:
        path: The .dat file.
        key: Lower-case block header prefix ('displacements' or 'forces').

    Returns:
        The values ordered by node id as an (n_node, 3) array, or None.
    """
    values: Dict[int, Tuple[float, float, float]] = {}
    active = started = False
    with open(path, errors='replace') as handle:
        for line in handle:
            s = line.strip()
            low = s.lower()
            if low.startswith(key):
                active, started, values = True, False, {}
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
                vals = (float(parts[1]), float(parts[2]), float(parts[3]))
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


def run_bvp(outdir: str, n: int, eps: float, ccx: str,
            inc: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """Generate, solve and read back the NLGEOM boundary-value problem.

    Args:
        outdir: Directory for the deck and results.
        n: Elements per edge.
        eps: Amplitude of the manufactured boundary field.
        ccx: CalculiX binary.
        inc: Initial pseudo-time increment for *STATIC, or 0 for the
            default single step.

    Returns:
        The pair (displacements, reactions), each (n_node, 3).
    """
    nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0,
                                      geom=vc.geom_single)
    props = {1: vc.iso(100.0)}
    fixed = vc.boundary_mask(nodes)
    ustar = vc.u_star(nodes) * eps
    lines = ['*HEADING', ' finite-strain U3 verification', '*NODE']
    for i, p in enumerate(nodes):
        lines.append('%d, %.12e, %.12e, %.12e' % (i + 1, p[0], p[1], p[2]))
    lines.append('*ELEMENT,TYPE=C3D4,ELSET=BODY')
    for e, t in enumerate(tets):
        lines.append('%d, %s' % (e + 1, ', '.join(str(int(x) + 1)
                                                   for x in t)))
    lines += mms_ccx.nset(range(1, len(nodes) + 1), 'ALLNODES')
    lines += ['*SOLID SECTION,ELSET=BODY,MATERIAL=Mat_Body',
              '*MATERIAL,NAME=Mat_Body', '*ELASTIC']
    ee, nu = mms_ccx.e_nu(*props[1])
    lines.append('%.12e, %.12e' % (ee, nu))
    lines += ['*STEP, NLGEOM', '*STATIC, SOLVER=PARDISO']
    if inc > 0.0:
        lines.append('%.6g, 1.0, 1e-6, 1.0' % inc)
    lines += ['*BOUNDARY']
    for i in np.flatnonzero(fixed):
        for c in range(3):
            lines.append('%d, %d, %d, %.12e'
                         % (i + 1, c + 1, c + 1, ustar[i, c]))
    lines += ['*NODE PRINT,NSET=ALLNODES', 'U,RF', '*END STEP']
    os.makedirs(outdir, exist_ok=True)
    src = os.path.join(outdir, 'm.inp')
    with open(src, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    subprocess.run([sys.executable, FBARES, src,
                    os.path.join(outdir, 'mf.inp'), '--cycles', '1',
                    '--elset', 'BODY', '--nlgeom'], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    env = dict(os.environ)
    with open(os.path.join(outdir, 'solve.log'), 'w') as log:
        subprocess.run([ccx, 'mf'], cwd=outdir, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    u = read_dat_block(os.path.join(outdir, 'mf.dat'), 'displacements')
    rf = read_dat_block(os.path.join(outdir, 'mf.dat'), 'forces')
    if u is None or rf is None:
        raise SystemExit('verify_nlccx: no U/RF written; see %s'
                         % os.path.join(outdir, 'solve.log'))
    return u, rf


def run_uniaxial(outdir: str, n: int, eps: float, ccx: str,
                 inc: float = 0.1) -> Tuple[np.ndarray, np.ndarray]:
    """Generate, solve and read a uniaxial-stretch NLGEOM problem.

    Rollers on the three symmetry faces and u_x = eps on x = 1.  This is
    the 20-50 % stretch path: the lateral faces are free, so the exact
    solution is a smooth uniaxial state and the load path is reachable
    where the multi-axial MMS boundary field tangles the mesh.

    Args:
        outdir: Directory for the deck and results.
        n: Elements per edge.
        eps: Applied stretch.
        ccx: CalculiX binary.
        inc: Initial *STATIC pseudo-time increment.

    Returns:
        The pair (displacements, reactions), each (n_node, 3).
    """
    nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0,
                                      geom=vc.geom_single)
    props = {1: vc.iso(100.0)}
    tol = 1e-9
    lines = ['*HEADING', ' uniaxial stretch', '*NODE']
    for i, p in enumerate(nodes):
        lines.append('%d, %.12e, %.12e, %.12e' % (i + 1, p[0], p[1], p[2]))
    lines.append('*ELEMENT,TYPE=C3D4,ELSET=BODY')
    for e, t in enumerate(tets):
        lines.append('%d, %s' % (e + 1, ', '.join(str(int(x) + 1)
                                                   for x in t)))
    lines += mms_ccx.nset(range(1, len(nodes) + 1), 'ALLNODES')
    lines += ['*SOLID SECTION,ELSET=BODY,MATERIAL=Mat_Body',
              '*MATERIAL,NAME=Mat_Body', '*ELASTIC']
    ee, nu = mms_ccx.e_nu(*props[1])
    lines.append('%.12e, %.12e' % (ee, nu))
    lines += ['*STEP, NLGEOM', '*STATIC, SOLVER=PARDISO']
    if inc > 0.0:
        lines.append('%.6g, 1.0, 1e-6, 1.0' % inc)
    lines += ['*BOUNDARY']
    for i, p in enumerate(nodes):
        if p[0] <= tol:
            lines.append('%d, 1, 1, 0.0' % (i + 1))
        if p[1] <= tol:
            lines.append('%d, 2, 2, 0.0' % (i + 1))
        if p[2] <= tol:
            lines.append('%d, 3, 3, 0.0' % (i + 1))
        if p[0] >= 1.0 - tol:
            lines.append('%d, 1, 1, %.12e' % (i + 1, eps))
    lines += ['*NODE PRINT,NSET=ALLNODES', 'U,RF', '*END STEP']
    os.makedirs(outdir, exist_ok=True)
    src = os.path.join(outdir, 'm.inp')
    with open(src, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    subprocess.run([sys.executable, FBARES, src,
                    os.path.join(outdir, 'mf.inp'), '--cycles', '1',
                    '--elset', 'BODY', '--nlgeom'], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    env = dict(os.environ)
    with open(os.path.join(outdir, 'solve.log'), 'w') as log:
        subprocess.run([ccx, 'mf'], cwd=outdir, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    u = read_dat_block(os.path.join(outdir, 'mf.dat'), 'displacements')
    rf = read_dat_block(os.path.join(outdir, 'mf.dat'), 'forces')
    if u is None or rf is None:
        raise SystemExit('verify_nlccx: no U/RF written; see %s'
                         % os.path.join(outdir, 'solve.log'))
    return u, rf


def uniaxial(n: int, eps: float, ccx: str, inc: float, out: str) -> None:
    """Check one uniaxial-stretch case.

    Args:
        n: Elements per edge.
        eps: Applied stretch.
        ccx: CalculiX binary.
        inc: Initial increment.
        out: Output directory.
    """
    nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0,
                                      geom=vc.geom_single)
    props = {1: vc.iso(100.0)}
    tol = 1e-9
    g, vol = proto.grads(nodes, tets)
    nl = proto.FbarNL(nodes, tets, mat, g, vol, props, 1)
    u, rf = run_uniaxial(out, n, eps, ccx, inc)
    if not np.isfinite(u).all():
        print('   stretch %.2f: non-finite displacements' % eps)
        fails.append('stretch %.2f finite' % eps)
        return
    fixed = ((nodes <= tol) | (nodes >= 1.0 - tol)).any(axis=1)
    fmask = np.repeat(fixed, 3)
    uflat = u.reshape(-1)
    f = nl.force(uflat)
    scale = max(np.abs(f).max(), 1e-30)
    back = np.zeros((3 * len(nodes),))
    check = np.flatnonzero(~fmask)
    print('   stretch %.2f: max|u_x| %.4f  free residual/max|f| %.3e'
          % (eps, np.abs(u[:, 0]).max(),
             np.abs(f[check]).max() / scale))
    rr = rf.reshape(-1)
    active = np.flatnonzero(fmask & (rr != 0.0))
    chk('stretch %.2f  applied stretch reached' % eps,
        np.abs(u[:, 0]).max(), eps, 1e-6 * eps)
    chk('stretch %.2f  free-dof residual / max|f|' % eps,
        np.abs(f[check]).max() / scale, 0.0, 1e-4)
    chk('stretch %.2f  reactions vs FbarNL.force' % eps,
        np.abs(rr[active] - f[active]).max() / scale, 0.0, 1e-4)


def main() -> int:
    """Run the NLGEOM force check.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--ccx', default=os.environ.get('CCX_MMS', 'ccx_fbar'))
    ap.add_argument('--out', default='out_nlccx')
    ap.add_argument('--n', type=int, default=4)
    ap.add_argument('--uniaxial', action='store_true',
                    help='run the 20/50%% uniaxial-stretch cases instead of '
                         'the MMS boundary-value problem')
    a = ap.parse_args()

    if a.uniaxial:
        for eps in (0.2, 0.5):
            uniaxial(a.n, eps, a.ccx, 0.1,
                     os.path.join(a.out, 'uniax_%g' % eps))
        print('\n%s' % ('all checks passed' if not fails
                        else 'FAILED: ' + ', '.join(fails)))
        return 1 if fails else 0

    n = a.n
    eps = float(os.environ.get('EPS', '1e-3'))
    inc = float(os.environ.get('INC', '0.0'))
    nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0, geom=vc.geom_single)
    props = {1: vc.iso(100.0)}
    fixed = vc.boundary_mask(nodes)
    g, vol = proto.grads(nodes, tets)
    faces, patch = proto.topology(tets, mat)
    K = proto.assemble('fbar_1', nodes, tets, mat, g, vol, faces, patch,
                       props)
    nl = proto.FbarNL(nodes, tets, mat, g, vol, props, 1)
    u, rf = run_bvp(a.out, n, eps, a.ccx, inc)
    uflat = u.reshape(-1)
    f_nl = nl.force(uflat)[np.repeat(fixed, 3)]
    f_lin = (K @ uflat)[np.repeat(fixed, 3)]
    rr = rf.reshape(-1)[np.repeat(fixed, 3)]

    den = np.abs(f_nl).max()
    print('n = %d, eps = %g, %d prescribed dofs (inc %g)'
          % (n, eps, len(rr), inc))
    chk('max|RF - FbarNL.force(u_ccx)| / max|f|', np.abs(rr - f_nl).max()
        / den, 0.0, 1e-5)
    active = np.abs(rr - f_lin).max() / den
    ok = active > 1e-3
    print('   %-58s %12.4e  %s'
          % ('the finite-strain force differs from the linear one', active,
             'PASS' if ok else 'FAIL'))
    if not ok:
        fails.append('finite-strain branch inactive')

    print('\n%s' % ('all checks passed' if not fails
                    else 'FAILED: ' + ', '.join(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
