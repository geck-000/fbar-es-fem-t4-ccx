"""Verify the finite-strain F-barES-FEM-T4 internal force against exact answers.

Every check below has a closed-form expected value, so a failure localises to
one equation rather than to 'the element'.

  N1  f(0) = 0 exactly.  F = I everywhere, so J = J~ = Jbar = 1, H = 0, T = 0.
  N2  rigid translation gives f = 0: the force depends on u only through
      gradients.
  N3  homogeneous uniform stretch.  Every element sees the same F, so every
      smoothing operator in eqs. (1), (6)-(8) returns it unchanged, F~ = Fbar
      = F, and the Hencky stress is analytic.  The affine field must therefore
      be an exact equilibrium state -- residual zero at q = 0, to round-off --
      and C1111 must match the analytic Hencky tangent as eps -> 0.  This is
      the patch test, and it exercises eqs. (1)-(18) end to end.
  N4  the finite-strain element must converge to its own small-strain
      reduction as eps -> 0, at O(eps).  This is what ties run_nl to run.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as P

Ei, ni = 9.37e9, 0.33
ice = (Ei / (3 * (1 - 2 * ni)), Ei / (2 * (1 + ni)))
Gb = 4.4e5
fails = []


def chk(name, got, want, tol):
    ok = abs(got - want) <= tol
    print('   %-56s %12.4e  %s' % (name, got, 'PASS' if ok else 'FAIL'))
    if not ok:
        fails.append(name)


nodes, tets, mat = P.mesh_box(6, 0.375, 0.625, 0.3, geom=P.GEOM['sphere'])
g, vol = P.grads(nodes, tets)
props = {0: ice, 1: (500 * Gb, Gb)}

print('N1  f(0) = 0')
for c in (0, 1, 2, 3):
    nl = P.FbarNL(nodes, tets, mat, g, vol, props, c)
    f = nl.force(np.zeros(3 * len(nodes)))
    chk('c=%d  max|f(0)|' % c, np.abs(f).max(), 0.0, 1e-6)

print('\nN2  rigid translation gives f = 0')
u = np.zeros(3 * len(nodes))
u[0::3] = 0.37; u[1::3] = -0.21; u[2::3] = 0.08
for c in (0, 2):
    nl = P.FbarNL(nodes, tets, mat, g, vol, props, c)
    chk('c=%d  max|f(rigid)|' % c, np.abs(nl.force(u)).max(), 0.0, 1e-6)

print('\nN3  homogeneous uniform stretch: affine field is exact equilibrium')
# ONE material.  The smoothing of eqs. (1) and (6)-(8) is not allowed to cross
# a material interface, so on a two-tag mesh the interface nodes see a
# truncated stencil and carry a genuine (and correct) residual even when the
# two tags share properties.  A patch test needs a single-material mesh.
#
# It also needs FLAT FACES to read a resultant off, and the periodic cell does
# not have them: mesh_box jitters the wrapped images too, so the cell boundary
# is ragged and spills outside [0,1]^3.  So N3a uses the unjittered box for the
# free-surface statement, and N3b makes the stronger statement -- exactness on
# a badly distorted mesh -- through the periodic constraint instead.
hom = {0: ice}

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
        # For F = diag(lam,1,1) the Hencky strain is diag(ln lam,0,0), so the
        # Cauchy stress is T_xx = (K + 4G/3) ln(lam).  By Nanson the x = 1 face
        # keeps unit CURRENT area (J = lam, F^-T = diag(1/lam,1,1)), so the
        # resultant its nodes carry is T_xx itself.
        Txx = (ice[0] + 4.0 * ice[1] / 3.0) * np.log(lam)
        fx = f.reshape(-1, 3)[face0, 0].sum()
        chk('c=%d lam=%.3f  face resultant vs analytic Hencky' % (c, lam),
            abs(fx / Txx - 1.0), 0.0, 1e-9)

print('  N3b  jittered periodic cell: affine field is an exact solution')
mat1 = np.zeros_like(mat)
for c in (0, 1, 2, 3):
    cn, fln, it, rn = P.run_nl('fbar_%d' % c, 6, (0.375, 0.625), {0: ice, 1: ice},
                               eps=1e-3, jitter=0.3, geomname='sphere')
    chk('c=%d  fluctuation of the homogeneous cell' % c, fln, 0.0, 1e-7)
    exact = ice[0] + 4.0 * ice[1] / 3.0
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
            chk('c=%d  |C1111_nl/C1111_lin - 1| at eps=1e-5' % c, rel, 0.0, 5e-4)

print('\nN5  frame indifference: a rigid rotation carries no force')
# Under u = (Q - I) x every element sees F = Q, so J = 1, F~ = Q, J~ = 1,
# Jbar = 1 and Fbar = Q.  Then B = Q Q^T = I and H = 0, so T = 0 exactly.
# This is the check that a finite-strain formulation is objective, and it is
# the one that would catch a missing rotation anywhere in eqs. (1)-(11).
th = 0.7
Q = np.array([[np.cos(th), -np.sin(th), 0.0],
              [np.sin(th), np.cos(th), 0.0],
              [0.0, 0.0, 1.0]])
ax = np.array([1.0, 2.0, 3.0]); ax /= np.linalg.norm(ax)
Kx = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
Q = np.eye(3) + np.sin(th) * Kx + (1 - np.cos(th)) * (Kx @ Kx)
urot = ((Q - np.eye(3)) @ nodes.T).T.ravel()
for c in (0, 1, 2, 3):
    nl = P.FbarNL(nodes, tets, mat, g, vol, props, c)
    fr = nl.force(urot)
    # scale by the force a comparable STRETCH produces, else 'small' is meaningless
    ustr = np.zeros(3 * len(nodes)); ustr[0::3] = 0.1 * nodes[:, 0]
    ref = np.abs(nl.force(ustr)).max()
    chk('c=%d  max|f(rigid rotation)| / max|f(stretch)|' % c,
        np.abs(fr).max() / ref, 0.0, 1e-11)

print('\nN6  the linear operator IS df/du at the reference configuration')
# This is the check that matters for a port: the small-strain stiffness that
# would be implemented in CalculiX must be exactly the consistent tangent of
# the paper's element at u = 0, column for column -- not merely something that
# gives similar answers.  Central differences on the finite-strain force.
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
        e = np.zeros(nd); e[j] = hstep
        Jnum[:, j] = (nlj.force(e) - nlj.force(-e)) / (2.0 * hstep)
    den = np.abs(Klin).max()
    chk('c=%d  max|df/du - K_lin| / max|K_lin|' % c,
        np.abs(Jnum - Klin).max() / den, 0.0, 2e-6)
    # and the structure the port has to reproduce: dev block symmetric,
    # volumetric block not (section 4 of the formulation)
    asym = np.abs(Klin - Klin.T).max() / den
    print('       c=%d  asymmetry max|K - K^T|/max|K| = %.3e %s'
          % (c, asym, '(non-symmetric, as eq. 17 requires)' if c else
             '(symmetric, as S = E at c=0 requires)'))
    if c == 0 and asym > 1e-12:
        fails.append('c=0 must be symmetric')
    if c > 0 and asym < 1e-6:
        fails.append('c=%d must be non-symmetric' % c)

print('\nN7  rigid-body content of the assembled operator')
# A free (unconstrained) assembly must have exactly 6 zero modes and no more:
# 3 translations + 3 rotations, and nothing spurious.
for c in (0, 1, 2, 3):
    Klin = P.assemble('fbar_%d' % c, n4, t4, m4, g4, v4, fc4, pt4, props,
                      0.0, False, None).toarray()
    Ksym = 0.5 * (Klin + Klin.T)
    ev = np.linalg.eigvalsh(Ksym)
    scale = ev.max()
    nz = int((np.abs(ev) < 1e-9 * scale).sum())
    print('       c=%d  zero modes = %d (expect 6)   lambda_7/lambda_max = %.3e'
          % (c, nz, ev[6] / scale))
    if nz != 6:
        fails.append('c=%d has %d zero modes' % (c, nz))

print('\n%s' % ('all checks passed' if not fails
                else 'FAILED: ' + ', '.join(fails)))
sys.exit(1 if fails else 0)
