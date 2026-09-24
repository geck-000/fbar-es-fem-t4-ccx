"""Element-level stability test: spurious-mode census on a homogeneous block.

The patch test says an element is CONSISTENT.  The cantilever sweep says it
does not LOCK.  Neither can see an inf-sup failure -- unstable elements pass
both.  This is the test that sees it, and it needs no RVE.

Construction.  A homogeneous cube of near-incompressible material, every
boundary node clamped, so there are no rigid modes and the only thing left is
the element's own kinematics.  Take the lowest eigenvalues of K and, for each
mode, measure how much of its strain is VOLUMETRIC:

    r = sum_e V_e (div u|_e)^2 / sum_e V_e (eps_e : eps_e)

A genuine soft mode at K/G = 5000 is a shear mode: it costs G-level energy and
has r ~ 0.  A SPURIOUS mode is one that dilates the elements -- r of order one
-- and still costs only G-level energy, because the scheme's volumetric
operator cannot see that dilatation.  That combination is the signature, and
it is absolute: it needs no reference element to compare against.

What each pathology looks like:

  locking     lambda_1 grows like K as nu -> 1/2; the element has too few soft
              modes and every remaining one is expensive.
  stable      lambda_1 stays at the G scale and every soft mode has r ~ 0.
  unstable    lambda_1 stays at the G scale but soft modes carry r ~ 1.
"""
import sys
import warnings

import numpy as np
import scipy.sparse.linalg as spl

warnings.filterwarnings('ignore')
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import smoothing_proto as S                                  # noqa: E402


def census(scheme, n, K, G, jitter=0.3, stab=0.0, bubble=False, nev=30,
           two_phase=False, ice=None, radius=0.30):
    """two_phase=True embeds a brine sphere in ice.

    A HOMOGENEOUS cube cannot show the instability at all: every node patch
    reaches the clamped boundary, so none of them is unanchored and the null
    mode has nowhere to live.  Measured -- ns_vol reads n_bad = 0 there while
    the same scheme reads a 19x fluctuation on the real cells.  The soft phase
    has to be BOUNDED BY THE STIFF ONE for the test to mean anything."""
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
    out = []
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


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    G = 4.4e5
    ice = (9.37e9 / (3 * (1 - 2 * 0.33)), 9.37e9 / (2 * 1.33))
    arms = [('c3d4', 0.0, False, 'c3d4'),
            ('ns_vol', 0.0, False, 'ns_vol (NS-FEM) s=0'),
            ('ns_vol', 0.07, False, 'ns_vol s=0.07'),
            ('fs_ns', 0.0, False, 'fs_ns (FS/NS-FEM)')]
    print('TWO-PHASE clamped cube: brine sphere r=0.30 in ice, n=%d, '
          'jitter 0.3' % n)
    print('r is measured over the SOFT PHASE only -- the only place the null '
          'mode can live.\n')
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


if __name__ == '__main__':
    main()
