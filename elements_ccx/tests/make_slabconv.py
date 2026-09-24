"""Emit a mesh-convergence cell with a thin compliant layer resolved.

    make_slabconv.py OUTSTEM n [state] [kg] [jitter] [bridge] [load] [confine]

WHY THIS EXISTS.  A thin compliant layer cannot be resolved at practical mesh
sizes: with the layer ~0.0125 thick, the meshes actually run give 0.35 to 1.04
elements through it, so there is no asymptotic regime and no convergence rate
to measure.  The layer is made THICK instead of the mesh fine.  This is a
question about meshes, not about a particular cell, so the geometry is free to
be a different one.

GEOMETRY.  A soft slab spanning x in [lo, hi], pierced by a square lattice of
stiff bridges running in x.  Deliberately NOT a plain slab: under load a plain
layered cell is one-dimensional, the strain is uniform inside each phase, the
nodal average of the divergence is exact, and every scheme returns
bit-identical answers.  The bridges make the soft layer three-dimensionally
confined.

EVERY PHASE BOUNDARY IS A MULTIPLE OF 0.1, AND THAT IS THE WHOLE POINT.

smoothing_proto's own `bridged` geometry puts the bridge edges at
1/6 +- 0.06, which lands on no mesh plane: the geometry then changes with the
mesh and a convergence study measures the geometry drifting rather than the
discretisation converging.  Here the slab faces sit at 0.4 and 0.6 and the
bridges span [0.2, 0.4] and [0.6, 0.8] in both y and z, so with n a multiple
of 10 every phase boundary is a mesh plane, the phase assignment is exact at
every n, and the only thing that changes between meshes is h.  n is REQUIRED
to be a multiple of 10.

    elements through the slab = 0.2 n

n = 10, 20, 30, 40, 50, 60 gives 2 to 12 elements through the slab at
6e3 .. 1.3e6 tets.

LOADING, AND WHY IT IS ALONG THE SLAB AND NOT ACROSS IT.

Confined compression ACROSS the slab (drive the slab normal) looks like the
obvious locking test and is not one: the soft layer then deforms in uniaxial
strain -- eps_xx only -- a uniform state with div u =/= 0, so there is no
isochoric constraint to violate and every scheme gets it right.

So the load runs ALONG the slab (drive y), with the cell confined in x and z:

    y = 0   u_y = 0        y = 1   u_y = eps
    x = 0, x = 1           u_x = 0
    z = 0, z = 1           u_z = 0

The soft layer is compressed in its own plane while the stiff phase on either
side blocks expansion through the slab normal and the z rollers block the
rest.  It has to deform at nearly constant volume, which is the constraint
that locks a displacement tet.  Well posed, no rigid modes, and the lateral
planes are symmetry planes of the bridge lattice, so the rollers are exact.

    C1111_eff = sum RF_y(y = 1) / eps      R = C1111(und) / C1111(drn)

TWO DECKS PER CALL.  <OUTSTEM>_abq.inp uses C3D4H in the inclusion;
<OUTSTEM>_ccx.inp uses C3D4 throughout and asks for PARDISO, and is what
fbares.py converts to F-barES-FEM-T4.  Same nodes, same elements, same
boundary conditions, so the two codes look at the same model.
"""
import os
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import smoothing_proto as S                                  # noqa: E402

G_SOFT = 440029.33528897085           # representative soft-phase shear modulus
E_STIFF, NU_STIFF = 9.37e9, 0.33
KG_DRAINED = 5.0                      # K = 2.2e6 against G = 4.4e5
EPS = 1.0e-3                          # applied axial strain


def iso(kg, g=G_SOFT):
    """E, nu from a bulk/shear ratio."""
    k = kg * g
    return 9.0 * k * g / (3.0 * k + g), (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g))


def main():
    stem, n = sys.argv[1], int(sys.argv[2])
    state = sys.argv[3] if len(sys.argv) > 3 else 'und'
    kg = float(sys.argv[4]) if len(sys.argv) > 4 else 500.0
    jit = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
    bridge = sys.argv[6] if len(sys.argv) > 6 else 'one'
    load = sys.argv[7] if len(sys.argv) > 7 else 'x'
    confine = sys.argv[8] if len(sys.argv) > 8 else 'sym'
    lo, hi = 0.4, 0.6

    if n % 10:
        raise SystemExit('make_slabconv: n must be a multiple of 10 so that '
                         'every phase boundary (x=%.1f, %.1f and the bridge '
                         'edges at 0.2/0.4/0.6/0.8) lands on a mesh plane.  A '
                         'staircased boundary makes the GEOMETRY change with '
                         'the mesh, which is a mesh-dependent geometry error '
                         'sitting inside the convergence study it is meant to '
                         'measure -- it read R = 1.1054, 1.5388, 1.1308 at '
                         'n = 10, 15, 20.' % (lo, hi))

    # HOW MUCH STIFF PHASE PIERCES THE SLAB decides whether the cell locks.
    # With four 0.2x0.2 bridges (16% of the slab area) the drained cell is
    # still carried by the stiff phase, R is only 1.18, and C3D4 and F-barES
    # agree to 0.05% -- a cell where the element does not matter cannot
    # discriminate between elements.  Thinning the bridges forces the load
    # through the soft phase, which is what confines it and what makes C3D4
    # lock.
    #
    # Every edge is a multiple of 0.1 either way, so n % 10 == 0 resolves the
    # geometry exactly at any n.
    PATTERNS = {
        'four': lambda t: (0.2 <= t < 0.4) or (0.6 <= t < 0.8),   # 16% stiff
        'one':  lambda t: 0.4 <= t < 0.6,                          #  4% stiff
        'none': lambda t: False,                                   #  0%, 1-D
    }
    inb = PATTERNS[bridge]

    def bridged(x, y, z, blo, bhi):
        if not (blo <= x < bhi):
            return 0
        return 0 if (inb(y) and inb(z)) else 1

    nodes, tets, mat = S.mesh_box(n, lo, hi, 0.0, geom=bridged)

    if jit:
        # INTERIOR NODES ONLY.  mesh_box displaces boundary nodes too (in
        # matched periodic pairs), but these decks select faces by coordinate,
        # so a displaced boundary node is no longer ON the face and never gets
        # constrained -- which silently unrestrains the cell.  Distorting the
        # interior is what the test needs anyway.
        qmin = float(os.environ.get('FBAR_BLOCK_QMIN', 0.15))
        rng = np.random.default_rng(11)
        tol0 = 1e-12
        inter = ~((nodes <= tol0) | (nodes >= 1.0 - tol0)).any(axis=1)
        d = (rng.random(nodes.shape) - 0.5) * 2.0 / n
        base, amp = nodes.copy(), jit
        for _ in range(40):
            nodes = base.copy()
            nodes[inter] += amp * d[inter]
            p4 = nodes[tets]
            vv = np.linalg.det(p4[:, 1:, :] - p4[:, :1, :]) / 6.0
            if vv.min() > 0 and vv.min() >= qmin * vv.mean():
                break
            amp *= 0.8
        else:
            raise SystemExit('make_slabconv: cannot jitter n=%d without '
                             'tangling' % n)
        sys.stderr.write('  jitter %.3f -> %.3f (min/mean tet volume %.3f)\n'
                         % (jit, amp, vv.min() / vv.mean()))
    S.grads(nodes, tets)                       # fixes any inverted tets

    e_soft, nu_soft = iso(kg if state == 'und' else KG_DRAINED)
    soft = [e for e in range(len(tets)) if mat[e] == 1]
    hard = [e for e in range(len(tets)) if mat[e] == 0]
    tol = 1e-9

    def nset(name, sel):
        ids = [str(i + 1) for i in np.flatnonzero(sel)]
        out = ['*NSET,NSET=' + name]
        out += [', '.join(ids[c:c + 12]) for c in range(0, len(ids), 12)]
        return out

    common = ['*NODE']
    common += ['%d, %.12e, %.12e, %.12e' % (i + 1, p[0], p[1], p[2])
               for i, p in enumerate(nodes)]

    def elblock(kind, ids, elset):
        out = ['*ELEMENT,TYPE=%s,ELSET=%s' % (kind, elset)]
        out += ['%d, %s' % (e + 1, ', '.join(str(x + 1) for x in tets[e]))
                for e in ids]
        return out

    tail = []
    # X0/X1 keep their names -- the driven pair, whatever axis it is on -- so
    # the extractor does not have to know the load direction.
    ax = {'x': 0, 'y': 1, 'z': 2}[load]
    o1, o2 = [c for c in (0, 1, 2) if c != ax]
    tail += nset('X0', (nodes[:, ax] <= tol))
    tail += nset('X1', (nodes[:, ax] >= 1.0 - tol))
    # CONFINEMENT IS WHAT DECIDES WHETHER THE CELL LOCKS AT ALL.
    #
    # 'full' rollers both lateral pairs, which puts the whole cell in uniaxial
    # STRAIN.  The soft phase then has to change volume and incompressibility
    # never
    # binds -- measured, C3D4 and F-barES agreed to 0.01-0.2% across every
    # bridge pattern, both load directions and both K/G, while R ranged over
    # 1.03 to 7.16.  A cell that does not lock cannot tell two elements apart.
    #
    # 'sym' rollers only the low faces, leaving the high ones free: uniaxial
    # STRESS.  Now the cell contracts laterally, the stiff phase at nu = 0.33
    # and the soft phase at nu -> 0.5 disagree about by how much, and the soft
    # phase is forced to deform at nearly constant volume by the stiff phase
    # around it.  That is the
    # constraint that locks a displacement tet.
    if confine == 'full':
        ysel = (nodes[:, o1] <= tol) | (nodes[:, o1] >= 1.0 - tol)
        zsel = (nodes[:, o2] <= tol) | (nodes[:, o2] >= 1.0 - tol)
    else:
        ysel = nodes[:, o1] <= tol
        zsel = nodes[:, o2] <= tol
    tail += nset('YFACE', ysel)
    tail += nset('ZFACE', zsel)
    tail += ['*SOLID SECTION,ELSET=Sphere_Only,MATERIAL=Mat_Inclusion',
             '*MATERIAL,NAME=Mat_Inclusion', '*ELASTIC',
             '%.12e, %.12e' % (e_soft, nu_soft)]
    if hard:
        tail += ['*SOLID SECTION,ELSET=Matrix_Only,MATERIAL=Mat_Matrix',
                 '*MATERIAL,NAME=Mat_Matrix', '*ELASTIC',
                 '%.12e, %.12e' % (E_STIFF, NU_STIFF)]
    # The driven X1 magnitude is non-zero, so *BOUNDARY must sit AFTER *STEP:
    # Abaqus rejects any prescribed magnitude other than zero in the model
    # definition ("PRESCRIBED *BOUNDARY MAGNITUDES MUST BE ZERO IN THE MODEL
    # DEFINITION").  CalculiX tolerates either placement, so emitting the block
    # in the step keeps both decks correct.
    bc = ['*BOUNDARY',
          'X0, %d, %d, 0.0' % (ax + 1, ax + 1),
          'X1, %d, %d, %.12e' % (ax + 1, ax + 1, EPS),
          'YFACE, %d, %d, 0.0' % (o1 + 1, o1 + 1),
          'ZFACE, %d, %d, 0.0' % (o2 + 1, o2 + 1)]

    for tag, kinc, step in (('abq', 'C3D4H', '*STATIC'),
                            ('ccx', 'C3D4', '*STATIC, SOLVER=PARDISO')):
        L = list(common)
        L += elblock(kinc, soft, 'Sphere_Only')
        if hard:
            L += elblock('C3D4', hard, 'Matrix_Only')
        L += tail
        L += ['*STEP', step] + bc
        L += ['*NODE PRINT,NSET=X1', 'U,RF',
              '*NODE PRINT,NSET=X0', 'U,RF',
              '*END STEP']
        open('%s_%s.inp' % (stem, tag), 'w').write('\n'.join(L) + '\n')

    k = kg if state == 'und' else KG_DRAINED
    print('n=%-3d  %s  K/G=%-6g  bridge=%s load=%s/%s  nodes=%d tets=%d (soft %d, %.1f%%)  '
          'elements through slab = %.0f'
          % (n, state, k, bridge, load, confine, len(nodes), len(tets), len(soft),
             100.0 * len(soft) / len(tets), (hi - lo) * n))
    print('  E_soft=%.6e nu_soft=%.9f   eps=%.1e   ->  %s_{abq,ccx}.inp'
          % (e_soft, nu_soft, EPS, stem))


if __name__ == '__main__':
    main()
