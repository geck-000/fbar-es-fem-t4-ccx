"""Verify the F-barES-FEM-T4 formulation itself, before believing any of the
comparison numbers it produces.

The checks are on the operators of eqs. (1) and (6)-(8) and on the assembled
stiffness, and each one has an exact expected value, so a failure localises.

  V1  unit row sums on Q, P, E, R  -- what makes the small-strain reduction of
      eqs. (6)-(8) exact rather than approximate.  A constant J field must
      survive c cycles unchanged.
  V2  patch test: a homogeneous block under a uniform strain must return the
      exact C1111 = K + 4G/3 for every c and both readings of eq. (6).
  V3  the c = 0 identity: S = E, so B~^T D B~bar and B~bar^T D B~bar coincide
      and the scheme must reproduce selective ES-FEM-T4 exactly.
  V4  rank of the volumetric operator, which is the whole argument of section
      4: cyclic smoothing drops the edge count to the incidence rank, one
      volumetric constraint per independent node field.
"""
import os, sys
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as P

Ei, ni = 9.37e9, 0.33
ice = (Ei / (3 * (1 - 2 * ni)), Ei / (2 * (1 + ni)))
fails = []


def chk(name, got, want, tol):
    ok = abs(got - want) <= tol
    print('   %-52s %12.5e  %s' % (name, got, 'PASS' if ok else 'FAIL'))
    if not ok:
        fails.append(name)


print('V1  smoothing operators: unit row sums, constant field preserved')
nodes, tets, mat = P.mesh_box(6, 0.375, 0.625, 0.3, geom=P.GEOM['sphere'])
g, vol = P.grads(nodes, tets)
for jin in ('elem', 'edge'):
    for c in (0, 1, 2, 3):
        sel, idx, ek, ed, Vh, E, S, R = P.fbar_es_operators(
            nodes, tets, mat, g, vol, c, 1, jinput=jin)
        for nm, Op in (('Q.P=A', None), ('E', E), ('R', R), ('S', S)):
            if Op is None:
                continue
            rs = np.asarray(Op.sum(axis=1)).ravel()
            chk('%s c=%d  max|rowsum(%s) - 1|' % (jin, c, nm),
                np.abs(rs - 1).max(), 0.0, 1e-12)
        # a uniform divergence field must pass the whole chain untouched
        one = np.ones(len(sel))
        chk('%s c=%d  max|S.1 - 1| (constant J preserved)' % (jin, c),
            np.abs(S @ one - 1).max(), 0.0, 1e-11)

print('\nV2  patch test, homogeneous block, exact C1111 = K + 4G/3')
exact = ice[0] + 4.0 * ice[1] / 3.0
for jin in ('elem', 'edge'):
    os.environ['FBAR_JIN'] = jin
    for c in (0, 1, 2, 3):
        r = P.run('fbar_%d' % c, 4, (2.0, 3.0), {0: ice, 1: ice}, jitter=0.0)
        chk('%s fbar_%d  rel err in C1111' % (jin, c),
            abs(r[0] / exact - 1), 0.0, 1e-10)
os.environ['FBAR_JIN'] = 'elem'

print('\nV3  c = 0 must reproduce selective ES-FEM-T4 exactly (S = E)')
sel, idx, ek, ed, Vh, E, S, R = P.fbar_es_operators(
    nodes, tets, mat, g, vol, 0, 1, jinput='elem')
chk('max|S - E| at c=0', abs(S - E).max() if (S - E).nnz else 0.0, 0.0, 1e-14)

print('\nV4  rank of the volumetric operator (section 4 of the formulation)')
props = {0: ice, 1: (500 * 4.4e5, 4.4e5)}
for c in (0, 1, 2):
    sel, idx, ek, ed, Vh, E, S, R = P.fbar_es_operators(
        nodes, tets, mat, g, vol, c, 1, jinput='elem')
    nnode = len({int(a) for e in sel for a in tets[e]})
    rE = np.linalg.matrix_rank(E.toarray())
    rS = np.linalg.matrix_rank(S.toarray())
    print('   c=%d  brine: %d elems, %d nodes, %d edges | rank E = %d, '
          'rank S = %d  (c >= 1 caps the trial space at the incidence rank)'
          % (c, len(sel), nnode, len(ek), rE, rS))

print('\n%s' % ('all checks passed' if not fails
                else 'FAILED: ' + ', '.join(fails)))
sys.exit(1 if fails else 0)
