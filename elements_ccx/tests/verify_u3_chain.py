"""Mirror u3vol.f's element-weight walk in Python and check it reproduces the
verified operator S = E A^c from smoothing_proto.fbar_es_operators.

The Fortran builds the chain as a ROW walk over element weight lists, because
that is what an element can do with only a node->element map.  The prototype
builds it as sparse matrix products.  They must agree exactly; if they do not,
the CalculiX element is a different method from the one that was validated.

Also checks the two properties the walk must have for the patch test to pass:
unit row sums (a constant J survives the chain) and the per-material
restriction of V_n in eq. (6).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as P

fails = []


def walk(nodes, tets, mat, g, vol, imat, edge, ncyc):
    """Exactly u3vol.f: eq. (8) ring, then ncyc cycles of A = P Q."""
    at = {}
    for e in range(len(tets)):
        if mat[e] != imat:
            continue
        for a in tets[e]:
            at.setdefault(int(a), []).append(e)
    na, nb = edge
    cur = {}
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
        nxt = {}
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
            # compare the chain row against S[h, :]
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
                fails.append('c=%d jit=%.1f %s' % (ncyc, jitter, nm))

print('\n%s' % ('all checks passed' if not fails
                else 'FAILED: ' + ', '.join(fails)))
sys.exit(1 if fails else 0)
