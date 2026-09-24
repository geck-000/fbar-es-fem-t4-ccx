"""Numerical comparison of strain-smoothing schemes on linear tets.

Written to settle, before any Fortran, whether a FACE-based smoothing fixes
the instability that node-based smoothing (NS-FEM) has on the inclusion.

The cell is the failure in miniature: a soft slab normal to x, embedded in a
stiff matrix, driven by a prescribed affine displacement u = eps.X on the
whole boundary.  Two numbers come out.

  C1111   2U/(eps^2 V), the confined stiffness.  This is the locking-sensitive
          modulus -- it contains K directly.  Too high = locking.
  fluc    max |u - eps.X| over the interior, divided by the applied
          displacement.  This is the mode detector used on the real cells:
          an affine boundary condition admits NO fluctuation in the exact
          solution of a homogeneous body, and in a two-phase body only a
          small one.  NS-FEM can read 19x on such cells.

Schemes:
  c3d4     standard displacement tet
  ns_vol   element deviatoric + NODE-smoothed volumetric  (NS-FEM)
  ns_full  node-smoothed full strain                      (textbook NS-FEM)
  fs_full  face-smoothed full strain                      (textbook FS-FEM)
  fs_ns    FACE-smoothed deviatoric + node-smoothed volumetric
           (the selective FS/NS-FEM recommended for near-incompressible)
  fs_vol   face-smoothed volumetric + element deviatoric  (completeness)

Smoothing never crosses a material interface: a face whose two tets have
different materials is treated as a boundary face, and node patches are keyed
by (node, material).
"""
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

M = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
T = np.diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])


def dmat(K, G, part='full'):
    """Voigt [xx,yy,zz,gxy,gyz,gzx]; energy = 1/2 eps^T D eps."""
    vol = K * np.outer(M, M)
    dev = 2.0 * G * (T - np.outer(M, M) / 3.0)
    if part == 'vol':
        return vol
    if part == 'dev':
        return dev
    return vol + dev


_JIT_SHRUNK = {}


def _tets_of(n):
    """Connectivity only, for the jitter quality check."""
    def nid(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k
    split = ((0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
             (0, 3, 2, 7), (0, 6, 4, 7), (0, 2, 6, 7))
    t = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = [nid(i + (v & 1), j + ((v >> 1) & 1), k + ((v >> 2) & 1))
                     for v in range(8)]
                for sp_ in split:
                    t.append([c[sp_[0]], c[sp_[1]], c[sp_[2]], c[sp_[3]]])
    return np.array(t)


def mesh_box(n, slab_lo, slab_hi, jitter=0.0, seed=7, geom=None):
    """n^3 hexes, 6 tets each (Freudenthal).  Returns nodes, tets, mat.

    `jitter` displaces interior nodes by up to jitter*h.  A perfectly
    structured mesh is the best case for a nodal-smoothing scheme; the real
    cells are unstructured Delaunay, and the spurious mode is known to be
    mesh-dependent, so the comparison has to be made on a disturbed mesh.
    Nodes on the periodic boundary are moved in matched pairs so that
    periodicity survives."""
    if geom is None:
        geom = GEOM['slab']
    g = np.linspace(0.0, 1.0, n + 1)
    X, Y, Z = np.meshgrid(g, g, g, indexing='ij')
    nodes = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    if jitter:
        rng = np.random.default_rng(seed)
        h = 1.0 / n
        # draw one displacement per PERIODIC IMAGE, so paired faces move alike
        key = np.round(np.mod(nodes * n, n)).astype(int)      # 0..n-1, wrapped
        uniq, inv = np.unique(key, axis=0, return_inverse=True)
        d = (rng.random((len(uniq), 3)) - 0.5) * 2.0 * jitter * h
        # QUALITY CONTROL, and it is not optional.  At jitter 0.3 this mesh
        # TANGLES: 1 inverted tet at n=6, 14 at n=16.  grads() silently flips
        # their node order, which makes the volumes positive again and hides
        # it, but the tets still overlap their neighbours, the mesh is no
        # longer a partition of the cell, and every scheme -- C3D4 included --
        # then fails its own patch test (measured: fluc 1.4e-3 and a residual
        # of 0.21 against the exact affine field, at interior nodes).
        # Shrink the amplitude until no tet is inverted and the worst is still
        # a reasonable fraction of the mean.
        base = nodes.copy()
        amp, qmin = 1.0, float(os.environ.get('FBAR_MESH_QMIN', '0.15'))
        for _ in range(40):
            nodes = base + amp * d[inv]
            p4 = nodes[_tets_of(n)]
            vv = np.linalg.det(p4[:, 1:, :] - p4[:, :1, :]) / 6.0
            if vv.min() > 0 and vv.min() >= qmin * vv.mean():
                break
            amp *= 0.8
        else:
            raise RuntimeError('mesh_box: cannot jitter n=%d without tangling'
                               % n)
        if amp < 1.0:
            _JIT_SHRUNK[(n, jitter, seed)] = amp

    def nid(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k

    # the 6-tet split of a cube, all with positive volume for this corner order
    split = ((0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
             (0, 3, 2, 7), (0, 6, 4, 7), (0, 2, 6, 7))
    tets, mat = [], []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = [nid(i + (v & 1), j + ((v >> 1) & 1), k + ((v >> 2) & 1))
                     for v in range(8)]
                xc, yc, zc = (i + 0.5) / n, (j + 0.5) / n, (k + 0.5) / n
                m = geom(xc, yc, zc, slab_lo, slab_hi)      # 1 = soft phase
                for s in split:
                    tets.append([c[s[0]], c[s[1]], c[s[2]], c[s[3]]])
                    mat.append(m)
    return nodes, np.array(tets), np.array(mat)


# --------------------------------------------------------------------------
# GEOMETRY.  A soft slab spanning the cell is useless as a test: under a
# periodic eps_xx it is a one-dimensional layered problem, the strain is
# uniform inside each phase, so the nodal average of the divergence is EXACT
# and every scheme returns bit-identical answers (measured: 5.9656e+09 for
# c3d4, ns_vol and fs_ns alike).  The spurious mode needs soft-phase nodes
# that touch no stiff element at all -- a three-dimensional pocket.
def _slab(x, y, z, lo, hi):
    return 1 if lo <= x < hi else 0


def _sphere(x, y, z, lo, hi):
    r = hi                                   # hi carries the radius
    return 1 if (x - .5) ** 2 + (y - .5) ** 2 + (z - .5) ** 2 < r * r else 0


def _bridged(x, y, z, lo, hi):
    """Soft slab pierced by a square lattice of stiff bridges, the geometry
    where the mode was measured."""
    if not (lo <= x < hi):
        return 0
    b = 0.18                                 # bridge half-width
    if (abs(((y * 3) % 1.0) - .5) < b) and (abs(((z * 3) % 1.0) - .5) < b):
        return 0                             # stiff bridge
    return 1


GEOM = {'slab': _slab, 'sphere': _sphere, 'bridged': _bridged}


def grads(nodes, tets):
    """Shape-function gradients and volumes for every tet."""
    p = nodes[tets]                                  # (ne,4,3)
    e = p[:, 1:, :] - p[:, :1, :]                    # (ne,3,3)
    det = np.linalg.det(e)
    vol = det / 6.0
    if np.any(vol <= 0):
        # flip the last two nodes of the inverted ones
        bad = vol <= 0
        tets[bad] = tets[bad][:, [0, 1, 3, 2]]
        p = nodes[tets]
        e = p[:, 1:, :] - p[:, :1, :]
        vol = np.linalg.det(e) / 6.0
    inv = np.linalg.inv(e)                           # (ne,3,3)
    g = np.zeros((len(tets), 4, 3))
    g[:, 1:, :] = np.transpose(inv, (0, 2, 1))
    g[:, 0, :] = -g[:, 1:, :].sum(axis=1)
    return g, vol


def bmat(g):
    """6x12 strain-displacement for one tet from its 4x3 gradients."""
    B = np.zeros((6, 12))
    for a in range(4):
        gx, gy, gz = g[a]
        c = 3 * a
        B[0, c + 0] = gx
        B[1, c + 1] = gy
        B[2, c + 2] = gz
        B[3, c + 0] = gy; B[3, c + 1] = gx
        B[4, c + 1] = gz; B[4, c + 2] = gy
        B[5, c + 0] = gz; B[5, c + 2] = gx
    return B


# --------------------------------------------------------------------------
# BUBBLE ENRICHMENT (bFS-FEM, Techscience CMES v127n2).
#
# u = u_lin + phi_b * u_b per element, phi_b = 256 L1 L2 L3 L4, 3 internal dof.
#
# Over a WHOLE element int grad(phi_b) dV = 0 (the bubble vanishes on every
# face), so element-averaged smoothing cannot see the bubble at all.  It is
# visible only to SUB-CELL averages, which is exactly why the S-FEM smoothing
# domains are built by splitting each tet at its centroid.  For the
# face-based domains the sub-cell is (3 face nodes + centroid), volume V/4.
#
# The bubble dofs are NOT condensed here.  With face smoothing an element's
# bubble couples to its four face-neighbours' bubbles, so the bubble block is
# not block-diagonal and element-local condensation is invalid.  Carrying the
# 3 dofs explicitly is exact and this is a prototype.
QP = np.array([[0.5854102, 0.1381966, 0.1381966, 0.1381966],
               [0.1381966, 0.5854102, 0.1381966, 0.1381966],
               [0.1381966, 0.1381966, 0.5854102, 0.1381966],
               [0.1381966, 0.1381966, 0.1381966, 0.5854102]])


def bub_pattern(gp):
    """6x3 strain matrix of the field phi*e_c given grad(phi) = gp."""
    B = np.zeros((6, 3))
    gx, gy, gz = gp
    B[0, 0] = gx
    B[1, 1] = gy
    B[2, 2] = gz
    B[3, 0] = gy; B[3, 1] = gx
    B[4, 1] = gz; B[4, 2] = gy
    B[5, 0] = gz; B[5, 2] = gx
    return B


def subcell_bubble_grads(nodes, tets, g):
    """Mean grad(phi_b) over each of a tet's 4 face sub-cells -> (ne,4,3)."""
    out = np.zeros((len(tets), 4, 3))
    for e, t in enumerate(tets):
        p = nodes[t]                       # (4,3) parent nodes
        cen = p.mean(axis=0)
        for drop in range(4):              # sub-cell opposite node `drop`
            verts = np.vstack([p[[i for i in range(4) if i != drop]], cen])
            acc = np.zeros(3)
            for w in QP:                   # 4-point degree-2 rule, weights 1/4
                x = w @ verts
                L = np.array([1.0 + g[e, k] @ (x - p[k]) for k in range(4)])
                gp = np.zeros(3)
                for k in range(4):
                    prod = np.prod([L[j] for j in range(4) if j != k])
                    gp += prod * g[e, k]
                acc += 0.25 * 256.0 * gp
            out[e, drop] = acc
    return out


def topology(tets, mat):
    """faces -> the tets sharing them; (node,material) -> its patch of tets."""
    faces = {}
    for e, t in enumerate(tets):
        for drop in range(4):
            f = tuple(sorted(int(t[i]) for i in range(4) if i != drop))
            faces.setdefault(f, []).append(e)
    patch = {}
    for e, t in enumerate(tets):
        for a in t:
            patch.setdefault((int(a), int(mat[e])), []).append(e)
    return faces, patch


class Assembler:
    def __init__(self, ndof):
        self.ndof = ndof
        self.r, self.c, self.v = [], [], []

    def add(self, dofs, Kloc):
        d = np.asarray(dofs)
        self.r.append(np.repeat(d, len(d)))
        self.c.append(np.tile(d, len(d)))
        self.v.append(Kloc.ravel())

    def tocsr(self):
        return sp.coo_matrix(
            (np.concatenate(self.v),
             (np.concatenate(self.r), np.concatenate(self.c))),
            shape=(self.ndof, self.ndof)).tocsr()


def dofs_of(nodelist):
    return np.array([3 * n + c for n in nodelist for c in range(3)])


def expand(B, tet, nodelist):
    """Re-index a tet's 6x12 B onto a wider node list."""
    idx = {n: i for i, n in enumerate(nodelist)}
    out = np.zeros((6, 3 * len(nodelist)))
    for a in range(4):
        j = 3 * idx[int(tet[a])]
        out[:, j:j + 3] = B[:, 3 * a:3 * a + 3]
    return out


def assemble(scheme, nodes, tets, mat, g, vol, faces, patch, props,
             stab=0.0, bubble=False, sbg=None):
    """props[m] = (K, G) for material m.  Bubble dofs follow the nodal ones."""
    nn = len(nodes)
    ndof = 3 * nn + (3 * len(tets) if bubble else 0)
    A = Assembler(ndof)
    B = [bmat(g[e]) for e in range(len(tets))]

    def bdofs(e):
        return np.array([3 * nn + 3 * e + c for c in range(3)])

    # ---- element-wise parts -------------------------------------------
    if scheme in ('c3d4', 'ns_vol', 'fs_vol'):
        for e, t in enumerate(tets):
            K, G = props[mat[e]]
            if scheme == 'c3d4':
                D = dmat(K, G, 'full')
            elif scheme == 'ns_vol':
                D = dmat(stab * K, G, 'full') if stab else dmat(K, G, 'dev')
            else:
                D = dmat(K, G, 'dev')
            A.add(dofs_of(t), vol[e] * B[e].T @ D @ B[e])

    # ---- face-smoothed parts ------------------------------------------
    if scheme in ('fs_full', 'fs_ns', 'fs_vol'):
        part = {'fs_full': 'full', 'fs_ns': 'dev', 'fs_vol': 'vol'}[scheme]
        for f, els in faces.items():
            # never smooth across a material interface
            groups = ([els] if len(els) == 2 and mat[els[0]] == mat[els[1]]
                      else [[e] for e in els])
            for grp in groups:
                nl = sorted(set(int(x) for e in grp for x in tets[e]))
                Vk = sum(vol[e] for e in grp) / 4.0
                cols = [sum((vol[e] / 4.0) * expand(B[e], tets[e], nl)
                            for e in grp) / Vk]
                dl = [dofs_of(nl)]
                if bubble:
                    for e in grp:
                        drop = [i for i in range(4)
                                if int(tets[e][i]) not in f][0]
                        cols.append((vol[e] / 4.0) * bub_pattern(sbg[e, drop])
                                    / Vk)
                        dl.append(bdofs(e))
                Bb = np.hstack(cols)
                K, G = props[mat[grp[0]]]
                A.add(np.concatenate(dl), Vk * Bb.T @ dmat(K, G, part) @ Bb)

    # ---- F-barES-FEM-T4 -------------------------------------------------
    if scheme.startswith('fbar'):
        cyc = int(scheme.split('_')[-1]) if scheme.split('_')[-1].isdigit() else 1
        for m in sorted(set(int(x) for x in mat)):
            sel, idx, ekeys, edges, Vh, Edev, Svol, _R = fbar_es_operators(
                nodes, tets, mat, g, vol, cyc, m)
            K, G = props[m]
            for h, k in enumerate(ekeys):
                els_h = edges[k]
                nl = sorted(set(int(x) for e in els_h for x in tets[e]))
                # deviatoric: edge-smoothed once -- a narrow stencil, so a
                # dense block per edge is fine
                Bd = sum((vol[e] / 6.0 / Vh[h]) * expand(B[e], tets[e], nl)
                         for e in els_h)
                A.add(dofs_of(nl), Vh[h] * Bd.T @ dmat(K, G, 'dev') @ Bd)

            # VOLUMETRIC, as a sparse triple product.
            #
            # E A^c spans a wide element neighbourhood -- roughly 2c+1 rings --
            # so building one dense outer product per edge over that support
            # costs ~(3*300)^2 per edge and ran the machine out of memory.
            # The operator is  K_vol = K * Gv^T diag(V_h) Gv  with
            # Gv = (E A^c) D_div, and both factors are sparse, so the whole
            # thing is one sparse product instead of 13k dense blocks.
            rd, cd, vd = [], [], []
            for j, e in enumerate(sel):
                row = M @ B[e]                       # 1 x 12, the divergence
                for a in range(4):
                    for cc in range(3):
                        w = row[3 * a + cc]
                        if w:
                            rd.append(j)
                            cd.append(3 * int(tets[e][a]) + cc)
                            vd.append(w)
            Ddiv = sp.coo_matrix((vd, (rd, cd)),
                                 shape=(len(sel), ndof)).tocsr()
            # F-bar is Petrov-Galerkin, and that is the whole point.  Eq. (17)
            # of the paper: the stress comes from the modified gradient Fbar
            # (eq. 11, volumetric part smoothed c times) but it is paired with
            # the stretching of the UNMODIFIED edge gradient Ftilde (eq. 1).
            # So the volumetric operator is Btilde^T D Bbar, not Bbar^T D Bbar.
            #
            # Symmetrising it changes the coercivity, not just the count.  Both
            # pairings have the same kernel, but the
            # Galerkin form Bbar^T D Bbar is a Gram matrix, positive semidefinite
            # for every K, while the symmetric part of the PG form turns
            # indefinite above a mesh-dependent K/G threshold: measured
            # kappa* = 7367, 5114, 3592 at n = 6, 8, 12 (verify_pairing.py), with
            # the MMS error then growing linearly in K/G.  Inside the definite
            # range the two agree closely: on the sphere cell at K/G = 5000 the
            # fluctuations are 1.07 (PG) and 0.53 (Galerkin), so the earlier
            # 17.53 figure is not reproduced by this cell and should not be
            # quoted (see verify_pairing.py).  Keeping the edge test space is
            # what the published method specifies; the Galerkin form is the
            # stable fallback.
            #
            # At c = 0, Svol is E and the two forms coincide, which is why
            # fbar_0 was the only row that looked sane.
            Gtrial = (Svol @ Ddiv).tocsr()
            Gtest = (Edev @ Ddiv).tocsr()
            W = sp.diags(K * Vh)
            # FBAR_SYM=1 assembles the Galerkin form Bbar^T D Bbar
            # instead.  That is NOT what eq. (17) specifies, but it is
            # symmetric, and ccx calls PARDISO with mtype = -2 -- so whether
            # the difference is measurable decides whether a CalculiX port has
            # to touch the solver path at all.  Measured, not assumed.
            if os.environ.get('FBAR_SYM') == '1':
                Gtest = Gtrial
            Kv = (Gtest.T @ W @ Gtrial).tocoo()
            A.r.append(Kv.row); A.c.append(Kv.col); A.v.append(Kv.data)

    # ---- ES-FEM deviatoric + node-patch volumetric ----------------------
    #
    # This is F-barES-FEM-T4(1) with eq. (8) -- the final elem->edge smoothing
    # of J -- DROPPED, so the volumetric energy sits on node patches instead of
    # edge domains.  It matters because that is the one variant CalculiX can
    # already deliver: the volumetric half is exactly a node-smoothed patch
    # stiffness, and only the deviatoric half is new.  The full
    # eq. (8) operator cannot be an element (173-node stencil at c=1, 494 at
    # c=2 against *USER ELEMENT's 255-node limit) and cannot be *EQUATION +
    # SPRING1 either (tried earlier: overlapping patches gave a
    # confined-compression reaction 1.55x the closed form).
    if scheme == 'esdev_nsvol':
        for m in sorted(set(int(x) for x in mat)):
            sel, idx, ekeys, edges, Vh, E, S, R = fbar_es_operators(
                nodes, tets, mat, g, vol, 0, m)
            K, G = props[m]
            for h, k in enumerate(ekeys):
                els_h = edges[k]
                nl = sorted(set(int(x) for e in els_h for x in tets[e]))
                Bd = sum((vol[e] / 6.0 / Vh[h]) * expand(B[e], tets[e], nl)
                         for e in els_h)
                A.add(dofs_of(nl), Vh[h] * Bd.T @ dmat(K, G, 'dev') @ Bd)
        for (_a, m), els in patch.items():
            nl = sorted(set(int(x) for e in els for x in tets[e]))
            Va = sum(vol[e] for e in els) / 4.0
            K, G = props[m]
            r = sum((vol[e] / 4.0) * (M @ expand(B[e], tets[e], nl))
                    for e in els) / Va
            A.add(dofs_of(nl), K * Va * np.outer(r, r))

    # ---- node-smoothed parts ------------------------------------------
    if scheme in ('ns_vol', 'ns_full', 'fs_ns'):
        part = 'full' if scheme == 'ns_full' else 'vol'
        for (_a, m), els in patch.items():
            nl = sorted(set(int(x) for e in els for x in tets[e]))
            Va = sum(vol[e] for e in els) / 4.0
            K, G = props[m]
            if part == 'full':
                Bb = sum((vol[e] / 4.0) * expand(B[e], tets[e], nl)
                         for e in els) / Va
                A.add(dofs_of(nl), Va * Bb.T @ dmat(K, G, 'full') @ Bb)
                continue
            # volumetric only: one scalar row, theta_a
            row = [sum((vol[e] / 4.0) * (M @ expand(B[e], tets[e], nl))
                       for e in els) / Va]
            dl = [dofs_of(nl)]
            if bubble:
                # int_e N_a div(phi_b u_b) dV = -(32 V_e/105) g_a . u_b
                for e in els:
                    la = [i for i in range(4) if int(tets[e][i]) == a][0]
                    row.append(-(32.0 * vol[e] / 105.0) * g[e, la] / Va)
                    dl.append(bdofs(e))
            r = np.concatenate(row)
            A.add(np.concatenate(dl),
                  K * (1.0 - stab) * Va * np.outer(r, r))

    return A.tocsr()


def run(scheme, n, slab, props, eps=1e-3, stab=0.0, bubble=False,
        jitter=0.0, geomname='slab'):
    nodes, tets, mat = mesh_box(n, slab[0], slab[1], jitter,
                                geom=GEOM[geomname])
    g, vol = grads(nodes, tets)
    faces, patch = topology(tets, mat)
    sbg = subcell_bubble_grads(nodes, tets, g) if bubble else None
    K = assemble(scheme, nodes, tets, mat, g, vol, faces, patch, props,
                 stab, bubble, sbg)

    # PERIODIC boundary condition, u(x + L e_i) = u(x) + eps.(L e_i).
    #
    # This is what the real cells use, and it matters: with an affine
    # Dirichlet condition on the whole boundary the soft slab is clamped on
    # every side and the bellows mode simply cannot form -- every scheme then
    # reads a fluctuation of 0.2-0.4 and the test discriminates nothing.
    nn = len(nodes)
    ndof = K.shape[0]
    tol = 1e-9
    # master of each node: its image in [0,1)^3
    base = np.mod(np.round(nodes * 1e9) / 1e9, 1.0)
    base[np.abs(base - 1.0) < tol] = 0.0
    key = np.round(base * 1e7).astype(np.int64)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    first = np.full(len(uniq), -1, dtype=int)
    for i in range(nn):
        if first[inv[i]] < 0:
            first[inv[i]] = i
    master = first[inv]
    shift = nodes - nodes[master]                    # 0 or +/-1 per component

    # u = T q + a ; q are the master nodal dofs plus every bubble dof
    nm = len(uniq)
    nq = 3 * nm + (ndof - 3 * nn)
    rows, cols, vals = [], [], []
    for i in range(nn):
        m = inv[i]
        for c in range(3):
            rows.append(3 * i + c); cols.append(3 * m + c); vals.append(1.0)
    for j in range(3 * nn, ndof):
        rows.append(j); cols.append(3 * nm + (j - 3 * nn)); vals.append(1.0)
    T = sp.coo_matrix((vals, (rows, cols)), shape=(ndof, nq)).tocsr()

    a = np.zeros(ndof)
    a[0:3 * nn:3] = eps * shift[:, 0]                # only eps_xx is applied

    Kr = (T.T @ K @ T).tocsr()
    fr = -(T.T @ (K @ a))
    # one master node pinned to kill the rigid translation
    keep = np.ones(nq, dtype=bool)
    keep[0:3] = False
    q = np.zeros(nq)
    q[keep] = spl.spsolve(Kr[keep][:, keep].tocsc(), fr[keep])
    u = T @ q + a

    # 0.5 a.Ku, not 0.5 u.Ku: the two agree exactly whenever K is
    # symmetric (T^T K u = 0 kills the q part), and only the first is
    # right for the non-symmetric F-bar operator, where `a` is the
    # macroscopic test field and K u the internal force.
    energy = 0.5 * a @ (K @ u)
    c1111 = 2.0 * energy / (eps ** 2)                # unit volume
    # u is the TOTAL displacement (the jump `a` makes it so), defined up to a
    # constant.  The fluctuation is what is left after the macroscopic affine
    # field and that constant are removed -- the same quantity the real cells
    # report as max|u| against the applied displacement.
    un = u[:3 * nn]
    aff = np.zeros_like(un)
    aff[0::3] = eps * nodes[:, 0]
    # Remove the rigid translation PER COMPONENT.  u is interleaved
    # [ux,uy,uz,...], so a single scalar mean mixes the three directions and
    # leaves part of the translation that pinning q[0:3] = 0 introduced --
    # which showed up as C3D4 apparently failing its own patch test on a
    # jittered mesh (fluc 2e-3 where the answer is 0).
    fluct = (un - aff).reshape(-1, 3)
    fluct = fluct - fluct.mean(axis=0)
    fl = np.abs(fluct).max() / eps

    # PRESSURE OSCILLATION, measured rather than assumed.
    #
    # p_e = K div(u)|_e is the element pressure the material actually feels,
    # whatever the scheme used to build the stiffness.  A stable element gives
    # neighbouring soft-phase elements nearly the same pressure, so the face
    # jump
    # is O(h) against the mean.  A checkerboard gives jumps of order the mean
    # or larger.  Reported as mean and max jump over soft-soft faces,
    # normalised by the mean |p| in the soft phase.
    soft = np.flatnonzero(mat == 1)
    theta = np.zeros(len(tets))
    for e in range(len(tets)):
        theta[e] = M @ (bmat(g[e]) @ un[dofs_of(tets[e])])
    Kb = props[1][0]

    # TWO pressure fields, and the distinction is not cosmetic.
    #
    # praw = K div(u) per element is what a C3D4 would report from the same
    # displacement field.  It is scheme-independent, which makes it a fair
    # cross-scheme number, but it is NOT the pressure a smoothed scheme
    # computes -- and for any scheme that smooths the volumetric strain it is
    # biased against it by construction: smoothing the stress while leaving u
    # free necessarily roughens the raw element divergence.  Reporting only
    # praw made F-barES-FEM-T4 look like it got worse with c, which is the
    # reverse of the paper's Fig. 6.
    #
    # psch is each scheme's OWN pressure, mapped to elements so the same
    # face-jump metric applies: K*theta per element for C3D4, K*theta_a
    # averaged back from the nodes for NS-FEM, K*(E A^c theta) averaged back
    # from the edges for F-barES-FEM-T4.  This is the field the paper plots.
    praw = Kb * theta
    psch = praw.copy()
    if scheme.startswith('fbar'):
        cyc = int(scheme.split('_')[-1]) if scheme.split('_')[-1].isdigit() else 1
        sel, idx, ekeys, edges, Vh, Ee, Sv, Rr = fbar_es_operators(
            nodes, tets, mat, g, vol, cyc, 1)
        psch[sel] = Kb * (Rr @ (Sv @ theta[sel]))
    elif scheme in ('ns_vol', 'ns_full', 'fs_ns', 'esdev_nsvol'):
        acc = {}
        for (a, mm), els in patch.items():
            if mm != 1:
                continue
            Va = sum(vol[e] for e in els) / 4.0
            acc[a] = sum(vol[e] / 4.0 * theta[e] for e in els) / Va
        for e in soft:
            psch[e] = Kb * np.mean([acc[int(a)] for a in tets[e]])

    def _osc(pf):
        pbar = np.abs(pf[soft]).mean()
        j = [abs(pf[els[0]] - pf[els[1]]) for f, els in faces.items()
             if len(els) == 2 and mat[els[0]] == 1 and mat[els[1]] == 1]
        j = np.array(j) if j else np.zeros(1)
        return ((j.mean() / pbar, j.max() / pbar) if pbar
                else (float('nan'), float('nan')))

    osc, oscmax = _osc(praw)
    oscs, _ = _osc(psch)
    return c1111, fl, osc, oscmax, len(nodes), len(tets), oscs


def runR(scheme, n, slab, Kb, Gb, stiff, bubble=False, stab=0.0,
         jitter=0.0):
    """R = C1111(undrained) / C1111(drained) of the soft phase, in miniature.

    It needs no external reference: R -> 1 means the scheme has lost the soft
    phase's bulk stiffness entirely, the failure mode seen on the real cells
    (R 6.29 -> 1.89)."""
    und = {0: stiff, 1: (Kb, Gb)}
    drn = {0: stiff, 1: (Kb / 1000.0, Gb)}
    cu, fl, osc, oscm, _, _, oscs = run(scheme, n, slab, und, stab=stab,
                                        bubble=bubble, jitter=jitter)
    # THE DENOMINATOR IS ALWAYS PLAIN C3D4, exactly as meshconv.py builds R:
    # Abaqus substitutes the hybrid element in the undrained cell only, so the
    # drained twin must stay the element both codes share.
    cd = run('c3d4', n, slab, drn, jitter=jitter)[0]
    return cu / cd, fl, osc, oscm, oscs


def main():
    Gb = 4.4e5
    Ei, ni = 9.37e9, 0.33
    stiff = (Ei / (3 * (1 - 2 * ni)), Ei / (2 * (1 + ni)))

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    global JIT
    JIT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    lo, hi = 0.375, 0.625
    print('stiff K/G = %.2f   slab x in [%.3f,%.3f]   n = %d   jitter = %.2f'
          % (stiff[0] / stiff[1], lo, hi, n, JIT))
    print('periodic cell, macroscopic eps_xx applied')

    exact = stiff[0] + 4.0 * stiff[1] / 3.0
    print('\n-- patch test, stiff phase only (exact C1111 = %.6e) --'
          % exact)
    for sc in ['c3d4', 'ns_vol', 'ns_full', 'fs_full', 'fs_ns', 'fs_vol']:
        for bub in (False, True):
            if bub and sc not in ('fs_full', 'fs_ns'):
                continue
            c, fl, _o, _om, _, _, _ = run(sc, 4, (2.0, 3.0), {0: stiff, 1: stiff},
                                       bubble=bub)
            print('   %-9s%-8s C1111=%.6e  rel err=%9.2e'
                  % (sc, ' +bubble' if bub else '', c, abs(c / exact - 1)))

    arms = [('c3d4', False), ('ns_vol', False), ('fs_full', False),
            ('fs_ns', False), ('fs_full', True), ('fs_ns', True)]
    print('\n-- R = C1111(undrained)/C1111(drained) vs bulk-to-shear ratio --')
    print('   the reference paper validates bFS-FEM at K/G = 25..100; the '
          'soft phase here is 5000')
    hdr = '   %-16s' % 'K/G' + ''.join('%11s' % ('%.0f' % r)
                                       for r in (50, 100, 500, 1000, 5000))
    print(hdr)
    for sc, bub in arms:
        tag = sc + (' +bub' if bub else '')
        row, extra = '   %-16s' % tag, []
        for ratio in (50, 100, 500, 1000, 5000):
            R, fl, osc, oscm = runR(sc, n, (lo, hi), ratio * Gb, Gb, stiff,
                                    bubble=bub, jitter=JIT)
            row += '%11.4f' % R
            extra.append((fl, osc))
        print(row)
        print('   %-16s' % '  fluc / p-osc'
              + ''.join('%11s' % ('%.1f/%.2f' % e) for e in extra))


if __name__ == '__main__':
    main()


# --------------------------------------------------------------------------
# F-barES-FEM-T4  (Onishi, Iida & Amaya, IJCM 15(7) 1845003, 2018)
#
# Small-strain reduction of eqs (1)-(11).  Two smoothing paths, deliberately
# different widths:
#
#   DEVIATORIC   edge-based, ONCE.  The smoothing domain of edge h is the set
#                of elements touching h, weighted V_e/6 (a T4 has 6 edges), so
#                sum_h V_h = sum_e V_e.                          [eqs 1-5]
#
#   VOLUMETRIC   J is smoothed c times through a node<->element cycle
#                  node:    J_n = sum_{e at n} J_e V_e/4 / sum V_e/4     (6)
#                  element: J_e = (1/4) sum_{n in e} J_n                 (7)
#                then edge-smoothed once                                 (8)
#                and recombined as F = J^(1/3) F~^iso                    (11)
#
# In small strain J - 1 = div u, so the cycle is a linear smoother on the
# element divergences and the whole thing is
#
#   eps_h = dev(eps~_h) + (1/3) dbar_h I,   dbar = E A^c D u
#
# with A = P Q the one-cycle element->node->element operator.  Energy is then
# G|dev eps~|^2 + (K/2) dbar^2 per edge domain.
#
# c = 0 recovers plain selective ES-FEM (edge deviatoric, edge volumetric).
# Smoothing never crosses a material interface, as everywhere else here.
def fbar_es_operators(nodes, tets, mat, g, vol, cycles, m, jinput=None):
    """Return (sel, idx, ekeys, edges, V_h, E, S, R) for material m.

    E is the elem->edge ES-FEM operator of eq. (1)/(8); S is the full
    volumetric chain of eqs. (6)-(8); R is the edge->elem average used
    both to recover the scheme's own element pressure and, under the
    'edge' reading of eq. (6), to seed the cyclic smoothing.

    jinput selects the reading of the tilde on the input to eq. (6):
      'elem' -- the raw element Jacobian det(eF).  S = E A^c.
      'edge' -- an element restriction of the already edge-smoothed
                J of eq. (5), which is what Fig. 1's ordering ('use
                ES-FEM once ... then use NS-FEM c times') suggests.
                S = E A^c R E.
    Defaults to FBAR_JIN, itself defaulting to 'elem'.
    """
    jinput = jinput or os.environ.get('FBAR_JIN', 'elem')
    sel = np.flatnonzero(mat == m)
    idx = {e: i for i, e in enumerate(sel)}
    ne = len(sel)

    # element -> node (weights V_e/4) and node -> element (1/4)
    rq, cq, vq = [], [], []
    Vn = {}
    for e in sel:
        for a in tets[e]:
            Vn[int(a)] = Vn.get(int(a), 0.0) + vol[e] / 4.0
    nid = {a: i for i, a in enumerate(sorted(Vn))}
    for e in sel:
        for a in tets[e]:
            rq.append(nid[int(a)]); cq.append(idx[e])
            vq.append(vol[e] / 4.0 / Vn[int(a)])
    Q = sp.coo_matrix((vq, (rq, cq)), shape=(len(nid), ne)).tocsr()
    rp, cp, vp = [], [], []
    for e in sel:
        for a in tets[e]:
            rp.append(idx[e]); cp.append(nid[int(a)]); vp.append(0.25)
    P = sp.coo_matrix((vp, (rp, cp)), shape=(ne, len(nid))).tocsr()

    # edges of this material
    edges = {}
    for e in sel:
        t = [int(x) for x in tets[e]]
        for i in range(4):
            for j in range(i + 1, 4):
                edges.setdefault((min(t[i], t[j]), max(t[i], t[j])), []).append(e)
    ekeys = sorted(edges)
    Vh = np.array([sum(vol[e] / 6.0 for e in edges[k]) for k in ekeys])
    re, ce, ve = [], [], []
    for h, k in enumerate(ekeys):
        for e in edges[k]:
            re.append(h); ce.append(idx[e]); ve.append(vol[e] / 6.0 / Vh[h])
    E = sp.coo_matrix((ve, (re, ce)), shape=(len(ekeys), ne)).tocsr()

    # edge -> element, the plain average over an element's 6 edges
    rr, cr, vr = [], [], []
    ehash = {k: h for h, k in enumerate(ekeys)}
    for e in sel:
        t = [int(x) for x in tets[e]]
        for i in range(4):
            for j in range(i + 1, 4):
                rr.append(idx[e]); cr.append(ehash[(min(t[i], t[j]),
                                                    max(t[i], t[j]))])
                vr.append(1.0 / 6.0)
    R = sp.coo_matrix((vr, (rr, cr)), shape=(ne, len(ekeys))).tocsr()

    # The volumetric chain of eqs. (6)-(8), as one operator on the element
    # divergence field.  'elem': seed with the raw element J, cycle, then
    # edge-smooth -> S = E A^c.  'edge': seed instead with the element
    # restriction of the edge-smoothed J -> S = E A^c R E.
    A = (P @ Q).tocsr()
    Ac = sp.identity(ne, format='csr')
    for _ in range(cycles):
        Ac = (Ac @ A).tocsr()
    S = (E @ Ac @ R @ E).tocsr() if jinput == 'edge' else (E @ Ac).tocsr()
    return sel, idx, ekeys, edges, Vh, E, S, R


# ==========================================================================
# F-barES-FEM-T4 AT FINITE STRAIN -- the paper's element, not its linearisation
# ==========================================================================
#
# Onishi, Iida & Amaya, IJCM 15(7) 1845003 (2018), eqs. (1)-(18), elastic
# (the Prony terms of eq. (13) dropped, g_i = 0, so T^dev = 2G H^dev).
# Read off the typeset equations, not a text extraction -- pdftotext drops
# every tilde and bar in this paper and the diacritics carry the meaning.
#
# Two readings that the glyphs settle, and that a text dump does not:
#
#   eq. (6)  the input to the first cycle is `Elem_e J`, with NO tilde: the
#            raw element Jacobian.  eq. (10) then states outright that the
#            whole chain is a weighted mean of the raw `Elem_e J`, which is
#            what makes it the single linear operator S = E A^c.
#   eq. (19) the push-forward is `Edge_h Ftilde^-1`, WITH a tilde, not Fbar^-1.
#            So the kinematics -- F~, the stretching D~, the push-forward and
#            the volume -- are the ES-FEM motion throughout, and the only
#            F-bar modification is that T is evaluated from the recombined F
#            of eq. (11).  That is de Souza Neto's F-bar exactly.
#
# The tangent of eq. (19) is "omitted here" in the paper (it forwards to
# Onishi et al. 2017, which we do not have), so deriving it would be the one
# place a mistake could hide.  It is not needed: the internal force is fully
# specified by eqs. (1)-(18), and a modified Newton iteration driven by the
# small-strain operator converges to the same displacement field that the
# consistent tangent would.  The converged answer is the paper's element; only
# the path to it differs.
class FbarNL:
    """Internal force of F-barES-FEM-T4(c), finite strain, Hencky elastic."""

    def __init__(self, nodes, tets, mat, g, vol, props, cycles):
        self.tets, self.g, self.props = tets, g, props
        self.nn = len(nodes)
        self.parts = []
        for m in sorted(set(int(x) for x in mat)):
            sel, idx, ekeys, edges, Vh, E, S, R = fbar_es_operators(
                nodes, tets, mat, g, vol, cycles, m)
            # (edge, element, weight) incidence of eq. (1): every element
            # hands V_e/6 to each of its 6 edges.
            hh, ee, ww = [], [], []
            Ec = E.tocoo()
            for h, j, w in zip(Ec.row, Ec.col, Ec.data):
                hh.append(h); ee.append(j); ww.append(w)
            hh = np.array(hh); ee = np.array(ee); ww = np.array(ww)
            gsel = g[sel][ee]                        # (L,4,3)
            nodesel = tets[sel][ee]                  # (L,4)
            dof = (nodesel[:, :, None] * 3 + np.arange(3)).ravel()
            part = {'sel': sel, 'E': E, 'S': S, 'Vh': Vh,
                    'K': props[m][0], 'G': props[m][1], 'hh': hh,
                    'ee': ee, 'ww': ww, 'gsel': gsel, 'dof': dof,
                    'gel': g[sel], 'nodesel': nodesel,
                    'tetsel': tets[sel], 'ne': len(sel),
                    'nh': E.shape[0], 'L': len(hh)}
            self._prepare_tangent(part)
            self.parts.append(part)

    def _prepare_tangent(self, p):
        """Build the u-independent operators of the consistent tangent.

        The tangent of `force` is the operator chain

            du -> dF -> dF~ -> (dJ, dJbar) -> dFbar -> dT -> dA -> df

        with one sparse matrix per arrow.  Everything that depends only on
        the reference mesh -- the gradients, the elem->edge smoothing, the
        support weights -- is built once here; the configuration-dependent
        blocks are rebuilt per call in `tangent`.

        Args:
            p: The part dictionary of one material.
        """
        ne, nh, nn = p['ne'], p['nh'], self.nn
        grid = np.indices((4, 3, 3)).reshape(3, -1)     # (a,k,l) of one block
        a_, k_, l_ = grid[0], grid[1], grid[2]
        # dF_e[i,j] = sum_a du[3 n_a + i] g_e[a,j]: (9 ne) x (3 nn)
        i_ = np.repeat(np.arange(ne), 12)
        a4_ = np.tile(np.repeat(np.arange(4), 3), ne)
        k4_ = np.tile(np.arange(3), 4 * ne)
        rows, cols, vals = [], [], []
        for j in range(3):
            rows.append(9 * i_ + 3 * k4_ + j)
            cols.append(3 * p['tetsel'][i_, a4_] + k4_)
            vals.append(p['gel'][i_, a4_, j])
        p['Gd'] = sp.coo_matrix(
            (np.concatenate(vals),
             (np.concatenate(rows), np.concatenate(cols))),
            shape=(9 * ne, 3 * nn)).tocsr()
        p['E9'] = sp.kron(p['E'], sp.identity(9), format='csr')
        # df = sum_he ww v g_e^T dA_h: (3 nn) x (9 nh), the adjoint of the
        # force accumulation, so no scatter loop is needed.
        e_ = np.repeat(np.arange(p['L']), 36)
        ar_ = np.tile(a_, p['L'])
        kr_ = np.tile(k_, p['L'])
        lr_ = np.tile(l_, p['L'])
        hh = p['hh'][e_]
        rows = 3 * p['nodesel'][e_, ar_] + lr_
        cols = 9 * hh + 3 * kr_ + lr_
        vals = p['ww'][e_] * p['Vh'][hh] * p['gsel'][e_, ar_, kr_]
        p['Fd'] = sp.coo_matrix((vals, (rows, cols)),
                                 shape=(3 * nn, 9 * nh)).tocsr()

    @staticmethod
    def _blockdiag9(blocks, nh):
        """Assemble (nh, 9, 9) dense blocks into a sparse (9 nh, 9 nh).

        Args:
            blocks: The blocks, shape (nh, 9, 9), indexed [h, row, col].
            nh: Number of blocks.

        Returns:
            The block-diagonal sparse matrix.
        """
        idx = np.indices((9, 9)).reshape(2, -1)
        rows = (9 * np.arange(nh)[:, None] + idx[0][None, :]).ravel()
        cols = (9 * np.arange(nh)[:, None] + idx[1][None, :]).ravel()
        return sp.coo_matrix((blocks.reshape(nh, 81).ravel(), (rows, cols)),
                              shape=(9 * nh, 9 * nh)).tocsr()

    def _configuration(self, p, F):
        """Per-edge kinematic and stress state at one configuration.

        Args:
            p: The part dictionary of one material.
            F: Element deformation gradients, shape (ne, 3, 3).

        Returns:
            The dictionary of per-edge arrays (F~, J~, Jbar, Fbar, V, w, T,
            F~^-1, A, v, alpha, beta).
        """
        sel, E, S, Vh = p['sel'], p['E'], p['S'], p['Vh']
        K, G = p['K'], p['G']
        Ftil = (E @ F.reshape(len(sel), 9)).reshape(-1, 3, 3)
        Jtil = np.linalg.det(Ftil)
        Jbar = S @ np.linalg.det(F)
        alpha = (Jbar / Jtil) ** (1.0 / 3.0)
        Fbar = alpha[:, None, None] * Ftil
        B = Fbar @ np.transpose(Fbar, (0, 2, 1))
        w, V = np.linalg.eigh(B)
        H = np.einsum('hij,hj,hkj->hik', V, 0.5 * np.log(w), V)
        trH = np.trace(H, axis1=1, axis2=2)
        T = (K * trH)[:, None, None] * np.eye(3) \
            + 2.0 * G * (H - (trH / 3.0)[:, None, None] * np.eye(3))
        Finv = np.linalg.inv(Ftil)
        return {'Ftil': Ftil, 'Jtil': Jtil, 'Jbar': Jbar, 'Fbar': Fbar,
                'V': V, 'w': w, 'T': T, 'Finv': Finv, 'A': Finv @ T,
                'v': Jtil * Vh, 'alpha': alpha,
                'beta': 1.0 / (3.0 * Jbar)}

    def _hencky_blocks(self, cfg, p):
        """The 9x9 blocks of dT/dFbar, one per edge.

        The tensor logarithm is differentiated in the eigenbasis of
        B = F~bar F~bar^T: with eigenvalues w_i and h_i = log(w_i)/2,

            dH_ij = f_ij dB_ij,   f_ij = (h_i - h_j)/(w_i - w_j),

        and f_ii = 1/(2 w_i), which is also the degenerate limit when two
        eigenvalues coincide (B = I at u = 0, for instance, where the
        divided difference is 0/0).  The nine columns are obtained by
        applying the map to the nine unit perturbations in one batched
        pass -- no hand-indexed fourth-order tensor.

        Args:
            cfg: The per-edge state from `_configuration`.
            p: The part dictionary of one material.

        Returns:
            The array of blocks, shape (nh, 9, 9).
        """
        Fbar, V, w = cfg['Fbar'], cfg['V'], cfg['w']
        nh = len(w)
        unit = np.zeros((9, 3, 3))
        for n in range(9):
            unit[n, n // 3, n % 3] = 1.0
        dF = unit[:, None, :, :]
        dB = dF @ np.transpose(Fbar, (0, 2, 1))[None] \
            + Fbar[None] @ np.transpose(dF, (0, 1, 3, 2))
        dBe = np.einsum('hki,nhkl,hlj->nhij', V, dB, V)
        h = 0.5 * np.log(w)
        dw = w[:, :, None] - w[:, None, :]
        dh = h[:, :, None] - h[:, None, :]
        tol = 1e-10 * np.abs(w).max(axis=1)[:, None, None]
        f = np.where(np.abs(dw) > tol, dh / np.where(dw == 0.0, 1.0, dw),
                     0.5 / w[:, None, :])
        dH = np.einsum('hik,nhkl,hjl->nhij', V,
                       f[None] * dBe, V)
        dtr = np.trace(dH, axis1=2, axis2=3)
        eye = np.eye(3)[None, None]
        dT = (p['K'] * dtr)[:, :, None, None] * eye \
            + 2.0 * p['G'] * (dH - (dtr / 3.0)[:, :, None, None] * eye)
        # dT[n] is the response to perturbation n, i.e. column n of the
        # block; transpose to [h, row, col] for `_blockdiag9`.
        return dT.transpose(1, 2, 3, 0).reshape(nh, 9, 9)

    def tangent(self, u):
        """Consistent tangent d f / d u of the finite-strain force.

        Chain rule through eqs. (1)-(17): the F-bar recombination of eq.
        (11), the Hencky tangent of eqs. (12)-(13) in the eigenbasis of
        B = F F^T, and the eq. (17) push-forward with the current volume
        J~ V_h.  Exact to round-off; no finite differences anywhere.

        Args:
            u: The displacement vector, length 3 * n_nodes.

        Returns:
            The tangent as a sparse (3 n_nodes, 3 n_nodes) matrix.
        """
        U = u[:3 * self.nn].reshape(-1, 3)
        F = np.einsum('eai,eaj->eij', U[self.tets], self.g) + np.eye(3)
        nn, kdof = self.nn, 3 * self.nn
        Ktan = sp.csr_matrix((kdof, kdof))
        for p in self.parts:
            cfg = self._configuration(p, F[p['sel']])
            Ftil, Finv, Jtil = cfg['Ftil'], cfg['Finv'], cfg['Jtil']
            nh, ne = cfg['Fbar'].shape[0], p['ne']
            # dJ_e = J_e F_e^{-T} : dF_e, as a (ne) x (9 ne) operator.
            # vec(dF)[3 n + m] = dF[n, m], and the contraction index is
            # F_e^{-1}[m, n] -- so the transpose is explicit here.
            n_, m_ = np.indices((3, 3)).reshape(2, -1)
            Fe = F[p['sel']]
            Finve = np.linalg.inv(Fe)
            cols = 9 * np.arange(ne)[:, None] + (3 * n_ + m_)[None, :]
            Jmap = sp.coo_matrix(
                ((np.linalg.det(Fe)[:, None] * Finve[:, m_, n_]).ravel(),
                 (np.repeat(np.arange(ne), 9), cols.ravel())),
                shape=(ne, 9 * ne)).tocsr()
            # dJbar = S dJ
            Sd = p['S'] @ Jmap @ p['Gd']
            # dFbar = alpha (dF~ - F~ (F~^-1:dF~)/3) + Fbar dJbar/(3 Jbar)
            fv = Ftil.reshape(nh, 9)
            # the contraction F~^-1 : dF~ reads vec(F~^-T), not vec(F~^-1)
            iv = np.transpose(Finv, (0, 2, 1)).reshape(nh, 9)
            eye9 = np.eye(9)[None]
            B1b = cfg['alpha'][:, None, None] * eye9 \
                - (cfg['alpha'] / 3.0)[:, None, None] * fv[:, :, None] \
                * iv[:, None, :]
            B1 = self._blockdiag9(B1b, nh)
            b2v = cfg['beta'][:, None] * cfg['Fbar'].reshape(nh, 9)
            B2 = sp.coo_matrix(
                (b2v.ravel(),
                 (np.repeat(9 * np.arange(nh), 9)
                  + np.tile(np.arange(9), nh), np.repeat(np.arange(nh), 9))),
                shape=(9 * nh, nh)).tocsr()
            Hmap = self._blockdiag9(self._hencky_blocks(cfg, p), nh)
            # dA = -F~^-1 dF~ A + F~^-1 dT.  With the row-major vec of
            # numpy, vec(M X N) = (M (x) N^T) vec(X), so the 9x9 blocks
            # are -F~^-1 (x) A^T and F~^-1 (x) I, assembled as
            # [h, i, k, j, l] -> row (i,k), column (j,l).
            AT = np.transpose(cfg['A'], (0, 2, 1))
            N1 = self._blockdiag9(
                (-Finv[:, :, None, :, None]
                 * AT[:, None, :, None, :]).reshape(nh, 9, 9), nh)
            N2 = self._blockdiag9(
                (Finv[:, :, None, :, None]
                 * np.eye(3)[None, None, :, None, :]).reshape(nh, 9, 9), nh)
            dFtil = p['E9'] @ p['Gd']
            dFbar = B1 @ dFtil + B2 @ Sd
            dA = N1 @ dFtil + N2 @ (Hmap @ dFbar)
            # The force also carries the CURRENT volume J~ V_h, whose
            # derivative is the second term: df = sum ww V_h [dJ~ (g A)
            # + J~ (g dA)].  At u = 0, T = 0 so that term vanishes and the
            # rest already is K_lin exactly; at finite strain both count.
            colsE = 9 * np.arange(nh)[:, None] + (3 * n_ + m_)[None, :]
            Td = sp.coo_matrix(
                ((Jtil[:, None] * Finv[:, m_, n_]).ravel(),
                 (np.repeat(np.arange(nh), 9), colsE.ravel())),
                shape=(nh, 9 * nh)).tocsr()
            gA = np.einsum('lak,lkc->lac', p['gsel'], cfg['A'][p['hh']])
            ac = np.indices((4, 3)).reshape(2, -1)
            e12_ = np.repeat(np.arange(p['L']), 12)
            a12_, c12_ = np.tile(ac[0], p['L']), np.tile(ac[1], p['L'])
            Fd_vol = sp.coo_matrix(
                ((p['ww'][e12_] * p['Vh'][p['hh'][e12_]]
                  * gA[e12_, a12_, c12_]).ravel(),
                 (3 * p['nodesel'][e12_, a12_] + c12_,
                  p['hh'][e12_])),
                shape=(3 * nn, nh)).tocsr()
            Dv = sp.diags(np.repeat(cfg['Jtil'], 9), format='csr')
            Ktan = Ktan + p['Fd'] @ (Dv @ dA) + Fd_vol @ (Td @ dFtil)
        return Ktan.tocsr()

    def force(self, u):
        U = u[:3 * self.nn].reshape(-1, 3)
        # eq. (1) input: the element deformation gradient, F = I + grad u
        F = np.einsum('eai,eaj->eij', U[self.tets], self.g)
        F += np.eye(3)
        J = np.linalg.det(F)
        f = np.zeros(3 * self.nn)
        for p in self.parts:
            sel, E, S, Vh = p['sel'], p['E'], p['S'], p['Vh']
            K, G = p['K'], p['G']
            # eq. (1): edge-smoothed F~
            Ftil = (E @ F[sel].reshape(len(sel), 9)).reshape(-1, 3, 3)
            Jtil = np.linalg.det(Ftil)                       # eq. (5)
            # eqs. (6)-(8): c cycles elem->node->elem, then elem->edge
            Jbar = S @ J[sel]
            # eqs. (4), (9), (11):  F = Jbar^(1/3) * Jtil^(-1/3) * F~
            Fbar = (Jbar / Jtil)[:, None, None] ** (1.0 / 3.0) * Ftil
            # Hencky strain of the recombined F, eq. (12) note
            B = Fbar @ np.transpose(Fbar, (0, 2, 1))
            w, V = np.linalg.eigh(B)
            H = np.einsum('hij,hj,hkj->hik', V, 0.5 * np.log(w), V)
            trH = np.trace(H, axis1=1, axis2=2)
            # eqs. (12)-(13), elastic: T = K tr(H) I + 2G dev(H)
            Tstr = (K * trH)[:, None, None] * np.eye(3) \
                + 2.0 * G * (H - (trH / 3.0)[:, None, None] * np.eye(3))
            # eq. (17): paired with D~, pushed forward by F~^-1, over the
            # current volume of the ES-FEM motion, v = J~ V^ini
            A = np.linalg.solve(Ftil, Tstr)                  # (F~^-1 T)[k,p]
            v = Jtil * Vh
            # f_{P:p} = N~^ini_{P,k} (F~^-1 T)_{kp} v, accumulated per (h,e)
            contrib = (p['ww'] * v[p['hh']])[:, None, None] \
                * (p['gsel'] @ A[p['hh']])
            f += np.bincount(p['dof'], weights=contrib.ravel(),
                             minlength=3 * self.nn)
        return f


def run_nl(scheme, n, slab, props, eps=1e-3, jitter=0.0, geomname='sphere',
           tol=1e-10, maxit=60, verbose=False, consistent=False, history=None,
           ninc=1, macro='xx'):
    """Periodic cell driven by macroscopic eps_xx, solved with the finite-strain
    F-barES-FEM-T4 internal force and a modified Newton iteration whose
    iteration matrix is the small-strain operator of the same scheme.

    With `consistent=True` the iteration matrix is instead the exact
    tangent `FbarNL.tangent` at the current iterate, refactored every
    iteration; the residual then falls quadratically.  `history`, if
    given, collects the relative residual of every iteration.  `ninc`
    applies `eps` in equal increments, which is what lets a 20-50 %
    macroscopic strain be reached at all; each increment still converges
    quadratically.  `macro` is 'xx' for uniaxial stretch or 'xz' for
    simple shear."""
    cyc = int(scheme.split('_')[-1])
    nodes, tets, mat = mesh_box(n, slab[0], slab[1], jitter,
                                geom=GEOM[geomname])
    g, vol = grads(nodes, tets)
    faces, patch = topology(tets, mat)
    Klin = assemble('fbar_%d' % cyc, nodes, tets, mat, g, vol, faces, patch,
                    props, 0.0, False, None)
    nl = FbarNL(nodes, tets, mat, g, vol, props, cyc)

    nn = len(nodes)
    ndof = Klin.shape[0]
    tolp = 1e-9
    base = np.mod(np.round(nodes * 1e9) / 1e9, 1.0)
    base[np.abs(base - 1.0) < tolp] = 0.0
    key = np.round(base * 1e7).astype(np.int64)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    first = np.full(len(uniq), -1, dtype=int)
    for i in range(nn):
        if first[inv[i]] < 0:
            first[inv[i]] = i
    master = first[inv]
    shift = nodes - nodes[master]
    nm = len(uniq)
    rows, cols, vals = [], [], []
    for i in range(nn):
        for c in range(3):
            rows.append(3 * i + c); cols.append(3 * inv[i] + c); vals.append(1.0)
    Tm = sp.coo_matrix((vals, (rows, cols)), shape=(ndof, 3 * nm)).tocsr()
    # macro = 'xx' is the uniaxial stretch eps_xx = eps; 'xz' is the simple
    # shear gamma_xz = eps (u_x = eps z), whose periodic images differ by
    # the z component of the shift.
    a = np.zeros(ndof)
    if macro == 'xx':
        a[0:3 * nn:3] = eps * shift[:, 0]
    elif macro == 'xz':
        a[0:3 * nn:3] = eps * shift[:, 2]
    else:
        raise ValueError("macro must be 'xx' or 'xz'")

    keep = np.ones(3 * nm, dtype=bool)
    keep[0:3] = False
    Kr = (Tm.T @ Klin @ Tm).tocsr()[keep][:, keep].tocsc()
    lu = spl.splu(Kr)

    # The macroscopic strain is applied in `ninc` equal increments.  A
    # single 50 % step from the affine field is simply too far for Newton
    # at this contrast -- the first correction tries to invert an element
    # -- and that is a property of the load step, not of the tangent.  The
    # last increment is where the quadratic convergence is measured.
    q = np.zeros(3 * nm)
    for inc in range(1, ninc + 1):
        a_inc = a * (inc / float(ninc))
        # normalise on the state the increment STARTS from, not on the bare
        # affine field: at a later increment the affine field alone can
        # already contain an inverted element, and it is not where the
        # iteration begins.
        f0 = np.linalg.norm((Tm.T @ nl.force(Tm @ q + a_inc))[keep]) or 1.0
        for it in range(maxit):
            u = Tm @ q + a_inc
            r = (Tm.T @ nl.force(u))[keep]
            rn = np.linalg.norm(r) / f0
            if history is not None:
                history.append(rn)
            if verbose:
                print('      inc %d it %2d  |r| = %.3e' % (inc, it, rn))
            if rn < tol:
                break
            alpha_print = 1.0
            if consistent:
                Kt = nl.tangent(u)
                Kr = (Tm.T @ Kt @ Tm).tocsr()[keep][:, keep].tocsc()
                lu = spl.splu(Kr)
            dq = lu.solve(r)
            alpha = 1.0
            if consistent:
                # A full Newton step can still overshoot near an inverted
                # element (Jbar/Jtil < 0), where the fractional power is
                # not defined.  Backtrack on the residual; at the solution
                # alpha = 1 is accepted and the convergence stays
                # quadratic.
                for _ in range(30):
                    q_try = q.copy()
                    q_try[keep] -= alpha * dq
                    try:
                        r_try = (Tm.T @ nl.force(Tm @ q_try + a_inc))[keep]
                        rn_try = np.linalg.norm(r_try) / f0
                    except (np.linalg.LinAlgError, ValueError):
                        rn_try = np.inf
                    if np.isfinite(rn_try) and rn_try < rn:
                        break
                    alpha *= 0.5
                alpha_print = alpha
            q[keep] -= alpha * dq
            if verbose:
                print('               alpha = %.2e' % alpha_print)
        else:
            raise RuntimeError('Newton did not converge, |r| = %.3e' % rn)

    u = Tm @ q + a
    c1111 = (a @ nl.force(u)) / eps ** 2
    un = u[:3 * nn]
    aff = np.zeros_like(un)
    if macro == 'xx':
        aff[0::3] = eps * nodes[:, 0]
    else:
        aff[0::3] = eps * nodes[:, 2]
    fl = (un - aff).reshape(-1, 3)
    fl = fl - fl.mean(axis=0)
    return c1111, np.abs(fl).max() / eps, it, rn
