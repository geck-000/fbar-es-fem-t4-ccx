#!/usr/bin/env python3
"""fbares.py -- rewrite a C3D4 deck as F-barES-FEM-T4.

Onishi, Iida & Amaya, Int. J. Comput. Methods 15(7) 1845003 (2018).  This
script is only the deck side: it builds the smoothing domains and rewrites the
deck, while the element operators live in the Fortran sources.

    K = K_dev   per EDGE (U2): V_h Bt^T D_dev Bt          eq. (1), (4), (13)
      + K_vol   per EDGE (U3): (K V_h) tbar^T sbar        eq. (6)-(11), (17)

The base tets stay in the deck retyped to U4 -- the null base tetrahedron --
contributing nothing: U2 and U3 carry the whole stiffness and both read the
tets for geometry through the 'U4' node->element map.

WHY THE CONNECTIVITY IS WRITTEN OUT AT ALL.  The elements recompute their own
weights from geometry -- which subsets of a ring form tets is not recoverable
from a node list -- but ccx has to know the stencil to size the element and to
build the matrix structure, so the node SET must be in the deck and must match
what u3vol walks exactly.  u3vol stops and names the element if it meets a
node the connectivity does not carry, so a mismatch is loud.

RING x SUPPORT: THE U3 ELEMENT MATRIX IS NOT DENSE OVER ITS STENCIL.

K_vol = (K V_h) tbar^T sbar is rank one, and its two factors have DIFFERENT
supports.  tbar is the UNSMOOTHED edge divergence of eq. (1) -- u3vol builds
it only from the tets that contain the edge -- so its support is exactly the
U2 ring, ~6.3 nodes.  Only sbar carries the wide E A^c support, ~33.7 nodes at
c = 1.  So s(ii,jj) can be nonzero only for i in the RING and j in the
SUPPORT: the (support - ring) x (support - ring) block is identically zero.

That block is most of the element, and mastruct.c was allocating structure for
all of it.  The connectivity is therefore written RING FIRST, the ring size is
carried in the type label, and mastruct.c / mafillsmas.f skip the outer block.
Measured at c = 1: insertions fall 3.3x and a finer production mesh
becomes buildable.
Nothing about the assembled matrix changes -- the entries dropped are zero.

STENCIL WIDTH IS THE BINDING CONSTRAINT.  Measured on the soft phase of a
production mesh (117437 tets, 36323 nodes, 184572 edges):

    U2 deviatoric      mean  6.3 nodes/edge   max  14   ->   42 DOF
    U3 volumetric c=1  mean 33.7              max 173   ->  519 DOF
    U3 volumetric c=2  mean 93.6              max 494   -> 1482 DOF

c=2 cannot be expressed at all: ccx carries a user element's node count in the
single byte lakon(8:8) and userelements.f rejects NODES > 255.  This script
refuses it rather than letting the solver truncate.  c=1 needs the element
matrix widened from 150 to 520 DOF along the whole e_c3d_u* path.
"""
import argparse
import os
import sys
from array import array
from collections import defaultdict

import numpy as np
import scipy.sparse as sp

LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
# Two suffix characters, LETTERS ONLY.  ccx's built-in element dispatch keys on
# digits in the label (elements.f: label(4:4).eq.'4' -> nope=4, '10' -> 10,
# '20' -> 20), so a type named U214 is claimed by the nope=4 rule before the
# *USER ELEMENT lookup and ccx reads 4 nodes instead of 14.  This showed up
# on only 3 of ~36000 elements.
# THE REAL CEILING ON c IS mastruct, NOT MEMORY AND NOT THE 255-NODE LIMIT.
#
# mastruct.c pushes one entry per off-diagonal (dof,dof) pair of the upper
# triangle of EVERY element onto `mast1`, with NO deduplication -- insert.c
# is a linked-list append, and the compression happens afterwards.  It costs
# 8 bytes each (mast1 and next, one ITG apiece) and the index is a 32-bit ITG
# in a stock build, so past 2^31 the 1.1x growth in insert.c overflows and
# ccx dies in u_realloc with a NEGATIVE allocation size.
#
# With the ring x support reduction above the count is
#
#     sum over edges of  T(3 n_h) - T(3 (n_h - r_h)) + T(3 r_h),  T(d)=d(d+1)/2
#
# still QUADRATIC in the stencil, so this stays the real ceiling on c --
# reached before memory is, and before the 255-node limit.
INS_LIMIT = 2 ** 31

NAMES = [a + b for a in LETTERS for b in LETTERS]


def final_structure(u2, u3, tcn, symmetric=False):
    """Exact number of distinct off-diagonal (dof,dof) pairs the matrix
    structure ends up holding.

    This is the UNION of the per-element patterns.  The INS_LIMIT count above
    is their SUM, which is ~20-25x larger because the same pair is produced by
    many different edges.  A ccx built with the deduplicating mastruct backend
    (CCX_MASTRUCT_DEDUP=chunk) never materialises that sum, so the union is
    the only thing that has to fit -- both under 2**31 and in RAM.

    Every U2 pattern (ring x ring) sits inside the U3 pattern of the same edge
    (ring x support, ring subset of support), so with R the edge-by-node ring
    incidence and S the edge-by-node support incidence, the F-bar part of the
    union is just the sparsity pattern of R^T S.  In the symmetric pairing the
    element matrix is dense over the support, so the pattern is S^T S instead
    and the ring is irrelevant.  Plain C3D4 elements outside the treated sets
    are in no stencil at all, so their own 4-cliques (T^T T) are unioned in on
    top.

    This counts the element patterns only.  ccx additionally expands DOFs that
    carry an SPC or MPC, and drops the constrained ones, which moves entries
    around.  Measured on an undrained production deck this lands 2% high --
    1.67e7 here against the 16 378 633 ccx goes on to report -- so treat it as
    a tight upper estimate rather than the exact final nnz.
    """
    # u2/u3 arrive as CSR (edge offsets + flat global node indices), which is
    # already the incidence pattern this needs -- no per-edge work at all.
    ring_e, ring_n, supp_e, supp_n = [], [], [], []
    h = 0
    for es in u3:
        ka, _, rptr, ridx = u2[es]
        _, _, sptr, sidx = u3[es]
        ne = len(ka)
        base = np.arange(h, h + ne, dtype=np.int64)
        ring_e.append(np.repeat(base, np.diff(rptr)))
        ring_n.append(ridx)
        supp_e.append(np.repeat(base, np.diff(sptr)))
        supp_n.append(sidx)
        h += ne
    if symmetric:
        ring_e, ring_n = supp_e, supp_n

    T = tcn

    nn = 1
    for a in (supp_n, ring_n):
        if a:
            nn = max(nn, max(int(x.max()) for x in a if x.size))
    if T.size:
        nn = max(nn, int(T.max()))
    nn += 1

    def inc(rows, cols, nr):
        rows = np.concatenate(rows) if rows else np.empty(0, dtype=np.int64)
        cols = (np.concatenate(cols).astype(np.int64, copy=False) if cols
                else np.empty(0, dtype=np.int64))
        return sp.csr_matrix((np.ones(rows.size, dtype=np.int32), (rows, cols)),
                             shape=(nr, nn))

    A = None
    if h:
        R = inc(ring_e, ring_n, h)
        S = inc(supp_e, supp_n, h)
        A = (R.T @ S).tocsr()
    if T.size:
        te = np.repeat(np.arange(T.shape[0], dtype=np.int64), 4)
        Ti = inc([te], [T.reshape(-1)], T.shape[0])
        P = (Ti.T @ Ti).tocsr()
        A = P if A is None else A + P
    if A is None:
        return 0

    # unordered pairs: symmetrise the pattern, then split diagonal from
    # off-diagonal (each off-diagonal pair is stored twice)
    A.data[:] = 1
    B = A + A.T
    B.data[:] = 1
    ndiag = int(B.diagonal().astype(bool).sum())
    noff = (B.nnz - ndiag) // 2
    # a node pair i!=j contributes a full 3x3 block (9 dof pairs); a node with
    # itself contributes the 3 off-diagonal entries of its own block
    return int(9 * noff + 3 * ndiag)


def _rstrip_lines(path):
    """The source deck line by line, newline stripped, without ever holding
    all of it -- read().splitlines() on a 400 MB deck is ~0.7 GB of str."""
    with open(path) as f:
        for ln in f:
            yield ln.rstrip('\n').rstrip('\r')


def iskw(line):
    s = line.lstrip()
    return s.startswith('*') and not s.startswith('**')


def parse(path):
    """tet ids, tet connectivity, elset membership, elset -> material,
    material -> (E, nu).

    Coordinates are deliberately NOT kept: the stencils are pure connectivity
    and nothing downstream reads a node position, so holding them costs a few
    hundred MB on a fine cell for nothing.

    The elements come back as two packed arrays -- ids (n,) and connectivity
    (n, 4) -- rather than a dict of tuples.  On the 0.0060 cell that dict is
    ~2.7M boxed ids and ~11M boxed node numbers, several GB, to carry data
    that fits in 100 MB flat.
    """
    tid, tcn = array('q'), array('q')
    elset_of, mat_of, elastic = defaultdict(lambda: array('q')), {}, {}
    mode, cur, curmat = None, None, None
    for ln in open(path):
        s = ln.strip()
        if not s or s.startswith('**'):
            continue
        if iskw(ln):
            u = s.upper().replace(' ', '')
            mode = None
            if u.startswith('*NODE') and 'OUTPUT' not in u:
                mode = 'n'
            elif u.startswith('*ELEMENT'):
                if 'TYPE=C3D4' in u:
                    mode = 'e'
                    cur = next((p.split('=')[1] for p in s.split(',')
                                if p.strip().upper().startswith('ELSET=')),
                               'ALL').strip()
            elif u.startswith('*SOLIDSECTION'):
                es = next((p.split('=')[1] for p in s.split(',')
                           if p.strip().upper().startswith('ELSET=')), None)
                mt = next((p.split('=')[1] for p in s.split(',')
                           if p.strip().upper().startswith('MATERIAL=')), None)
                if es and mt:
                    mat_of[es.strip()] = mt.strip()
            elif u.startswith('*MATERIAL'):
                curmat = next((p.split('=')[1] for p in s.split(',')
                               if p.strip().upper().startswith('NAME=')),
                              None)
                if curmat:
                    curmat = curmat.strip()
            elif u.startswith('*ELASTIC'):
                mode = 'el'
            continue
        f = [x.strip() for x in s.split(',') if x.strip()]
        if mode == 'e' and len(f) >= 5:
            e = int(f[0])
            tid.append(e)
            tcn.append(int(f[1]))
            tcn.append(int(f[2]))
            tcn.append(int(f[3]))
            tcn.append(int(f[4]))
            elset_of[cur].append(e)
        elif mode == 'el' and curmat and len(f) >= 2:
            elastic.setdefault(curmat, (float(f[0]), float(f[1])))
            mode = None
    tid = np.frombuffer(tid, dtype=np.int64).copy()
    tcn = np.frombuffer(tcn, dtype=np.int64).copy().reshape(-1, 4)
    elset_of = {k: np.frombuffer(v, dtype=np.int64).copy()
                for k, v in elset_of.items()}
    return tid, tcn, elset_of, mat_of, elastic


def stencils(conn, ncyc, chunk=20000):
    """Edge list and, per edge, the U2 ring and the U3 = E A^c node support.

    Built with boolean sparse products rather than a Python walk: on the real
    cells this is ~500k edges and the walk is the whole runtime.  The PATTERN
    is all that is needed -- the elements recompute the weights themselves --
    and verify_u8_chain.py has already checked that the weighted walk agrees
    with the prototype operator to 1e-15.

    conn is (ne, 4) LOCAL node indices for one material's tets.
    """
    ne = len(conn)
    nn = int(conn.max()) + 1 if ne else 0
    inc = sp.coo_matrix((np.ones(4 * ne, dtype=np.int8),
                         (np.repeat(np.arange(ne), 4), conn.ravel())),
                        shape=(ne, nn)).tocsr()
    inc.data[:] = 1

    # edges of this material, and the elements at each.
    #
    # Vectorised: a dict keyed by (int, int) tuples costs a few hundred bytes
    # per edge and a Python loop over 6*ne of them, which at 0.0060 is ~4M
    # edges and ~1 GB of dict before anything else is allocated.  Combining each
    # sorted pair into one int64 and calling np.unique gives the same edge
    # list in the same order -- lexicographic by (lo, hi), since hi < nn.
    ii = np.array([0, 0, 0, 1, 1, 2])
    jj = np.array([1, 2, 3, 2, 3, 3])
    ea = conn[:, ii].ravel()
    eb = conn[:, jj].ravel()
    lo = np.minimum(ea, eb)
    hi = np.maximum(ea, eb)
    upair, inv = np.unique(lo * nn + hi, return_inverse=True)
    del ea, eb, lo, hi
    keys = list(zip((upair // nn).tolist(), (upair % nn).tolist()))
    rows = inv.ravel()
    cols = np.repeat(np.arange(ne, dtype=np.int64), 6)
    E = sp.coo_matrix((np.ones(len(cols), dtype=np.int8), (rows, cols)),
                      shape=(len(keys), ne)).tocsr()
    E.data[:] = 1

    # A = P Q as a PATTERN: elements sharing a node
    A = (inc @ inc.T).tocsr()
    A.data[:] = 1

    # Returned FLAT, in CSR form (offsets + one index array), not as a list
    # of one small array per edge: at 0.0060 that list is ~4M ndarray objects
    # per side, whose ~112 bytes of object header each cost more than the
    # node indices they carry.
    rcnt, ridx, scnt, sidx = [], [], [], []
    for lo in range(0, len(keys), chunk):
        blk = E[lo:lo + chunk]
        r = (blk @ inc).tocsr()          # U2: nodes of the edge's own tets
        r.data[:] = 1
        rcnt.append(np.diff(r.indptr))
        ridx.append(r.indices[:r.indptr[-1]])
        s = blk
        for _ in range(ncyc):
            s = (s @ A).tocsr()
            s.data[:] = 1
        w = (s @ inc).tocsr()            # U3: nodes of the E A^c support
        w.data[:] = 1
        scnt.append(np.diff(w.indptr))
        sidx.append(w.indices[:w.indptr[-1]])

    def csr(cnt, idx):
        cnt = np.concatenate(cnt) if cnt else np.zeros(0, dtype=np.int64)
        ptr = np.zeros(len(cnt) + 1, dtype=np.int64)
        np.cumsum(cnt, out=ptr[1:])
        return ptr, (np.concatenate(idx) if idx
                     else np.zeros(0, dtype=np.int64))

    ka = np.array([k[0] for k in keys], dtype=np.int64)
    kb = np.array([k[1] for k in keys], dtype=np.int64)
    rptr, ridx = csr(rcnt, ridx)
    sptr, sidx = csr(scnt, sidx)
    return ka, kb, rptr, ridx, sptr, sidx


def card(eid, conn_nodes):
    """One *ELEMENT card, wrapped.

    At most 15 fields per line and NO trailing comma.  textpart in elements.f
    is dimensioned (16); a trailing comma adds a field and overflows it, after
    which ccx reads the continuation as a fresh element and reports
    'element N is already defined'.  ccx continues on node count instead.
    """
    row = [str(eid)] + [str(x) for x in conn_nodes]
    return [','.join(row[c:c + 15]) for c in range(0, len(row), 15)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('dst')
    ap.add_argument('--elset', action='append', default=None,
                    help='element set to treat; repeat for more than one. '
                         'Default: every C3D4 set in the deck.')
    ap.add_argument('--cycles', type=int,
                    default=int(os.environ.get('CCX_FBAR_C', 1)),
                    help='c, the number of cyclic smoothings of J, eq. (6)-(7)')
    ap.add_argument('--solver', default=os.environ.get('FBAR_SOLVER',
                                                       'PARDISO'))
    ap.add_argument('--symmetric', action='store_true',
                    help='gate the structure for the Galerkin pairing '
                         '(dense over the support); run ccx with '
                         'CCX_FBAR_SYM=1')
    ap.add_argument('--audit', action='store_true',
                    help='print both structure counts and exit')
    ap.add_argument('--nlgeom', action='store_true',
                    help='write *STATIC, NLGEOM: the U3 internal force is the '
                         'finite-strain one and the step is nonlinear')
    a = ap.parse_args()

    tid, tcn, elset_of, mat_of, elastic = parse(a.src)
    sets = a.elset or sorted(elset_of)
    for es in sets:
        if es not in elset_of:
            raise SystemExit('fbares: no C3D4 elements in ELSET=%s '
                             '(deck has %s)' % (es, sorted(elset_of)))
        if es not in mat_of:
            raise SystemExit('fbares: ELSET=%s has no *SOLID SECTION' % es)

    if a.cycles < 0 or a.cycles > 3:
        raise SystemExit('fbares: --cycles must be 0..3')

    maxel = int(tid.max()) if tid.size else 0
    tid_order = np.argsort(tid)
    tid_sorted = tid[tid_order]

    # --- build the stencils, per material -------------------------------
    #
    # PER MATERIAL, always.  The smoothing of eqs. (1) and (6)-(8) must not
    # cross a phase boundary: averaging a divergence across a 1000x modulus
    # contrast is not the method, and u6patch.f records that a shared patch
    # took equilibrium_gap from 2.5e-3 to 2.6e-1 on the finer layered cell.
    # An interface edge legitimately gets one smoothing domain per phase.
    u2, u3 = {}, {}
    stats = {}
    for es in sets:
        # local numbering is FIRST-ENCOUNTER, walking elements in ascending
        # id and nodes in connectivity order -- np.unique(return_index) gives
        # each value's first position, so sorting by that reproduces it
        # exactly, without the dict or the 4*ne interpreter iterations.
        ids = np.sort(elset_of[es])
        raw = tcn[tid_order[np.searchsorted(tid_sorted, ids)]]
        flat = raw.ravel()
        uniq, first, inv = np.unique(flat, return_index=True,
                                     return_inverse=True)
        o = np.argsort(first)
        # global node ids fit int32 comfortably; these arrays are the biggest
        # thing the stencils carry, so do not spend 8 bytes on each
        back = uniq[o].astype(np.int32)
        rank = np.empty(uniq.size, dtype=np.int64)
        rank[o] = np.arange(uniq.size, dtype=np.int64)
        conn = rank[inv.ravel()].reshape(-1, 4)
        del raw, flat, uniq, first, inv, o, rank
        ka, kb, rptr, ridx, sptr, sidx = stencils(conn, a.cycles)
        # one fancy-index over the whole flat array, not one per edge
        u2[es] = (back[ka], back[kb], rptr, back[ridx])
        u3[es] = (back[ka], back[kb], sptr, back[sidx])
        rs = np.diff(rptr)
        ss = np.diff(sptr)
        stats[es] = (len(ids), len(back), len(ka), rs, ss)

    # --- the mastruct insertion wall ------------------------------------
    #
    # mastruct.c builds the matrix structure by pushing ONE entry per
    # (dof, dof) pair of every element onto `mast1` and compressing
    # afterwards, and its counter is a 32-bit ITG in a stock build.
    #
    # A U3 element is NOT dense over its whole stencil, and exploiting that
    # is what makes c = 1 reachable on a production cell -- see RING x
    # SUPPORT at the top of this file.  Only the first nring rows of the
    # element matrix can be nonzero, so the deck costs
    #
    #     sum over edges of  T(3 n_h) - T(3 (n_h - r_h))    (U3)
    #                      + T(3 r_h)                       (U2)
    #
    # with T(d) = d(d+1)/2, instead of T(3 n_h) + T(3 r_h) for the dense
    # block.  Still quadratic in the stencil, and still the real ceiling on
    # c, but a factor ~3.3 lower at c = 1.
    # THE GATE IS ON THE FINAL STRUCTURE, ALWAYS.
    #
    # There used to be a second path here that gated on the TRANSIENT insertion
    # count -- the sum above, what legacy mastruct.c pushes onto mast1 before
    # sorting. That path is gone, and with it the --direct-structure flag that
    # selected between them, because the transient count is a property of a
    # backend nobody should now be building against: the deduplicating
    # mastruct (CCX_MASTRUCT_DEDUP=chunk) streams the transient away and never
    # materialises it, so gating on it rejected decks that build perfectly.
    #
    # Concretely it rejected 2.5-element undrained cells at 5.20e9
    # transient pairs against a 2.15e9 limit, when their
    # deduplicated structure fits comfortably. Keeping a flag whose default
    # answer was wrong for the supported build is worse than having no flag.
    nnz = final_structure(u2, u3, tcn, a.symmetric)
    if a.audit:
        print('AUDIT  PG structure         %.4e pairs' % nnz)
        print('AUDIT  symmetric structure  %.4e pairs'
              % final_structure(u2, u3, tcn, True))
        return
    if nnz >= INS_LIMIT:
        raise SystemExit(
            'fbares: the deduplicated structure of this deck holds %.2e '
            '(dof,dof) pairs, past the %.2e a 32-bit ITG can index.  Use '
            '--cycles %d or coarsen the mesh.'
            % (nnz, float(INS_LIMIT), max(a.cycles - 1, 0)))

    # --- the 255-node wall ----------------------------------------------
    worst = max((s[4].max(), es) for es, s in stats.items())
    if worst[0] > 255:
        raise SystemExit(
            'fbares: the c=%d volumetric stencil reaches %d nodes in '
            'ELSET=%s. ccx stores a user element\'s node count in the single '
            'byte lakon(8:8) and userelements.f rejects NODES > 255, so this '
            'cannot be expressed as an element at all. Use --cycles 1, or '
            'assemble the volumetric term outside the element loop (a '
            'direct global-assembly pass).'
            % (a.cycles, worst[0], worst[1]))

    # --- connectivity, ORDERED: ring first ------------------------------
    #
    # nl ALREADY contains both edge nodes -- it is the node set of the tets at
    # the edge (U2) or of the E A^c support (U3), and both contain the edge
    # itself.  Declaring len(nl)+2 made ccx read two fields past the end of
    # every card and swallow the next element's id, which showed up as
    # duplicate ids and a connectivity carrying its own successor.
    #
    # For U3 the ORDER now matters as well as the set: the first nring nodes
    # must be exactly the edge ring, because that is what mastruct.c and
    # mafillsmas.f use to skip the identically-zero outer block.  Getting it
    # wrong is not silent -- e_c3d_u3 checks that tbar vanishes past nring and
    # stops with the element number if it does not.
    groups = {}

    def push(key, row):
        g = groups.get(key)
        if g is None:
            g = groups[key] = array('i')
        g.extend(row)

    for es in sets:
        ka, kb, rptr, ridx = u2[es]
        _, _, sptr, sidx = u3[es]
        for h in range(len(ka)):
            na = int(ka[h])
            nb = int(kb[h])
            rl = ridx[rptr[h]:rptr[h + 1]]
            sl = sidx[sptr[h]:sptr[h + 1]]
            rest = sorted(int(x) for x in rl if x != na and x != nb)
            push(('U2', es, len(rl), len(rl)), [na, nb] + rest)
            ring = set(int(x) for x in rl)
            supp = set(int(x) for x in sl)
            if not ring <= supp:
                raise SystemExit(
                    'fbares: the edge ring of %d-%d is not contained in its '
                    'E A^c support -- the ring-first ordering the element '
                    'relies on is not well defined' % (na, nb))
            inner = [na, nb] + sorted(ring - {na, nb})
            outer = sorted(supp - ring)
            push(('U3', es, len(supp), len(inner)), inner + outer)

    # --- type names -----------------------------------------------------
    #
    # U2:  U2<xy>              xy just disambiguates (kind, elset, size).
    # U3:  U3<r><xy>           r ENCODES nring: LETTERS[nring-1], read back by
    #                          mastruct.c as lakon(3:3) and by e_c3d_u3 and
    #                          mafillsmas.f the same way.  elements.f keys the
    #                          *USER ELEMENT lookup on label(2:5), so 'U' plus
    #                          four characters is all ccx can carry and three
    #                          suffix letters is the most that fits.
    suffix = {}
    for k in sorted(x for x in groups if x[0] == 'U2'):
        i = len(suffix)
        if i >= len(NAMES):
            raise SystemExit('fbares: more U2 (elset, size) groups than the '
                             '%d type names available' % len(NAMES))
        suffix[k] = NAMES[i]
    per = defaultdict(int)
    sym_names = {}
    for k in sorted(x for x in groups if x[0] == 'U3'):
        if a.symmetric:
            # '_' is the Galerkin marker read by mastruct.c, mafillsmas.f and
            # e_c3d_u3.f: no ring x support reduction, the element is dense
            # over the support.  One type per support size, merged across
            # elsets (the material is read per element), named by the
            # two-letter table.
            nr = k[2]
            if nr not in sym_names:
                if len(sym_names) >= len(NAMES):
                    raise SystemExit('fbares: more than %d support sizes'
                                     % len(NAMES))
                sym_names[nr] = '_' + NAMES[len(sym_names)]
            suffix[k] = sym_names[nr]
            continue
        nr = k[3]
        if not 1 <= nr <= len(LETTERS):
            raise SystemExit(
                'fbares: an edge ring spans %d nodes, and the ring size is '
                'carried in ONE label character (lakon(3:3), A..Z = 1..%d) so '
                'that mastruct.c can skip the zero outer block. Beyond that '
                'the reduction cannot be expressed; remesh, or drop the '
                'reduction and accept the dense structure.'
                % (nr, len(LETTERS)))
        j = per[nr]
        per[nr] += 1
        if j >= len(NAMES):
            raise SystemExit('fbares: more than %d U3 (elset, size) groups at '
                             'ring size %d' % (len(NAMES), nr))
        suffix[k] = LETTERS[nr - 1] + NAMES[j]

    # --- emit -------------------------------------------------------------
    fh = open(a.dst, 'w')

    def emit(s):
        fh.write(s)
        fh.write('\n')

    u5decl, done_step = False, False
    # two streaming passes over the source rather than holding the whole deck
    # as a list of lines: on the 0.0060 cell that list is ~0.7 GB
    drop, sawnp = False, any(
        iskw(l) and l.upper().replace(' ', '').startswith(('*NODEPRINT',
                                                           '*NODEFILE'))
        for l in open(a.src))
    eid = maxel
    for ln in _rstrip_lines(a.src):
        if not iskw(ln):
            # a data line of a dropped block goes with it
            if drop:
                continue
        else:
            u = ln.upper().replace(' ', '')
            # NO ELEMENT IN AN F-bar DECK CARRIES STRESS.  The base tets are
            # U4 and return a null matrix; U2 and U3 are smoothing domains
            # with no shape function of their own and write nothing to stx.
            # Left in, *EL PRINT would report a column of exact zeros as if
            # it were the answer.  Read the result from displacements and
            # reactions instead -- which is why a *NODE PRINT is added below
            # if the deck has none.
            drop = u.startswith('*ELPRINT') or u.startswith('*ELFILE')
            if drop:
                emit('** FBAR fbares: dropped '
                           + ln.split(',')[0].strip()
                           + ' -- no element in an F-bar deck carries stress')
                continue
            if u.startswith('*ENDSTEP') and not sawnp:
                sawnp = True
                # *NODE FILE, not *NODE PRINT: printing needs an NSET and
                # ccx does not define NALL itself (cgx does), so a
                # *NODE PRINT,NSET=NALL is answered with 'node set NALL does
                # not exist' and an EMPTY .dat.  *NODE FILE needs no set and
                # writes the .frd the post-processing already reads.
                emit('*NODE FILE')
                emit('U,RF')
            if u.startswith('*ELEMENT') and any(
                    ('ELSET=' + e).upper() in u for e in sets):
                if not u5decl:
                    emit('*USER ELEMENT,TYPE=U4,NODES=4,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3')
                    u5decl = True
                emit(ln.replace('C3D4', 'U4').replace('c3d4', 'U4'))
                continue
            if u.startswith('*STATIC'):
                # The volumetric operator of eq. (17) is NOT symmetric, so the
                # assembled matrix is not either and incomplete-Cholesky PCG
                # does not apply.  PARDISO takes it through the asymmetric
                # path (nasym -> mafillsmas -> mtype=11).
                ln = '*STATIC, SOLVER=' + a.solver
            if u.startswith('*STEP') and a.nlgeom:
                # CalculiX carries NLGEOM on the *STEP card, not on
                # *STATIC: a NLGEOM parameter on *STATIC is read and
                # silently ignored with a warning.
                ln = '*STEP, NLGEOM'
            if u.startswith('*STEP') and not done_step:
                done_step = True
                # ALL *USER ELEMENT declarations first, then all *ELEMENT
                # blocks.  ccx sorts the deck into per-keyword chains
                # (keystart.f: *USER ELEMENT is position 3, *ELEMENT is 4), so
                # interleaving them puts a multi-line element's continuation
                # at a chain boundary, where ccx reads it as a fresh element
                # and reports 'element N is already defined'.
                emit('** FBAR F-barES-FEM-T4(c=%d): %d edge domains'
                           % (a.cycles, sum(len(u2[e][0]) for e in sets)))
                seen_types = set()
                for k in sorted(groups):
                    key = (k[0], suffix[k])
                    if key in seen_types:
                        continue
                    seen_types.add(key)
                    emit('*USER ELEMENT,TYPE=%s%s,NODES=%d,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3'
                               % (k[0], suffix[k], k[2]))
                for k in sorted(groups):
                    kind, es, sz, nr = k
                    tag = 'FBD' if kind == 'U2' else 'FBV'
                    emit('*ELEMENT,TYPE=%s%s,ELSET=%s_%s'
                               % (kind, suffix[k], tag, es))
                    g, w = groups[k], k[2]
                    for off in range(0, len(g), w):
                        conn = g[off:off + w]
                        eid += 1
                        # konl(1), konl(2) ARE THE EDGE NODES, in both U2 and
                        # U3; u2edge/u3vol identify the edge from them and
                        # find the tets themselves.  For U3, konl(1..nring)
                        # is the edge ring -- see above.
                        for _c in card(eid, conn):
                            emit(_c)
                # Their OWN elset + *SOLID SECTION, never the phase's: an
                # element in the phase elset joins its *EL PRINT set, and
                # printoutelem.f would try to integrate a smoothing domain
                # that has no material volume of its own -- and it would
                # corrupt the volume-averaged stress the homogenisation reads.
                for k in sorted(groups):
                    kind, es, sz, nr = k
                    tag = 'FBD' if kind == 'U2' else 'FBV'
                    emit('*SOLID SECTION,ELSET=%s_%s,MATERIAL=%s'
                               % (tag, es, mat_of[es]))
        emit(ln)

    fh.close()

    print('fbares: %s -> %s   (c = %d)' % (a.src, a.dst, a.cycles))
    for es in sets:
        nel, nnd, ned, rs, ss = stats[es]
        E, nu = elastic.get(mat_of[es], (float('nan'), float('nan')))
        kg = (1.0 / (3.0 * (1.0 - 2.0 * nu))) / (1.0 / (2.0 * (1.0 + nu)))
        print('  %-14s %7d tets  %6d nodes  %7d edges   material %-14s '
              'K/G = %.0f' % (es, nel, nnd, ned, mat_of[es], kg))
        print('      U2 ring   mean %5.1f  p99 %4d  max %4d nodes '
              '(%d DOF)' % (rs.mean(), np.percentile(rs, 99), rs.max(),
                            3 * (rs.max() + 0)))
        print('      U3 stencil mean %5.1f  p99 %4d  max %4d nodes '
              '(%d DOF)' % (ss.mean(), np.percentile(ss, 99), ss.max(),
                            3 * (ss.max() + 0)))
    ndof = 3 * max(s[4].max() for s in stats.values())
    print('  %d element types; widest element %d DOF' % (len(groups), ndof))
    print('  deduplicated structure %.2e pairs = %.1f GB of irow -- '
          'run ccx with CCX_MASTRUCT_DEDUP=chunk'
          % (nnz, nnz * 4.0 / 2.0 ** 30))
    if ndof > 765:
        print('  NOTE: the e_c3d_u* family and mafillsm.f hold 765 DOF '
              '(255 nodes, the lakon(8:8) encoding limit).  %d DOF cannot '
              'be a ccx user element at all.' % ndof)
    print('  run with CCX_FBAR_C=%d -- the same value the deck was '
          'generated with' % a.cycles)
    if a.symmetric:
        print('  the U3 label carries "_": the symmetric (Galerkin) '
              'volumetric pairing, dense over each support, solved on the '
              'symmetric path.  No environment switch is needed.')


if __name__ == '__main__':
    main()
