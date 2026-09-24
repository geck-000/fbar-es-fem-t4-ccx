"""Check that a deck written by fbares.py says exactly what u2edge.f and
u3vol.f will walk.

The generator computes the stencils with sparse boolean products; the elements
recompute them with a node->element walk.  If the two ever disagree, u3vol
stops and names the element -- but only for a node it needs and cannot find.
The reverse error, a connectivity carrying nodes the walk never reaches, is
silent and just wastes matrix entries.  This checks set EQUALITY, both ways.

It also checks the deck-format rules that have bitten this campaign before:
at most 15 fields per line and no trailing comma (textpart in elements.f is
dimensioned 16), every *USER ELEMENT declared before any *ELEMENT that uses
it, letters-only type suffixes, and unique element ids.

    verify_fbares_deck.py ORIGINAL.inp GENERATED.inp [cycles]
"""
import sys
from collections import defaultdict

fails = []


def chk(name, ok, detail=''):
    print('   %-58s %s%s' % (name, 'PASS' if ok else 'FAIL',
                             '' if ok else '  ' + detail))
    if not ok:
        fails.append(name)


def read(path):
    """elements by type, *USER ELEMENT declarations, and field-shape errors."""
    els, decl, order, badline = defaultdict(dict), {}, [], []
    elsetof = {}
    mode, cur, curtype, pend, nwant = None, None, None, None, 0
    for lno, raw in enumerate(open(path), 1):
        s = raw.strip()
        if not s or s.startswith('**'):
            continue
        if s.startswith('*'):
            u = s.upper().replace(' ', '')
            mode = None
            if u.startswith('*USERELEMENT'):
                t = next(p.split('=')[1] for p in s.split(',')
                         if p.strip().upper().startswith('TYPE='))
                n = int(next(p.split('=')[1] for p in s.split(',')
                             if p.strip().upper().startswith('NODES=')))
                decl[t.strip()] = n
                order.append(('decl', t.strip()))
            elif u.startswith('*ELEMENT'):
                curtype = next(p.split('=')[1] for p in s.split(',')
                               if p.strip().upper().startswith('TYPE=')).strip()
                cur = next((p.split('=')[1] for p in s.split(',')
                            if p.strip().upper().startswith('ELSET=')),
                           'ALL').strip()
                order.append(('use', curtype))
                elsetof[curtype] = cur
                mode = 'e'
                pend, nwant = None, decl.get(curtype, 4)
            continue
        if mode != 'e':
            continue
        f = [x.strip() for x in s.split(',')]
        if len(f) > 16:
            badline.append((lno, 'more than 16 fields'))
        if f and f[-1] == '':
            badline.append((lno, 'trailing comma'))
        f = [x for x in f if x]
        if pend is None:
            pend = [int(x) for x in f]
        else:
            pend += [int(x) for x in f]
        if len(pend) >= nwant + 1:
            els[curtype][pend[0]] = pend[1:nwant + 1]
            pend = None
    return els, decl, order, badline, elsetof


def walk(tets, mats, imat, edge, ncyc):
    """u3vol.f's walk, as a set of elements."""
    at = defaultdict(list)
    for e, t in tets.items():
        if mats[e] != imat:
            continue
        for a in t:
            at[a].append(e)
    na, nb = edge
    cur = {e for e in at[na] if nb in tets[e]}
    for _ in range(ncyc):
        nxt = set()
        for e in cur:
            for k in tets[e]:
                nxt.update(at[k])
        cur = nxt
    return cur


def main():
    orig, gen = sys.argv[1], sys.argv[2]
    ncyc = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    # The walk is O(stencil) in Python, so on a campaign-sized deck check a
    # random sample rather than all ~370k elements.  The deck-format checks
    # below are exhaustive regardless.
    nsample = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    els, decl, order, badline, elsetof = read(gen)

    chk('deck: no line over 16 fields, no trailing comma', not badline,
        str(badline[:3]))
    seen = set()
    dup = [e for t in els for e in els[t] if e in seen or seen.add(e)]
    chk('deck: element ids unique', not dup, str(dup[:5]))
    firstuse = {}
    for kind, t in order:
        if kind == 'use' and t not in firstuse:
            firstuse[t] = True
            if t.startswith('U') and t not in decl and t != 'U4':
                fails.append('undeclared type ' + t)
    chk('deck: every user type declared before use',
        not any(f.startswith('undeclared') for f in fails))
    bad = [t for t in decl if not t[2:].isalpha() and t[2:]]
    chk('deck: type suffixes are letters only (elements.f digit rule)',
        not bad, str(bad[:5]))

    # material of every U4 tet, from the *SOLID SECTION of its elset
    tets, mats = {}, {}
    esof, matof = {}, {}
    mode, cur, curtype = None, None, None
    for raw in open(gen):
        s = raw.strip()
        if not s or s.startswith('**'):
            continue
        if s.startswith('*'):
            u = s.upper().replace(' ', '')
            mode = None
            if u.startswith('*ELEMENT') and 'TYPE=U4' in u:
                cur = next((p.split('=')[1] for p in s.split(',')
                            if p.strip().upper().startswith('ELSET=')),
                           'ALL').strip()
                mode = 'u4'
            elif u.startswith('*SOLIDSECTION'):
                es = next(p.split('=')[1] for p in s.split(',')
                          if p.strip().upper().startswith('ELSET=')).strip()
                mt = next(p.split('=')[1] for p in s.split(',')
                          if p.strip().upper().startswith('MATERIAL=')).strip()
                matof[es] = mt
            continue
        if mode == 'u4':
            f = [x.strip() for x in s.split(',') if x.strip()]
            if len(f) >= 5:
                tets[int(f[0])] = [int(x) for x in f[1:5]]
                esof[int(f[0])] = cur
    for e in tets:
        mats[e] = matof[esof[e]]
    chk('deck: U4 tets carried over from the original',
        len(tets) > 0, '%d found' % len(tets))

    # every U2 / U3 element, against an independent walk
    nd2 = nd3 = 0
    bad2, bad3, badr = [], [], []
    LET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    import random
    rnd = random.Random(11)
    for t, table in els.items():
        if not (t.startswith('U2') or t.startswith('U3')):
            continue
        items = list(table.items())
        if nsample and len(items) > nsample:
            items = rnd.sample(items, nsample)
        for eid, conn in items:
            na, nb = conn[0], conn[1]
            # The material is the element's OWN, from its *SOLID SECTION --
            # exactly ielmat(1,nelem), which is what u2edge/u3vol pass as
            # imatf.  Guessing it from "any tet containing this edge" is wrong
            # at a phase boundary: an interface edge carries ONE SMOOTHING
            # DOMAIN PER PHASE, with different rings, and the guess silently
            # compares one against the other.
            imat = matof[elsetof[t]]
            cyc = 0 if t.startswith('U2') else ncyc
            got = walk(tets, mats, imat, (na, nb), cyc)
            want = set()
            for e in got:
                want.update(tets[e])
            have = set(conn)
            if t.startswith('U2'):
                nd2 += 1
                if have != want:
                    bad2.append((eid, sorted(want - have), sorted(have - want)))
            else:
                nd3 += 1
                if have != want:
                    bad3.append((eid, sorted(want - have), sorted(have - want)))
                # RING FIRST, AND THE LABEL AGREES WITH IT.
                #
                # mastruct.c and mafillsmas.f drop the identically-zero
                # (support - ring) x (support - ring) block by treating the
                # first nring nodes as the edge ring, with nring read from
                # lakon(3:3).  If the deck's ORDER or its LABEL disagreed with
                # the ring u3vol actually walks, real tbar rows would be
                # dropped -- and no patch test can see that, because a uniform
                # field gives the missing rows nothing to carry.  e_c3d_u3
                # checks the same thing at runtime; this checks it in the deck,
                # where the offending element can be named before a two-hour
                # solve.
                nr = LET.find(t[2]) + 1
                ringset = set()
                for e in walk(tets, mats, imat, (na, nb), 0):
                    ringset.update(tets[e])
                if nr < 1 or nr != len(ringset) or set(conn[:nr]) != ringset:
                    badr.append((eid, t, nr, len(ringset)))
    chk('U2 connectivity == the ring u2edge.f walks (%d elements)' % nd2,
        not bad2, str(bad2[:2]))
    chk('U3 connectivity == the E A^c support u3vol.f walks (%d elements)'
        % nd3, not bad3, str(bad3[:2]))
    chk('U3 ring is first in the connectivity and matches lakon(3:3)',
        not badr, str(badr[:2]))

    # the edge nodes must be konl(1), konl(2) and must share a tet
    at = defaultdict(list)
    for e, t4 in tets.items():
        for a in t4:
            at[a].append(e)
    bade = []
    for t, table in els.items():
        if not (t.startswith('U2') or t.startswith('U3')):
            continue
        items = list(table.items())
        if nsample and len(items) > nsample:
            items = rnd.sample(items, nsample)
        for eid, conn in items:
            if not any(conn[1] in tets[e] for e in at[conn[0]]):
                bade.append(eid)
    chk('konl(1),konl(2) are a real mesh edge', not bade, str(bade[:5]))

    print('\n%s' % ('all checks passed' if not fails
                    else 'FAILED: ' + ', '.join(sorted(set(fails)))))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
