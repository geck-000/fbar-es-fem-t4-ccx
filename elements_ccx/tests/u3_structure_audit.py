#!/usr/bin/env python3
"""Measure, for a converted ccx deck, the three quantities that decide whether
the 'direct U3 assembly' idea can help:

  1. transient  = number of insert() calls mastruct.c makes for U3 elements
                  (what fills mast1/next and overflows the 32-bit ITG)
  2. plan_alloc = sum(3*nring * 3*nope), the size the assembly plan proposes
                  to preallocate as the "compact" structure
  3. final      = the true deduplicated U3 dof-pair count (what jq/irow hold
                  after mastruct.c's sort + 'removing duplicate entries' pass)

If plan_alloc >= transient, a per-element "direct" build saves nothing: the
duplication is across elements, not inside one.

usage: u3_structure_audit.py <deck.inp>
"""
import re, sys, numpy as np

f = sys.argv[1] if len(sys.argv) > 1 else 'Job-LCOL4-utx-ccx.inp'

types = {}
for line in open(f):
    if line.startswith('*USER ELEMENT'):
        m = re.search(r'TYPE=(\w+).*NODES=(\d+)', line)
        if m and m.group(1).startswith('U3'):
            nm = m.group(1)
            types[nm] = (ord(nm[2]) - ord('A') + 1, int(m.group(2)))

conn = {}; cur = None; toks = []
def flush():
    if cur and toks:
        conn.setdefault(cur, []).extend(toks)
for line in open(f):
    if line[0] == '*':
        flush(); toks = []
        m = re.match(r'\*ELEMENT,TYPE=(\w+)', line.strip())
        cur = m.group(1) if m and m.group(1) in types else None
        continue
    if cur:
        toks.extend(line.split(','))
flush()

M = 1 << 21
keys = []; ringnodes = []
nel = transient = plan_alloc = 0
for t, tk in conn.items():
    r, s = types[t]
    a = np.array(tk, dtype=np.int64).reshape(-1, s + 1)
    c = a.shape[0]; nel += c
    rr, ss = 3 * r, 3 * s
    transient  += c * (rr * ss - rr * (rr - 1) // 2)   # the jj/ll loop
    plan_alloc += c * rr * ss                          # plan section 3.1 step 2
    node = a[:, 1:]; ring = node[:, :r]
    ringnodes.append(ring.ravel())
    p = np.broadcast_to(ring[:, :, None], (c, r, s)).reshape(-1)
    q = np.broadcast_to(node[:, None, :], (c, r, s)).reshape(-1)
    keys.append(np.unique(np.minimum(p, q) * M + np.maximum(p, q)))
    del a, node, ring, p, q

allk = np.unique(np.concatenate(keys))
lo, hi = allk // M, allk % M
ndist, nself = int((lo != hi).sum()), int((lo == hi).sum())
final = 9 * ndist + 3 * nself

print("deck                    : %s" % f)
print("U3 elements             : %d" % nel)
print("transient insert() calls: %.4e" % transient)
print("plan preallocation      : %.4e  (%.2fx the transient)"
      % (plan_alloc, plan_alloc / transient))
print("true final structure    : %.4e  (%.1fx smaller than transient)"
      % (final, transient / final))
