"""Report slabconv.sh over finished output, without re-solving.

    report_slabconv.py <root> <K/G> <load axis> <n>...

Split out of the shell driver for the same reason report_abaqus_ratio.py was:
a sweep that takes an hour must not have to be repeated to re-read it, and
editing the reporting must not be able to disturb a running solve.
"""
import os
import re
import sys

EPS = 1.0e-3


def faces(dat, ax):
    """(C1111, equilibrium residual) from the driven face reactions."""
    if not os.path.isfile(dat) or os.path.getsize(dat) == 0:
        return float('nan'), float('nan')
    blocks, cur, kind = {}, None, None
    for ln in open(dat):
        m = re.match(r'\s*(forces|displacements).*for set (\S+)', ln)
        if m:
            kind, cur = m.group(1), m.group(2).upper()
            blocks.setdefault((kind, cur), [])
            continue
        if cur and ln.strip() and not ln.lstrip().startswith(('*', 'S T')):
            f = ln.split()
            if len(f) >= 4:
                try:
                    blocks[(kind, cur)].append([float(x) for x in f[1:4]])
                except ValueError:
                    cur = None
    f1 = blocks.get(('forces', 'X1'), [])
    f0 = blocks.get(('forces', 'X0'), [])
    if not f1:
        return float('nan'), float('nan')
    r1 = sum(r[ax] for r in f1)
    bal = abs(r1 + sum(r[ax] for r in f0)) / abs(r1) if r1 else float('nan')
    return r1 / EPS, bal


def abq_faces(dat, ax):
    """(C1111, equilibrium residual) from an ABAQUS .dat.

    `*NODE PRINT, NSET=X1 / U,RF` prints one table per set with six value
    columns -- U1 U2 U3 RF1 RF2 RF3 -- so the driven reaction is column 3+ax.
    The rows are summed here rather than trusting Abaqus's own TOTAL line,
    which is not emitted for every variable combination.
    """
    if not os.path.isfile(dat) or os.path.getsize(dat) == 0:
        return float('nan'), float('nan')
    tot, cur = {}, None
    for ln in open(dat):
        m = re.search(r'NODE SET\s+(\S+)', ln)
        if m:
            cur = m.group(1).upper()
            tot.setdefault(cur, [0.0] * 6)
            continue
        if cur is None:
            continue
        f = ln.split()
        # a data row: node number then six values.  Abaqus interleaves
        # MAXIMUM/MINIMUM/TOTAL summary lines, which do not start with an int.
        if len(f) >= 7 and f[0].isdigit():
            try:
                v = [float(x) for x in f[1:7]]
            except ValueError:
                continue
            tot[cur] = [a + b for a, b in zip(tot[cur], v)]
    if 'X1' not in tot:
        return float('nan'), float('nan')
    r1 = tot['X1'][3 + ax]
    r0 = tot.get('X0', [0.0] * 6)[3 + ax]
    bal = abs(r1 + r0) / abs(r1) if r1 else float('nan')
    return r1 / EPS, bal


def main():
    root, kg, load = sys.argv[1], sys.argv[2], sys.argv[3]
    ns = [a for a in sys.argv[4:] if not a.startswith('--')]
    abqroot = next((a.split('=', 1)[1] for a in sys.argv[4:]
                    if a.startswith('--abq=')), None)
    ax = {'x': 0, 'y': 1, 'z': 2}[load]
    print('\nK/G = %s   load along %s   R = C1111(und)/C1111(drn), '
          'denominator always plain C3D4' % (kg, load))
    cols = ['n', 'el/slab', 'C1111 drn', 'und C3D4', 'und fbar c=1',
            'R C3D4', 'R fbar', 'C3D4 exc']
    if abqroot:
        cols[5:5] = ['und C3D4H']
        cols += ['R C3D4H']
    print(('%-4s %-8s' + ' %13s' * (3 + (1 if abqroot else 0))
           + ' %8s %8s %9s' + (' %8s' if abqroot else '')) % tuple(cols))
    for n in ns:
        b = os.path.join(root, 'n' + n)
        d, bd = faces(os.path.join(b, 'drn', 'm_ccx.dat'), ax)
        u, bu = faces(os.path.join(b, 'und', 'm_ccx.dat'), ax)
        fb, bf = faces(os.path.join(b, 'und_fbar1', 'm_ccx.dat'), ax)
        # How much stiffer plain C3D4 reads than the unlocked element on the
        # SAME mesh: the locking, isolated.  It is the quantity that has to
        # fall as h -> 0 if both elements are consistent.
        exc = 100.0 * (u - fb) / fb if fb == fb and fb else float('nan')
        vals = [n, 0.2 * int(n), d, u, fb]
        fmt = '%-4s %-8.0f %13.6e %13.6e %13.6e'
        ah = float('nan')
        if abqroot:
            # Abaqus C3D4H on the SAME mesh -- make_slabconv.py emits
            # m_abq.inp beside m_ccx.inp from one geometry, so this column is
            # a fourth element on identical nodes, not a foreign reference.
            ab = os.path.join(abqroot, 'n' + n)
            ah, bh = abq_faces(os.path.join(ab, 'und', 'm_abq.dat'), ax)
            vals.append(ah)
            fmt += ' %13.6e'
        vals += [u / d if d else float('nan'), fb / d if d else float('nan'),
                 exc]
        fmt += ' %8.4f %8.4f %8.2f%%'
        if abqroot:
            vals.append(ah / d if d and ah == ah else float('nan'))
            fmt += ' %8.4f'
        print(fmt % tuple(vals))
        for tag, v in (('drn', bd), ('und', bu), ('fbar', bf)):
            if v == v and v > 1e-6:
                print('     WARNING %s: the two faces carry %.2e of the driven '
                      'reaction out of balance' % (tag, v))


if __name__ == '__main__':
    main()
