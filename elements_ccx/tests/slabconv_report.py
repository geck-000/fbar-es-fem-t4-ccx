"""Mesh-convergence report + Richardson extrapolation for the slabconv cell.

Reads the CalculiX `.dat` files (plain C3D4 und/drn, F-barES-FEM-T4 c=0/c=1
und) and the Abaqus C3D4H CSV (`slabconv_abq.csv`, und/drn from Roihu), and
produces

  * a per-arm table of C1111 und/drn and R = C1111(und)/C1111(drn),
  * Richardson h->0 extrapolation of R and of C1111_und (linear tets, p=2),
  * a figure of C1111 und and drn against h, with the extrapolated limits.

The drained twin is plain C3D4 in both codes, so the denominator of R is the
same mesh/geometry everywhere and only the undrained element differs.

    python3 analysis/slabconv_report.py [root] [n ...]
"""
import os
import re
import sys

import numpy as np

ROOT = 'out_slabconv/kg500_one_xsym'
NS = [10, 20, 30, 40, 50, 60]
EPS = 1.0e-3
AX = 0  # load along x -> reaction component RF1


def faces(dat, ax):
    """(C1111, equilibrium residual) from a ccx `.dat` driven-face reaction."""
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


def read_ccx(n, arm):
    """C1111 from out_slabconv/<arm>/m_ccx.dat for arm in {und,drn,und_fbar0,und_fbar1}."""
    dat = os.path.join(ROOT, 'n%d' % n, arm, 'm_ccx.dat')
    c, _ = faces(dat, AX)
    return c


def read_abq():
    """{ (n, state): C1111 } from the Roihu C3D4H extraction CSV."""
    out = {}
    with open(os.path.join(ROOT, 'slabconv_abq.csv')) as f:
        for line in f:
            if line.startswith('n,'):
                continue
            n, state, c = line.strip().split(',')[:3]
            out[(int(n), state)] = float(c)
    return out


def richardson(xs, ys, p=2):
    """Fit y = y_inf + C h^p (h = 1/n) by least squares over the last k points.

    Returns (y_inf, C, residual).  The two-point fit is exact and is what a
    two-point Richardson step reduces to; the three-point fit has one degree of
    freedom and exposes whether the tail is actually at the assumed rate p.
    """
    h = 1.0 / np.asarray(xs, dtype=float)
    A = np.column_stack([np.ones(len(xs)), h ** p])
    coef, *_ = np.linalg.lstsq(A, np.asarray(ys, dtype=float), rcond=None)
    y_inf, C = coef
    resid = np.abs(A @ coef - np.asarray(ys)).max()
    return y_inf, C, resid


def main():
    global ROOT, NS
    if len(sys.argv) > 1:
        ROOT = sys.argv[1]
    if len(sys.argv) > 2:
        NS = [int(a) for a in sys.argv[2:]]

    abq = read_abq()

    arms = {
        'C3D4H (abq)': lambda n, st: abq.get((n, st), float('nan')),
        'F-bar c=1': lambda n, st: read_ccx(n, 'und_fbar1') if st == 'und'
        else read_ccx(n, 'drn'),
        'F-bar c=0': lambda n, st: read_ccx(n, 'und_fbar0') if st == 'und'
        else read_ccx(n, 'drn'),
        'C3D4 (ccx)': lambda n, st: read_ccx(n, 'und' if st == 'und' else 'drn'),
    }

    print('ROOT %s, load x, eps %.1e, R = C1111(und)/C1111(drn) [drn = plain C3D4]' % (ROOT, EPS))
    print()
    hdr = '%-6s %-8s' + ' %12s' * len(arms) + ' %12s'
    print(hdr % (('n', 'el/slab') + tuple(arms.keys()) + ('R',)))
    series = {a: [] for a in arms}
    rseries = {a: [] for a in arms}
    for n in NS:
        els = 0.2 * n
        row = []
        for a, fn in arms.items():
            u, d = fn(n, 'und'), fn(n, 'drn')
            series[a].append(u)
            r = u / d if d and d == d else float('nan')
            rseries[a].append(r)
            row.append(u)
        # R values
        rs = [rseries[a][-1] for a in arms]
        print(('%-6d %-8.0f' + ' %12.4e' * len(arms) + ' %12.4e')
              % ((n, els) + tuple(row) + (rs[0],)))
        # absolute und table
    print()
    print('R series')
    print('%-6s' % 'n' + ''.join(' %14s' % a for a in arms))
    for i, n in enumerate(NS):
        print('%-6d' % n + ''.join(' %14.4f' % rseries[a][i] for a in arms))

    print()
    print('Richardson extrapolation h -> 0 (h = 1/n, linear-tet rate p = 2)')
    print('%-12s %-10s %-12s %-12s %-10s' % ('arm', 'pts', 'R_inf', 'C', 'maxresid'))
    for a in arms:
        ys = rseries[a]
        # use the finest points that are finite
        idx = [i for i, y in enumerate(ys) if y == y and y > 0]
        fin = idx[-3:]
        if len(fin) >= 3:
            x = [NS[i] for i in fin]
            y = [ys[i] for i in fin]
            y_inf, C, resid = richardson(x, y)
            print('%-12s %-10s %-12.4f %-12.3g %-10.2e' % (a, '+'.join(map(str, x)), y_inf, C, resid))
        else:
            print('%-12s  (insufficient points)' % a)

    print()
    print('agreement of C1111_und at the finest point (n = %d)' % NS[-1])
    ref = series['C3D4H (abq)'][-1]
    for a in arms:
        v = series[a][-1]
        print('  %-12s %12.4e  rel vs C3D4H %+.3f%%' % (a, v, 100.0 * (v - ref) / ref))

    # ---------------- figure ----------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
        import figstyle
        figstyle.apply()

        h = 1.0 / np.asarray(NS, dtype=float)
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        colors = {'C3D4H (abq)': figstyle.VERM, 'F-bar c=1': figstyle.BLUE,
                  'F-bar c=0': figstyle.PURPLE, 'C3D4 (ccx)': figstyle.ORANGE}
        marks = {'C3D4H (abq)': 'o', 'F-bar c=1': 's', 'F-bar c=0': '^',
                 'C3D4 (ccx)': 'D'}
        for st, ax, ti in (('und', axes[0], 'undrained (locking)'),
                           ('drn', axes[1], 'drained')):
            for a in arms:
                y = [arms[a](n, st) for n in NS]
                ax.plot(h, y, marker=marks[a], color=colors[a], label=a,
                        ls='-', ms=6)
                # extrapolated limit for the und panel only
                if st == 'und':
                    idx = [i for i, v in enumerate(y) if v == v]
                    if len(idx) >= 3:
                        fin = idx[-3:]
                        y_inf, _, _ = richardson([NS[i] for i in fin], [y[i] for i in fin])
                        ax.axhline(y_inf, color=colors[a], ls=':', lw=1.0, alpha=0.6)
            ax.set_xlabel('h = 1/n')
            ax.set_title(ti)
            ax.invert_xaxis()
        axes[0].set_ylabel('C1111 [Pa]')
        axes[1].legend(fontsize=10, loc='best')
        axes[0].legend(fontsize=10, loc='best')
        fig.suptitle('slabconv cell: F-barES-FEM-T4 vs Abaqus C3D4H, identical meshes')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            out = os.path.join(ROOT, 'slabconv_convergence.%s' % ext)
            fig.savefig(out)
        print()
        print('wrote %s/slabconv_convergence.{png,pdf}' % ROOT)
    except Exception as e:  # pragma: no cover - plotting is best-effort
        print('figure skipped: %s' % e)


if __name__ == '__main__':
    main()
