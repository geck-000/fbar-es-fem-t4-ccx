# -*- coding: utf-8 -*-
"""Extract C1111_eff from a slabconv Abaqus ODB via the driven-face reaction.

    abaqus python slabconv_extract.py <odb> <n> <state> <out_csv>

The slabconv deck drives the X1 node set by eps = 1e-3 along x and fixes its
mate X0, so

    C1111_eff = sum RF_x(X1) / eps

and the equilibrium residual is |sum RF_x(X1) + sum RF_x(X0)| / |sum RF_x(X1)|.
This is the Abaqus-side complement to report_slabconv.py's abq_faces(), which
reads the same quantities from the *NODE PRINT block in the .dat; the ODB route
is what the other Roihu campaigns use, so it is the one run there.

The deck is a flat (no-*PART) input, so Abaqus folds every *NSET under a single
default instance; the sets are found there rather than on the assembly.
"""
import os
import sys

from odbAccess import openOdb

EPS = 1.0e-3
AX = 0  # load is along x -> reaction component RF1


def _node_set(odb, name):
    """The node set called `name`, searching instances then the assembly."""
    asm = odb.rootAssembly
    for key in asm.nodeSets.keys():
        if key.strip().upper() == name:
            return asm.nodeSets[key]
    for iname in asm.instances.keys():
        for key in asm.instances[iname].nodeSets.keys():
            if key.strip().upper() == name:
                return asm.instances[iname].nodeSets[key]
    return None


def main():
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    odb_path, n, state, out_csv = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

    odb = openOdb(path=odb_path, readOnly=True)
    step = odb.steps[odb.steps.keys()[-1]]
    frame = step.frames[-1]

    def sum_rf(name):
        region = _node_set(odb, name)
        if region is None:
            return float('nan')
        vals = frame.fieldOutputs['RF'].getSubset(region=region).values
        return sum(v.data[AX] for v in vals)

    r1 = sum_rf('X1')
    r0 = sum_rf('X0')
    c1111 = r1 / EPS
    bal = abs(r1 + r0) / abs(r1) if r1 else float('nan')

    new = not os.path.exists(out_csv)
    with open(out_csv, 'a') as f:
        if new:
            f.write('n,state,C1111,RF_X1,RF_X0,balance\n')
        f.write('%s,%s,%.12e,%.12e,%.12e,%.12e\n'
                % (n, state, c1111, r1, r0, bal))

    print('n=%s %s: C1111=%.6e  RF_X1=%.6e  balance=%.3g'
          % (n, state, c1111, r1, bal))
    odb.close()


if __name__ == '__main__':
    main()
