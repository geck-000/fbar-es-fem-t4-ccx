#!/bin/bash -l
# Is the CalculiX-vs-Abaqus gap on rve_gas.csv a mesh-density difference?
#
# The regenerated GAS cells carry about twice the elements the campaign's did
# -- 46462 for GAS_v00 here against the ~20700 matrix elements recorded for the
# same run in validation_ccx/results/results_scf.csv. A porous cell gets MORE COMPLIANT as the
# mesh is refined (a coarse mesh cannot resolve the stress concentration around
# a void, and over-stiffens), so a denser mesh would push E_eff down -- the
# right sign for the 0.5-3.6% deficit, and it grows with void content, which is
# also what was observed.
#
# The test holds the geometry EXACTLY fixed with SPAX_SAVE_PACKING /
# SPAX_LOAD_PACKING -- the same packed spheres remeshed at several L_mesh -- so
# nothing varies but the discretisation. If E_eff climbs toward the stored
# Abaqus value as L_mesh coarsens toward theirs, the gap is mesh, not solver.
set -eu
# Resolved before the cd below, while a relative $0 still points here.
VALDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# This validation drives SpaX: it generates cells with SpaX_Standalone.py, solves
# them through SpaX_CalculiX.py, and compares against stored Abaqus results. SpaX
# is a separate repository (https://github.com/geck-000/SpaX), so point SPAX_ROOT
# at a checkout of it. There is no sensible default from inside this repository.
if [ -z "${SPAX_ROOT:-}" ]; then
    echo "ERROR: set SPAX_ROOT to a SpaX checkout (https://github.com/geck-000/SpaX)." >&2
    exit 2
fi
if [ ! -f "$SPAX_ROOT/SpaX_Standalone.py" ]; then
    echo "ERROR: \$SPAX_ROOT=$SPAX_ROOT has no SpaX_Standalone.py." >&2
    exit 2
fi
cd "$SPAX_ROOT"

PY=${PY:-python3}
CCX=${CCX:-ccx_spax}
ROOT=${ROOT:-out_gasmesh}
ROW=${ROW:-GAS_v10}
MESHES=${MESHES:-"0.055 0.045 0.033 0.026"}
CPUS=${CPUS:-8}

export SPAX_CCX="$CCX"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_MESH_ORDER=2
export SPAX_SEED=${SPAX_SEED:-11}
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT/pack"
head -1 $VALDIR/params/rve_gas.csv > "$ROOT/row.csv"
grep "^$ROW," $VALDIR/params/rve_gas.csv >> "$ROOT/row.csv"
$PY - "$ROOT/row.csv" <<'PYEOF'
import csv, sys
p = sys.argv[1]
rows = list(csv.DictReader(open(p)))
fn = list(rows[0].keys())
rows[0]['Mode2'] = ''          # only uniaxial X is read here
w = csv.DictWriter(open(p, 'w', newline=''), fieldnames=fn)
w.writeheader(); w.writerows(rows)
PYEOF

echo "==== freeze one packing (seed $SPAX_SEED) ===="
if [ ! -f "$ROOT/pack/$ROW.npy" ]; then
    SPAX_SAVE_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
        "$ROOT/row.csv" "$ROOT/seed" > "$ROOT/seed.log" 2>&1
fi
ls "$ROOT/pack"

for lm in $MESHES; do
    tag="m${lm/./}"
    echo "==== L_mesh = $lm ===="
    $PY - "$ROOT/row.csv" "$ROOT/$tag.csv" "$lm" <<'PYEOF'
import csv, sys
src, dst, lm = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(src)))
fn = list(rows[0].keys())
rows[0]['L_mesh'] = lm
w = csv.DictWriter(open(dst, 'w', newline=''), fieldnames=fn)
w.writeheader(); w.writerows(rows)
PYEOF
    if [ ! -f "$ROOT/$tag/Job-$ROW-utx.inp" ]; then
        SPAX_LOAD_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
            "$ROOT/$tag.csv" "$ROOT/$tag" > "$ROOT/$tag.gen.log" 2>&1 \
            || { echo "  generation FAILED"; continue; }
    fi
    grep -m1 "Done: " "$ROOT/$tag.gen.log" | sed 's/^ */  /'
    python3 SpaX_CalculiX.py convert "$ROOT/$tag" > /dev/null
    python3 SpaX_CalculiX.py solve "$ROOT/$tag" --cpus "$CPUS" --jobs 1 | sed 's/^/  /'
    python3 SpaX_PostProcess.py "$ROOT/$tag.csv" "$ROOT/$tag" "$ROOT/$tag.csv.out" \
        > "$ROOT/$tag.post.log" 2>&1
done

echo
echo "==== mesh convergence on one frozen packing, $ROW ===="
python3 - "$ROOT" <<'PYEOF'
import csv, glob, os, re, sys
root = sys.argv[1]
ABQ = 7638133908.0        # $VALDIR/results/results_gas.csv, GAS_v10, Abaqus
print('%-8s %10s %12s %10s %9s %9s' % (
    'L_mesh', 'elements', 'E_x', 'phi_soft', 'vs Abaqus', 'gap'))
rows = []
for f in sorted(glob.glob(os.path.join(root, 'm*.csv.out'))):
    tag = os.path.basename(f).split('.')[0]
    gen = os.path.join(root, tag + '.gen.log')
    nel = '?'
    if os.path.isfile(gen):
        m = re.search(r'Done: \S+ \(\d+ nodes, (\d+) elements\)', open(gen).read())
        if m:
            nel = m.group(1)
    for r in csv.DictReader(open(f)):
        try:
            E = float(r['E_x']); phi = float(r['phi_soft_total'])
        except (KeyError, ValueError):
            continue
        lm = float(r['L_mesh'])
        rows.append((lm, nel, E, phi, r.get('equilibrium_gap', '')))
for lm, nel, E, phi, gap in sorted(rows, reverse=True):
    try:
        g = '%.1e' % float(gap)
    except (TypeError, ValueError):
        g = '-'
    print('%-8.3f %10s %12.5g %10.4f %+8.2f%% %9s' % (
        lm, nel, E, phi, 100 * (E - ABQ) / ABQ, g))
print('\nAbaqus GAS_v10: E_x = %.5g at ~20700 matrix elements' % ABQ)
PYEOF
