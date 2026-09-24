#!/bin/bash -l
# Compare CalculiX against the stored Abaqus results on a real campaign deck.
#
# validation_ccx/params/rve_gas.csv is the campaign's dominant configuration -- L=0.50 at
# L_mesh=0.033, the size shared by 859 of the deck rows -- and its Abaqus
# results are in validation_ccx/results/results_gas.csv: E_x, nu_x, E_z, nu_z, and the
# ACHIEVED porosity and inclusion fraction.
#
# Unlike the homogeneous cube, this cell has a randomly packed microstructure,
# so a fresh generation does not reproduce the campaign's exact geometry: the
# comparison bundles packing, meshing and solver. The achieved porosity is the
# control that says how much. It is measured from the mesh rather than copied
# from the deck, so if it lands on the stored value the geometry is comparable
# and the modulus difference is about the solver; if it does not, the modulus
# difference is about the packing and says nothing about CalculiX.
#
# Runs the shipped default solver -- no SPAX_CCX_SOLVER -- so this is also the
# timing and memory record at campaign size.
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
ROOT=${ROOT:-out_gasccx}
CSV=${CSV:-$VALDIR/params/rve_gas.csv}
CPUS=${CPUS:-8}
JOBS=${JOBS:-2}

export SPAX_CCX="$CCX"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12
# Order 1 is what the campaign used for first-order decks -- see
# validate_gas_order1.sh. Overridable so both can be run and compared.
export SPAX_MESH_ORDER=${SPAX_MESH_ORDER:-2}
# The campaign's own generation seed is not recorded in the tree, so this run
# packs its own microstructures. Fixing the seed at least makes THIS comparison
# reproducible.
export SPAX_SEED=${SPAX_SEED:-20260821}
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT"

echo "==== [1/4] generate (SPAX_SEED=$SPAX_SEED) ===="
if [ ! -f "$ROOT/Job-GAS_v10-utx.inp" ]; then
    SPAX_RESUME=1 "$PY" -u SpaX_Standalone.py "$CSV" "$ROOT" > "$ROOT/gen.log" 2>&1
fi
grep -c 'Written:' "$ROOT/gen.log" | sed 's/^/  decks: /'
grep 'Done: ' "$ROOT/gen.log" | sed 's/^/  /'

echo "==== [2/4] translate (shipped default solver) ===="
python3 SpaX_CalculiX.py convert "$ROOT"
grep -m1 '^\*STATIC' "$ROOT"/Job-GAS_v00-utx-ccx.inp | sed 's/^/  deck says: /'

echo "==== [3/4] solve ===="
/usr/bin/time -v -o "$ROOT/solve.time" \
    python3 SpaX_CalculiX.py solve "$ROOT" --cpus "$CPUS" --jobs "$JOBS"
sed -n 's/.*Maximum resident set size (kbytes): /  peak RSS over all jobs: /p' "$ROOT/solve.time" \
    | $PY -c "
import sys
for l in sys.stdin:
    kb=float(l.split(':')[1]); print('  peak RSS (largest concurrent job): %.0f MB'%(kb/1024))
"
grep -h '# of iterations' "$ROOT"/*.log 2>/dev/null | sort -u | sed 's/^/  /' | head -4

echo "==== [4/4] post-process and compare ===="
python3 SpaX_PostProcess.py "$CSV" "$ROOT" "$ROOT/results_gas_ccx.csv" | tail -12
echo
python3 $VALDIR/compare_ccx.py $VALDIR/results/results_gas.csv "$ROOT/results_gas_ccx.csv" \
    porosity phi_inclusion E_x nu_x E_z nu_z equilibrium_gap
