#!/bin/bash -l
# Reproduce results_homog_qxy.csv with CalculiX instead of Abaqus.
#
# The homogeneous cube is the one cell whose first-order answer is known in
# closed form: no inclusions, so E_eff = E_matrix, nu_eff = nu_matrix and
# G_eff = E/(2(1+nu)) whatever the mesh. Any disagreement between the two
# solvers there is a defect, not discretisation.
#
# Generation flags are copied from hpc/run_channels_q.sh, the script that
# produced the reference table -- SPAX_MESH_ORDER=2 above all, since the
# reference is a quadratic (C3D10) run.
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
OUT=${OUT:-out_ccxval}
CSV=${CSV:-$VALDIR/params/rve_homog_qxy.csv}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_ORDER=2 SPAX_MESH_TIMEOUT=1800 SPAX_OPT_PASSES=2
export SPAX_SLIVER_START=1 SPAX_MAX_RETRIES=12 SPAX_SLIVER_MULT_Q=1.0

echo "==== [1/4] generate decks ===="
rm -rf "$OUT"
mkdir -p "$OUT"
"$PY" -u SpaX_Standalone.py "$CSV" "$OUT" > "$OUT/gen.log" 2>&1
grep -c 'Written:' "$OUT/gen.log" | sed 's/^/  decks written: /'

echo "==== [2/4] translate for CalculiX ===="
python3 SpaX_CalculiX.py convert "$OUT"

echo "==== [3/4] solve with ccx ===="
python3 SpaX_CalculiX.py solve "$OUT" --cpus "${CPUS:-4}" --jobs "${JOBS:-2}"

echo "==== [4/4] post-process ===="
python3 SpaX_PostProcess.py "$CSV" "$OUT" "$OUT/results_homog_qxy_ccx.csv"
