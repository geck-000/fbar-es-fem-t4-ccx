#!/bin/bash -l
# Does the iterative solver reach campaign size?
#
# The decks in params/ are dominated by L=0.50 at L_mesh=0.033 (859 rows), and
# run out to L=1.28. scaling_bench.sh established that the direct solver is not
# a candidate at that size -- 4.4 GB at 488k equations, and no finish on a
# 1.4M-equation cell in an hour on 30 GB -- so this ladder drops it and pushes
# the iterative solver alone, at the campaign's own L_mesh.
#
# Reported per size: equations, wall time, peak resident memory, iteration
# count, and equilibrium_gap. The gap is the accuracy control -- there is no
# direct solve to compare against at these sizes, which is exactly why a check
# that needs no reference solve is the one worth having.
set -eu
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
ROOT=${ROOT:-out_scaleit}
SIZES=${SIZES:-"0.50 0.64 0.80"}
LMESH=${LMESH:-0.033}
TIMEOUT=${TIMEOUT:-7200}
CPUS=${CPUS:-8}
ITER_TOL=${ITER_TOL:-1e-5}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821
export SPAX_MESH_ORDER=2

mkdir -p "$ROOT"
RESULTS="$ROOT/scaling_iterative.tsv"
if [ ! -f "$RESULTS" ]; then
    printf 'L\tL_mesh\telements\tnodes\tequations\tstatus\twall_s\tpeak_MB\titers\tE_eff\tgap\n' > "$RESULTS"
fi

HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion'

for L in $SIZES; do
    tag="L${L/./}"
    echo "======== L = $L, L_mesh = $LMESH ========"
    disp=$($PY -c "print(0.01 * $L)")
    printf '%s\nSC,%s,%s,Composite,9.43e9,0.33,0.15,0.045,0.008,Uniaxial Tension X,%s,Simple Shear S13,%s,0.0,0.15,2.2e9,0.48,1.0,0.0,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,4.43e7\n' \
        "$HDR" "$L" "$LMESH" "$disp" "$disp" > "$ROOT/$tag.csv"

    if [ ! -f "$ROOT/$tag/Job-SC-utx.inp" ]; then
        echo "  meshing..."
        "$PY" -u SpaX_Standalone.py "$ROOT/$tag.csv" "$ROOT/$tag" \
            > "$ROOT/$tag.gen.log" 2>&1 || { echo "  generation FAILED"; continue; }
    fi
    nel=$(sed -n 's/.*Done: SC ([0-9]* nodes, \([0-9]*\) elements).*/\1/p' "$ROOT/$tag.gen.log")
    nno=$(sed -n 's/.*Done: SC (\([0-9]*\) nodes.*/\1/p' "$ROOT/$tag.gen.log")
    echo "  mesh: $nno nodes, $nel elements"

    work="$ROOT/${tag}_it"
    mkdir -p "$work"
    cp "$ROOT/$tag/Job-SC-utx.inp" "$work/"
    # No SPAX_CCX_SOLVER: this exercises the shipped default.
    python3 SpaX_CalculiX.py convert "$work/Job-SC-utx.inp" > /dev/null
    grep -m1 '^\*STATIC' "$work/Job-SC-utx-ccx.inp" | sed 's/^/  deck says: /'

    job="$work/Job-SC-utx-ccx"
    set +e
    ( cd "$work" && CCX_ITER_TOL="$ITER_TOL" OMP_NUM_THREADS=$CPUS \
        /usr/bin/time -v -o Job-SC-utx-ccx.time \
        timeout "$TIMEOUT" "$CCX" Job-SC-utx-ccx > Job-SC-utx-ccx.log 2>&1 )
    rc=$?
    set -e

    wall=$(sed -n 's/.*Elapsed (wall clock) time.*: //p' "$job.time" | head -1)
    wall=$($PY -c "
s='''$wall'''.strip().split(':')
print('%.1f'%(float(s[0])*60+float(s[1]) if len(s)==2 else float(s[0])*3600+float(s[1])*60+float(s[2])))" 2>/dev/null || echo nan)
    peak=$(sed -n 's/.*Maximum resident set size (kbytes): //p' "$job.time" | head -1)
    peak=$($PY -c "print('%.0f'%(${peak:-0}/1024.0))" 2>/dev/null || echo nan)
    neq=$(grep -A1 'number of equations' "$job.log" | tail -1 | tr -d ' ' || echo 0)
    iters=$(sed -n 's/^# of iterations = *//p' "$job.log" | tail -1)
    if [ -z "$iters" ]; then iters=-; fi

    if [ $rc -eq 124 ]; then
        status=TIMEOUT; eeff=nan; gap=nan
    elif grep -q 'Job finished' "$job.log" 2>/dev/null; then
        status=ok
        read -r eeff gap <<< "$(python3 -c "
import SpaX_CalculiX as c
try:
    r = c.extract_first_order('$job.dat','S11',0.01,$L)
    print('%.8e %.3e' % (r['E_eff'], r.get('equilibrium_gap', float('nan'))))
except Exception:
    print('nan nan')
" 2>/dev/null | tail -1)"
    else
        status=FAILED; eeff=nan; gap=nan
    fi

    echo "  $status  ${wall}s  ${peak}MB  eq=$neq  iters=$iters  E_eff=$eeff  gap=$gap"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$L" "$LMESH" "$nel" "$nno" "$neq" "$status" "$wall" "$peak" "$iters" "$eeff" "$gap" \
        >> "$RESULTS"
    rm -f "$job.dat" "$job.frd"
done

echo
column -t -s $'\t' "$RESULTS"
