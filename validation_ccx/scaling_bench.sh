#!/bin/bash -l
# How far does this ccx build actually go, and does the iterative solver help?
#
# The production cells run to millions of elements. This build is linked
# against SPOOLES alone (check: ldd $(which ccx) shows no PARDISO, no PaStiX),
# a direct solver whose fill-in grows much faster than the model. ccx also
# carries two built-in iterative solvers that need far less memory, but this is
# an ill-conditioned system -- a ~70x phase contrast and a near-incompressible
# phase -- which is where iterative methods are supposed to fail.
#
# So measure rather than assume: the same cell at growing size, solved by both,
# recording wall time, peak resident memory, and E_eff. E_eff is the control --
# an iterative solve that stops at a loose tolerance returns a plausible number
# that is simply wrong, and only comparing it against the direct answer on the
# same deck catches that.
#
# Cell size is grown through L at fixed L_mesh, because that is how the
# campaign scales and because refining L_mesh alone barely moves the element
# count on a cell whose mesh is dominated by inclusion surfaces.
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
ROOT=${ROOT:-out_scale}
SIZES=${SIZES:-"0.30 0.42 0.55"}
TIMEOUT=${TIMEOUT:-3600}
CPUS=${CPUS:-8}
CCX=${CCX:-ccx_spax}
export SPAX_CCX=$CCX
# The stock iterative criterion (c1=0.005) is 0.15% wrong on this problem;
# 1e-5 matches the direct solver to five decimals. See iterative_tolerance.sh.
ITER_TOL=${ITER_TOL:-1e-5}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=3600 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821
export SPAX_MESH_ORDER=2

rm -rf "$ROOT"
mkdir -p "$ROOT"
RESULTS="$ROOT/scaling.tsv"
printf 'L\telements\tnodes\tequations\tsolver\tstatus\twall_s\tpeak_MB\tE_eff\titers\tgap\n' > "$RESULTS"

HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion'

for L in $SIZES; do
    tag="L${L/./}"
    echo "======== L = $L ========"
    # Disp is 1% engineering strain, so it has to track L.
    disp=$($PY -c "print(0.01 * $L)")
    printf '%s\nSC,%s,0.035,Composite,9.43e9,0.33,0.15,0.045,0.008,Uniaxial Tension X,%s,Simple Shear S13,%s,0.0,0.15,2.2e9,0.48,1.0,0.0,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,4.43e7\n' \
        "$HDR" "$L" "$disp" "$disp" > "$ROOT/$tag.csv"

    "$PY" -u SpaX_Standalone.py "$ROOT/$tag.csv" "$ROOT/$tag" \
        > "$ROOT/$tag.gen.log" 2>&1 || { echo "  generation FAILED"; continue; }
    nel=$(sed -n 's/.*Done: SC ([0-9]* nodes, \([0-9]*\) elements).*/\1/p' "$ROOT/$tag.gen.log")
    nno=$(sed -n 's/.*Done: SC (\([0-9]*\) nodes.*/\1/p' "$ROOT/$tag.gen.log")
    echo "  mesh: $nno nodes, $nel elements"

    for solver in SPOOLES "ITERATIVE CHOLESKY"; do
        sname=$(echo "$solver" | tr ' ' '_')
        work="$ROOT/${tag}_${sname}"
        mkdir -p "$work"
        cp "$ROOT/$tag/Job-SC-utx.inp" "$work/"
        SPAX_CCX_SOLVER="$solver" python3 SpaX_CalculiX.py convert \
            "$work/Job-SC-utx.inp" > /dev/null

        job="$work/Job-SC-utx-ccx"
        echo "  -- $solver"
        # ccx writes its results into the working directory, so run it from
        # there on the bare job name rather than on a path.
        set +e
        ( cd "$work" && CCX_ITER_TOL="$ITER_TOL" OMP_NUM_THREADS=$CPUS \
              /usr/bin/time -v -o Job-SC-utx-ccx.time \
              timeout "$TIMEOUT" "$CCX" Job-SC-utx-ccx > Job-SC-utx-ccx.log 2>&1 )
        rc=$?
        set -e

        wall=$(sed -n 's/.*Elapsed (wall clock) time.*: //p' "$job.time" | head -1)
        wall=$($PY -c "
s='''$wall'''.strip().split(':')
print('%.1f' % (float(s[0])*60+float(s[1]) if len(s)==2 else float(s[0])*3600+float(s[1])*60+float(s[2])))
" 2>/dev/null || echo nan)
        peak=$(sed -n 's/.*Maximum resident set size (kbytes): //p' "$job.time" | head -1)
        peak=$($PY -c "print('%.0f' % (${peak:-0}/1024.0))" 2>/dev/null || echo nan)
        neq=$(grep -A1 'number of equations' "$job.log" | tail -1 | tr -d ' ' || echo 0)

        if [ $rc -eq 124 ]; then
            status=TIMEOUT; eeff=nan; egap=nan
        elif grep -q 'Job finished' "$job.log" 2>/dev/null; then
            status=ok
            # tail -1: an under-converged solve prints a warning first.
            read -r eeff egap <<< "$(python3 -c "
import SpaX_CalculiX as c
try:
    r = c.extract_first_order('$job.dat','S11',0.01,$L)
    print('%.8e %.3e' % (r['E_eff'], r.get('equilibrium_gap', float('nan'))))
except Exception:
    print('nan nan')
" 2>/dev/null | tail -1)"
        else
            status=FAILED; eeff=nan; egap=nan
        fi
        iters=$(sed -n 's/^# of iterations = *//p' "$job.log" | tail -1)
        if [ -z "$iters" ]; then iters=-; fi
        echo "     $status  ${wall}s  ${peak}MB  eq=$neq  iters=$iters  E_eff=$eeff  gap=$egap"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$L" "$nel" "$nno" "$neq" "$solver" "$status" "$wall" "$peak" "$eeff" "$iters" "$egap" \
            >> "$RESULTS"
        # The .dat and the factorisation scratch are the bulk of the disk here.
        rm -f "$job.dat" "$job.frd" "$work"/spooles.out "$work"/*.12d
    done
done

echo
echo "==== $RESULTS ===="
column -t -s $'\t' "$RESULTS"
