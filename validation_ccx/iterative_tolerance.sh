#!/bin/bash -l
# Does tightening the conjugate-gradient tolerance recover the direct-solver
# answer, and what does it cost?
#
# ccx's iterative solvers are the only ones whose memory fits a large periodic
# cell, but their tolerance is a compile-time constant (eps=1.e-4 in preiter.c)
# that no input deck can reach. calculix/patches/0001-iterative-tolerance.patch
# exposes it as CCX_ITER_TOL; this sweeps it against the SPOOLES answer on the
# SAME deck, so the only thing changing is the linear solve.
#
# Two things are measured per tolerance: the error in E_eff against the direct
# solve, and the equilibrium gap -- the disagreement between the
# volume-averaged stress and the reference-point reaction. The second matters
# because it needs no reference solve, so it is the check that remains
# available on a production cell too big to solve directly even once.
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

CCX=${CCX:-ccx_spax}
ROOT=${ROOT:-out_itertol}
DECK=${DECK:-out_scale/L030/Job-SC-utx.inp}
L=${L:-0.30}
CPUS=${CPUS:-8}
EPSLIST=${EPSLIST:-"5e-3 1e-3 1e-4 1e-5 1e-6"}

if [ ! -f "$DECK" ]; then
    echo "need $DECK -- run calculix/scaling_bench.sh with SIZES=\"0.30\" first" >&2
    exit 1
fi
if ! command -v "$CCX" > /dev/null; then
    echo "$CCX not found -- apply calculix/patches/0001-iterative-tolerance.patch" >&2
    exit 1
fi

rm -rf "$ROOT"; mkdir -p "$ROOT"
TSV="$ROOT/tolerance.tsv"
printf 'solver\teps\tstatus\titerations\twall_s\tpeak_MB\tE_eff\tE_err_pct\tequilibrium_gap\n' > "$TSV"

extract () {   # dat  -> "E_eff gap"
    python3 -c "
import SpaX_CalculiX as c
try:
    r = c.extract_first_order('$1', 'S11', 0.01, $L)
    print('%.10e %.3e' % (r['E_eff'], r.get('equilibrium_gap', float('nan'))))
except Exception:
    print('nan nan')
" 2>/dev/null | tail -1
}

run () {       # name  solver  eps
    local name=$1 solver=$2 eps=$3
    local work="$ROOT/$name"
    mkdir -p "$work"
    cp "$DECK" "$work/"
    SPAX_CCX_SOLVER="$solver" python3 SpaX_CalculiX.py convert \
        "$work/$(basename "$DECK")" > /dev/null

    set +e
    ( cd "$work" && CCX_ITER_TOL="$eps" OMP_NUM_THREADS=$CPUS \
        /usr/bin/time -v -o t.time "$CCX" Job-SC-utx-ccx > run.log 2>&1 )
    local rc=$?
    set -e

    local wall peak iters status
    wall=$(sed -n 's/.*Elapsed (wall clock) time.*: //p' "$work/t.time" | head -1)
    wall=$(python3 -c "
s='''$wall'''.strip().split(':')
print('%.1f'%(float(s[0])*60+float(s[1]) if len(s)==2 else float(s[0])*3600+float(s[1])*60+float(s[2])))" 2>/dev/null || echo nan)
    peak=$(sed -n 's/.*Maximum resident set size (kbytes): //p' "$work/t.time" | head -1)
    peak=$(python3 -c "print('%.0f'%(${peak:-0}/1024.0))" 2>/dev/null || echo nan)
    iters=$(sed -n 's/^# of iterations = *//p' "$work/run.log" | tail -1)
    [ -z "$iters" ] && iters=-

    if [ $rc -eq 0 ] && grep -q 'Job finished' "$work/run.log"; then
        status=ok
        read -r eeff gap <<< "$(extract "$work/Job-SC-utx-ccx.dat")"
    else
        status=FAILED; eeff=nan; gap=nan
    fi

    local err=nan
    if [ -n "${REF:-}" ] && [ "$eeff" != nan ]; then
        err=$(python3 -c "print('%+.4f'%(100*($eeff-$REF)/$REF))")
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$solver" "$eps" "$status" "$iters" "$wall" "$peak" "$eeff" "$err" "$gap" >> "$TSV"
    echo "  $solver eps=$eps -> $status  ${wall}s ${peak}MB iters=$iters E=$eeff err=$err% gap=$gap"
    rm -f "$work"/*.dat "$work"/*.frd "$work"/*.12d "$work"/spooles.out
    LAST_EEFF=$eeff
}

echo "==== reference: SPOOLES (direct) ===="
run ref SPOOLES 1e-4
REF=$LAST_EEFF
echo "  reference E_eff = $REF"

echo "==== ITERATIVE CHOLESKY, tolerance sweep ===="
for eps in $EPSLIST; do
    run "chol_$eps" "ITERATIVE CHOLESKY" "$eps"
done

echo
column -t -s $'\t' "$TSV"
