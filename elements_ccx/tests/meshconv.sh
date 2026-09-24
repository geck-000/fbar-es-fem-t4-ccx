#!/bin/bash -l
# Same-geometry mesh sweep for the F-barES-FEM-T4 arm.
#
# The Abaqus reference is NOT a convergence study: every one of its 16 rows is
# a separate packing (und and drn differ in porosity even within a seed), one
# seed is non-monotone, and neither seed's R increments shrink.  So it cannot
# say whether R converges or wanders.
#
# This can, for OUR element: make_slabconv.py rebuilds the SAME cell at every
# n -- the slab faces and bridge edges are fixed at multiples of 0.1, so every
# phase boundary lands on a mesh plane for any n multiple of 10 -- and the
# drained twin is the SAME MESH with one elastic card rewritten.  R therefore
# carries NO packing noise at all -- the only thing that changes between
# points is h.
#
#   NS="10 15 20 30 40" KG=500 elements_ccx/tests/meshconv.sh
#
# Per mesh: plain C3D4 drained and undrained, and F-barES-FEM-T4 c=1 undrained.
# R = C1111(und)/C1111(drn), denominator always plain C3D4 -- the same
# convention report_abaqus_ratio.py uses, so the arms differ only in the
# element under test.  The _abq.inp decks are for Abaqus C3D4H on Roihu and are
# written but not solved here.
set -eu
cd "$(dirname "$0")/../.."
PY=${PY:-python3}
NS=${NS:-"10 20 30 40"}
KG=${KG:-500}
JIT=${JIT:-0.3}
BRIDGE=${BRIDGE:-one}
LOAD=${LOAD:-x}
CONF=${CONF:-sym}
ROOT=${ROOT:-out_meshconv/kg${KG}_${BRIDGE}_${LOAD}$CONF}
mkdir -p "$ROOT"

for n in $NS; do
  for st in und drn; do
    d="$ROOT/n$n/$st"; mkdir -p "$d"
    $PY elements_ccx/tests/make_slabconv.py "$d/m" "$n" "$st" "$KG" "$JIT" "$BRIDGE" "$LOAD" "$CONF" \
        > "$d/gen.log" 2>&1
  done
  # F-barES-FEM-T4 c=1 on the undrained cell only
  w="$ROOT/n$n/und_fbar1"; mkdir -p "$w"
  $PY elements_ccx/fbares.py "$ROOT/n$n/und/m_ccx.inp" "$w/m_ccx.inp" \
      --elset Sphere_Only --cycles 1 > "$w/gen.log" 2>&1

  for arm in und/m_ccx drn/m_ccx und_fbar1/m_ccx; do
    j="$ROOT/n$n/$arm"
    ( cd "$(dirname "$j")" && CCX_FBAR_C=1 FBAR_CCX_REAL=ccx_fbar \
        FBAR_CCX_MEMMAX=${MEM:-24G} OMP_NUM_THREADS=${NT:-8} \
        ccx_capped "$(basename "$j")" > solve.log 2>&1 ) || true
  done
done

$PY elements_ccx/tests/report_slabconv.py "$ROOT" "$KG" "$LOAD" $NS
