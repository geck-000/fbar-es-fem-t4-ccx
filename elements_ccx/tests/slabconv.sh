#!/bin/bash -l
# Mesh convergence on a cell whose brine slab is RESOLVED.  See make_slabconv.py
# for why the campaign cells cannot answer this (0.35-1.04 elements through the
# slab; no asymptotic regime exists there).
#
#   NS="10 15 20 30 40" KG=500 elements_ccx/tests/slabconv.sh
#
# Per mesh: plain C3D4 drained and undrained, and F-barES-FEM-T4 c=1 undrained.
# R = C1111(und)/C1111(drn), denominator always plain C3D4 -- the same
# convention report_abaqus_ratio.py uses, so the arms differ only in the
# element under test.  The _abq.inp decks are for Abaqus C3D4H on Roihu and are
# written but not solved here.
set -eu
cd "$(dirname "$0")/../.."
PY=/home/giacomo/venvs/sci/bin/python3
NS=${NS:-"10 20 30 40"}
KG=${KG:-500}
JIT=${JIT:-0.3}
BRIDGE=${BRIDGE:-one}
LOAD=${LOAD:-x}
CONF=${CONF:-sym}
ROOT=${ROOT:-out_slabconv/kg${KG}_${BRIDGE}_${LOAD}$CONF}
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
