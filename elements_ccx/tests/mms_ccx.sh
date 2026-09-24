#!/bin/bash -l
# Manufactured-solution convergence through the real CalculiX element.
#
#   NS="8 12 16 24" KGS="10 100 1000 5000" elements_ccx/tests/mms_ccx.sh
#
# Single phase sweeps {c3d4, fbar1} plus fbar0 at the top contrast; two
# phase sweeps {c3d4, fbar1_incl, fbar1_all} (the production configuration
# is fbar1_incl: F-bar in the inclusion only).  Every row appends to
# $OUT/mms.csv, and the report prints the rates and the K/G summary.
set -eu
cd "$(dirname "$0")/../.."
PY=${PY:-python3}
OUT=${OUT:-out_mms}
NS=${NS:-"8 12 16 24"}
KGS=${KGS:-"10 100 1000 5000"}
TWO_KGS=${TWO_KGS:-"10 100 1000 5000"}
mkdir -p "$OUT"
echo "n,layout,kg,arm,h,l2,h1,l2_incl,h1_incl" > "$OUT/mms.csv"

for n in $NS; do
  for kg in $KGS; do
    $PY elements_ccx/tests/mms_ccx.py --out "$OUT" --n "$n" \
        --layout single --kg "$kg" --arms c3d4 fbar1
  done
  $PY elements_ccx/tests/mms_ccx.py --out "$OUT" --n "$n" \
      --layout single --kg 5000 --arms fbar0
  for kg in $TWO_KGS; do
    $PY elements_ccx/tests/mms_ccx.py --out "$OUT" --n "$n" \
        --layout two --kg "$kg" --arms c3d4 fbar1_incl fbar1_all
  done
done

$PY elements_ccx/tests/mms_report.py "$OUT/mms.csv"
