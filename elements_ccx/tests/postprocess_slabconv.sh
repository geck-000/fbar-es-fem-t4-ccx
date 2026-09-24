#!/bin/bash -l
#SBATCH --job-name=slabconv_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Assemble C1111_eff (C3D4H arm) from every solved slabconv ODB.  The Abaqus
# side of the mesh-convergence study; the CalculiX arm is read locally by
# elements_ccx/tests/report_slabconv.py.  One extractor pass per ODB appends a
# row to OUT; R = C1111(und)/C1111(drn) is then the same quantity the report
# tabulates.  Chains after submit_slabconv_abaqus.sh (--dependency=afterok).
# Env: WORKDIR, ROOT (sweep dir), OUT (csv).
set -e
unset SLURM_GTIDS
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load abaqus/2026

WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
ROOT=${ROOT:-out_slabconv/kg500_one_xsym}
OUT=${OUT:-out_slabconv/slabconv_abq.csv}
cd "$WORKDIR" || exit 1
mkdir -p logs

rm -f "$OUT"
for n in 10 20 30 40 50 60; do
  for st in und drn; do
    d="$ROOT/n$n/$st"
    odb="$d/m_abq.odb"
    sta="$d/m_abq.sta"
    [ -f "$odb" ] || { echo "MISSING n$n/$st: no ODB"; continue; }
    # Gate on the .sta, not the .odb: an OOM-killed solve leaves a truncated
    # ODB that opens but has no completed step (see lamellar_extract.py).
    if ! grep -q "COMPLETED SUCCESSFULLY" "$sta" 2>/dev/null; then
      echo "SKIP n$n/$st: solve did not complete"
      continue
    fi
    abaqus python slabconv_extract.py "$odb" "$n" "$st" "$OUT"
  done
done

echo "===== assembled C3D4H table ====="
cat "$OUT" 2>/dev/null || true
