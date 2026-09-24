#!/bin/bash -l
#SBATCH --job-name=slabconv
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=0-5
#
# Abaqus C3D4H on the mesh-convergence cell, CSC Roihu.
#
# The companion to elements_ccx/tests/slabconv.sh, which runs the SAME decks
# through CalculiX (plain C3D4 and F-barES-FEM-T4 c=1).  make_slabconv.py emits
# both variants from one geometry -- m_abq.inp with C3D4H in the inclusion,
# m_ccx.inp with C3D4 -- so the only difference between what the two codes see
# is the element keyword.  That is the whole point: it makes C3D4H a fourth arm
# on an identical mesh rather than a reference from another campaign with
# another packing.
#
# Decks live under WORKDIR (default .../test_rve/out_slabconv), one sweep
# directory per CASE, as n*/<und|drn>/m_abq.inp.  Stage them from the laptop:
#   rsync -a out_slabconv/ roihu:/scratch/project_XXXXXX/test_rve/out_slabconv/
#
# then, from the SHARED test_rve dir (so the array logs land in logs/ alongside
# the other campaigns):
#   cd /scratch/project_XXXXXX/test_rve
#   sbatch --array=0-5 submit_slabconv_abaqus.sh
#
# WORKDIR   where the n*/ directories live (default .../test_rve/out_slabconv)
# CASE      which sweep directory under WORKDIR (default kg500_one_xsym)
set -e

# CSC's /etc/profile.d/zz-csc-env.sh bails out for non-interactive shells
# (leaving `module` as the Tcl env-modules that cannot parse Abaqus's .lua Lmod
# file).  The same non-interactive init csc_solve_array.sh uses makes
# `module load abaqus/2026` work however the job was submitted.
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load abaqus/2026

WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve/out_slabconv}
CASE=${CASE:-kg500_one_xsym}
NS=(10 20 30 40 50 60)
n=${NS[${SLURM_ARRAY_TASK_ID:-0}]}

# Abaqus scratch on node-local NVMe when Slurm grants it, else a WORKDIR dir.
SCR="${LOCAL_SCRATCH:-$WORKDIR/scratch_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$SCR"

for st in und drn; do
    d="$WORKDIR/$CASE/n$n/$st"
    [ -f "$d/m_abq.inp" ] || { echo "missing $d/m_abq.inp"; exit 1; }
    cd "$d"
    # scratch= keeps the .odb and temporaries local; interactive so the job
    # does not return before the solve finishes and the array element exits
    # while Abaqus is still writing.
    abaqus job=m_abq input=m_abq.inp \
           cpus="${SLURM_CPUS_PER_TASK:-16}" \
           scratch="$SCR" \
           interactive
    echo "== $CASE n=$n $st done: $(ls -la m_abq.dat 2>/dev/null | awk '{print $5}') bytes"
done

[ "$SCR" != "$LOCAL_SCRATCH" ] && rm -rf "$SCR"
