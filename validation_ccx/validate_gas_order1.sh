#!/bin/bash -l
# The element-order control, and almost certainly the whole story.
#
# validate_gas.sh put CalculiX 0.5-3.6% below the stored Abaqus E_x, with the
# gap growing monotonically with void content. That was a like-for-unlike
# comparison: it ran SPAX_MESH_ORDER=2 (quadratic C3D10) against a reference
# solved with LINEAR tets.
#
# The campaign's own submitter says so. Every first-order campaign in the table
# in hpc/spax_submit.sh carries mesh order 1 (weibull, weibull_layer,
# nlgeom_layer, layercol -- all -utx); only the bending and torsion campaigns
# carry order 2. hpc/generate_array.sh documents the variable as
# "SPAX_MESH_ORDER (2 for bending)", i.e. 1 is the default. rve_gas.csv is a
# first-order deck -- Kappa=0, Mode utx / Mode2 utz -- so it was solved on
# C3D4.
#
# And the size is already measured: hybrid_locking_test.sh found order 1 reads
# E_eff +4.03% stiffer than order 2 on the same frozen geometry. That is the
# direction and the magnitude of the "gap", and its growth with void content is
# what linear-tet over-stiffening does when there is more geometry to resolve.
#
# So rerun the same deck at order 1 and compare like with like.
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

export SPAX_MESH_ORDER=1
ROOT=${ROOT:-out_gasccx_o1} exec bash calculix/validate_gas.sh
