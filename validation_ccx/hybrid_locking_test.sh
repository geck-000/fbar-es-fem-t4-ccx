#!/bin/bash -l
# What does CalculiX's lack of hybrid elements cost?
#
# The generator gives a phase hybrid elements when its nu >= SPAX_HYBRID_NU
# (0.45 by default). In these decks that is the brine: K=2.2 GPa, G=44.3 MPa,
# so nu = 0.490. CalculiX has no hybrid element -- it rejects C3D4H outright --
# so the translated deck runs that phase on the plain displacement element,
# which is the formulation that volumetrically locks as nu -> 0.5.
#
# Locking cannot be measured by comparing against Abaqus here (no licence), and
# it cannot be measured on a homogeneous cube either: uniform strain is
# represented exactly by any conforming element, so a single-phase cell shows
# nothing. It needs an inhomogeneous deformation and a control. This builds
# both, on ONE frozen packing so that geometry is held exactly fixed:
#
#   brine  nu=0.490  order 1 (C3D4)  vs order 2 (C3D10)   -> gap A
#   twin   nu=0.300  order 1 (C3D4)  vs order 2 (C3D10)   -> gap B
#
# The twin has the SAME inclusion Young's modulus as the brine
# (E = 9KG/(3K+G) = 1.320e8) and differs only in compressibility, so anything
# in A beyond B is the near-incompressibility, not the stiffness contrast or
# the mesh. A third run refines the brine order-2 mesh: if E_eff does not move,
# C3D10 is converged and the missing hybrid formulation costs nothing
# measurable at this nu.
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
ROOT=${ROOT:-out_lock}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=1800 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821

rm -rf "$ROOT"
mkdir -p "$ROOT/pack"

HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion'
# Same run_id in both rows: SPAX_LOAD_PACKING keys the frozen array on it, which
# is what makes the two materials share one geometry exactly.
BRINE='LOCK,0.30,0.050,Composite,9.43e9,0.33,0.15,0.045,0.008,Uniaxial Tension X,0.003,Simple Shear S13,0.003,0.0,0.15,2.2e9,0.48,1.0,0.0,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,4.43e7'
TWIN='LOCK,0.30,0.050,Composite,9.43e9,0.33,0.15,0.045,0.008,Uniaxial Tension X,0.003,Simple Shear S13,0.003,0.0,0.15,1.3201e8,0.30,1.0,0.0,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Solid,2.2e9,4.43e7'

printf '%s\n%s\n' "$HDR" "$BRINE" > "$ROOT/brine.csv"
printf '%s\n%s\n' "$HDR" "$TWIN"  > "$ROOT/twin.csv"
sed 's/,0.30,0.050,/,0.30,0.035,/' "$ROOT/brine.csv" > "$ROOT/brine_fine.csv"

echo "==== freeze one packing ===="
SPAX_SAVE_PACKING="$ROOT/pack" SPAX_MESH_ORDER=1 \
    "$PY" -u SpaX_Standalone.py "$ROOT/brine.csv" "$ROOT/seed" > "$ROOT/seed.log" 2>&1
ls -la "$ROOT/pack"

run_case () {          # name  csv  order
    local name=$1 csv=$2 order=$3
    echo "==== $name (order $order) ===="
    SPAX_LOAD_PACKING="$ROOT/pack" SPAX_MESH_ORDER="$order" \
        "$PY" -u SpaX_Standalone.py "$csv" "$ROOT/$name" > "$ROOT/$name.gen.log" 2>&1
    grep -m1 'Element types' "$ROOT/$name.gen.log" || true
    grep -m1 'Done: LOCK' "$ROOT/$name.gen.log" || true
    python3 SpaX_CalculiX.py convert "$ROOT/$name" > /dev/null
    python3 SpaX_CalculiX.py solve   "$ROOT/$name" --cpus "${CPUS:-4}" --jobs 2
    python3 SpaX_PostProcess.py "$csv" "$ROOT/$name" "$ROOT/$name.csv" \
        > "$ROOT/$name.post.log" 2>&1
    grep -E '^LOCK' "$ROOT/$name.post.log" || tail -3 "$ROOT/$name.post.log"
}

run_case brine_o1      "$ROOT/brine.csv"      1
run_case brine_o2      "$ROOT/brine.csv"      2
run_case twin_o1       "$ROOT/twin.csv"       1
run_case twin_o2       "$ROOT/twin.csv"       2
run_case brine_o2_fine "$ROOT/brine_fine.csv" 2

echo
echo "==== summary ===="
python3 - "$ROOT" <<'PYEOF'
import csv, os, sys
root = sys.argv[1]
rows = {}
for name in ('brine_o1', 'brine_o2', 'twin_o1', 'twin_o2', 'brine_o2_fine'):
    p = os.path.join(root, name + '.csv')
    if not os.path.isfile(p):
        continue
    with open(p) as f:
        r = list(csv.DictReader(f))
    if r:
        rows[name] = r[0]

def g(name, col):
    try:
        return float(rows[name][col])
    except (KeyError, ValueError, TypeError):
        return float('nan')

print('%-16s %14s %14s %10s' % ('case', 'E_eff', 'G_eff', 'phi_incl'))
for name in ('brine_o1', 'brine_o2', 'brine_o2_fine', 'twin_o1', 'twin_o2'):
    if name in rows:
        print('%-16s %14.6e %14.6e %10.4f' % (
            name, g(name, 'E_eff'), g(name, 'G_eff'), g(name, 'phi_inclusion')))

def gap(a, b, col):
    x, y = g(a, col), g(b, col)
    return 100.0 * (x - y) / y if y else float('nan')

print()
print('order 1 vs order 2, same frozen geometry:')
print('  A  brine (nu=0.490)  E_eff %+7.3f %%   G_eff %+7.3f %%'
      % (gap('brine_o1', 'brine_o2', 'E_eff'), gap('brine_o1', 'brine_o2', 'G_eff')))
print('  B  twin  (nu=0.300)  E_eff %+7.3f %%   G_eff %+7.3f %%'
      % (gap('twin_o1', 'twin_o2', 'E_eff'), gap('twin_o1', 'twin_o2', 'G_eff')))
print()
print('order-2 mesh refinement (0.050 -> 0.035), brine:')
print('     E_eff %+7.3f %%   G_eff %+7.3f %%'
      % (gap('brine_o2', 'brine_o2_fine', 'E_eff'),
         gap('brine_o2', 'brine_o2_fine', 'G_eff')))
PYEOF
