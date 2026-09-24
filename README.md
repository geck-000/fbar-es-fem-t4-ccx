# F-barES-FEM-T4 for CalculiX

A locking-free 4-node tetrahedral element for [CalculiX](http://www.calculix.de/),
implementing the F-barES-FEM-T4 formulation of Onishi, Iida and Amaya
(*Accurate viscoelastic large deformation analysis using F-bar aided edge-based
smoothed finite element method for 4-node tetrahedral meshes (F-BarES-FEM-T4)*,
Int. J. Comput. Methods **15**(7) 1845003, 2018,
DOI [10.1142/S0219876218450032](https://doi.org/10.1142/S0219876218450032)).

The element is a displacement formulation with no pressure degree of freedom and
no saddle-point solve. It suppresses volumetric locking by combining edge-based
strain smoothing (ES-FEM) with the F-bar treatment of the dilatation. The cycle
count `c` of the cyclic smoothing is the free parameter of the method, and
`c = 1` is the setting this repository delivers. On a cell that resolves the
near-incompressible inclusion layer, F-barES-FEM-T4 (`c = 1`) and Abaqus
`C3D4H` share a Richardson-extrapolated limit to 0.06%, while plain `C3D4`
remains 4% stiff at the finest mesh.

Three `*USER ELEMENT` types run over the underlying T4 connectivity:

| type | nodes | role |
|---|---|---|
| `U2` | edge ring | deviatoric edge-smoothed stiffness (symmetric ES-FEM block) |
| `U3` | `E A^c` support | volumetric chain `(K V_h) tbar^T sbar` (the F-bar term) |
| `U4` | 4 | null base tetrahedron, geometry only |

## Layout

| path | contents |
|---|---|
| `elements_ccx/` | Fortran element sources, the `fbares.py` deck generator, and the tests |
| `patches_ccx/` | patches to a stock `ccx_2.23` source tree |
| `validation_ccx/` | validation of the CalculiX backend against stored Abaqus results |

## Build

1. Copy the `*.f` and `*.c` sources from `elements_ccx/` into the CalculiX
   `src/` directory and add them to `SCCXF`/`SCCXC` in `Makefile.inc`.
2. Apply the patches in `patches_ccx/` (see its README).
3. Rebuild.

Patches `0001` to `0005` carry the element: the widened user-element matrix, the
element wiring, the asymmetric Petrov-Galerkin path, the ring-by-support
reduction, and the out-of-core PARDISO option.

## Generating a deck

`elements_ccx/fbares.py` rewrites a plain `C3D4` deck as F-barES-FEM-T4; see its
docstring and `elements_ccx/tests/slabconv.sh` for an end-to-end driver. Two
limits apply: `*USER ELEMENT` caps connectivity at 255 nodes, so `c = 1` is
deliverable while `c >= 2` needs a direct global-assembly pass, and the
non-symmetric tangent requires PARDISO's general-matrix path (`mtype = 11`).

## Verification

The Python checks need Python 3 with numpy and scipy, and `verify_nltan.py`
compiles a driver with gfortran. The tests verify the element against
closed-form answers, the compiled Fortran against the Python operator, and the
deck generator against the original deck. The main entry points are
`verify_fbar.py` (operator checks), `verify_fbar_nl.py` (finite strain),
`verify_u3_chain.py` (the volumetric walk against `S = E A^c`),
`verify_fbares_deck.py` (generator), `stability_modes.py` (locking and spurious
modes), and `u3_structure_audit.py` (assembled matrix structure).
`make_slabconv.py`, `slabconv.sh` and `report_slabconv.py` run the
mesh-convergence cell against the Abaqus `C3D4H` tables.

## Validation

`validation_ccx/` asks a separate question: does a real cell solved through
CalculiX match what Abaqus gave? The scripts drive
[SpaX](https://github.com/geck-000/SpaX) to generate cells and solve them, and
they compare the results with the stored Abaqus tables in
`validation_ccx/results/`. Point `SPAX_ROOT` at a SpaX checkout:

```bash
SPAX_ROOT=/path/to/SpaX  SPAX_CCX=/path/to/ccx  bash validation_ccx/validate_ccx.sh
```

Every script fails with a clear message when `SPAX_ROOT` is unset or does not
look like a SpaX checkout. The decks and the reference tables travel with the
folder, so nothing outside the two repositories is needed.
