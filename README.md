# F-barES-FEM-T4 for CalculiX

A locking-free 4-node tetrahedral element for [CalculiX](http://www.calculix.de/),
implementing the F-barES-FEM-T4 formulation of Onishi, Iida and Amaya
(*Accurate viscoelastic large deformation analysis using F-bar aided edge-based
smoothed finite element method for 4-node tetrahedral meshes (F-BarES-FEM-T4)*,
Int. J. Comput. Methods **15**(7) 1845003, 2018, DOI 10.1142/S0219876218450032).

The element is delivered as two `*USER ELEMENT` types over the underlying T4
connectivity:

| type | role |
|---|---|
| `U2` | deviatoric edge-ring stiffness (symmetric ES-FEM block) |
| `U3` | volumetric chain `E A^c` (the cyclically smoothed F-bar term) |

It is a **displacement formulation** --- no pressure degree of freedom, hence no
saddle-point solve --- that suppresses volumetric locking by combining
edge-based strain smoothing (ES-FEM) with the F-bar dilatation treatment and a
cyclic smoothing of the volumetric strain whose count `c` is the method's free
parameter. On a purpose-built cell that resolves the near-incompressible
inclusion layer, F-barES-FEM-T4 (`c=1`) and Abaqus `C3D4H` share a
Richardson-extrapolated limit to **0.06%**, while plain `C3D4` remains 4% stiff
at the finest mesh.

## Layout

| path | contents |
|---|---|
| `elements_ccx/` | Fortran element sources, the `fbares.py` deck generator, and the verification/validation tests |
| `patches_ccx/` | patches to a stock `ccx_2.23` source tree |
| `validation_ccx/` | validation of the CalculiX backend against stored Abaqus results, and what its missing hybrid elements cost |

## Building CalculiX with the element

1. Copy the `*.f` sources from `elements_ccx/` into the CalculiX `src/`
   directory and add them to `SCCXF` in `Makefile.inc`.
2. Apply the patches in `patches_ccx/` (see its `README.md`).
3. Rebuild.

The F-barES-FEM-T4 element is patches `0001`--`0005` (wide user-element
matrix, element wiring, asymmetric Petrov--Galerkin path, ring-x support
reduction, out-of-core PARDISO).

## Generating a deck and running

`elements_ccx/fbares.py` rewrites a plain `C3D4` deck as F-barES-FEM-T4; see its
docstring and `elements_ccx/tests/slabconv.sh` for the end-to-end driver. Two
hard limits follow from the stencil widths: `*USER ELEMENT` caps connectivity at
255 nodes, so `c=1` is deliverable but `c>=2` is not without a direct
global-assembly pass; and the non-symmetric tangent needs PARDISO's
general-matrix path (`mtype = 11`).

## Verification and validation

The Python checks need Python 3 with numpy and scipy; `verify_nltan.py` also
compiles a driver with gfortran, and the mesh-convergence comparison reads
Abaqus tables through `slabconv_extract.py`.

* `elements_ccx/tests/verify_fbar.py` --- operator checks V1--V4 (unit row sums,
  patch test recovering `C1111 = K + 4G/3`, the `c=0` collapse to selective
  ES-FEM-T4, the volumetric-constraint rank).
* `elements_ccx/tests/verify_fbar_nl.py` --- finite-strain checks N1--N4.
* `elements_ccx/tests/verify_u3_chain.py` --- the Fortran volumetric walk
  against the reference Python operator `S = E A^c`.
* `elements_ccx/tests/stability_modes.py` --- locking and spurious-mode
  behaviour.
* `elements_ccx/tests/make_slabconv.py`, `slabconv.sh`, `report_slabconv.py`,
  `slabconv_extract.py` --- the mesh-convergence cell and the Abaqus `C3D4H`
  comparison.

* `elements_ccx/tests/u3_structure_audit.py` --- the assembled matrix structure
  the volumetric chain actually needs, against what stock `ccx` reserves.

## Validating the backend against Abaqus

`validation_ccx/` is a second, coarser question than the element checks above:
not "is the operator right" but "does solving a real cell through CalculiX give
what Abaqus gave". It covers the plain-`ccx` comparison on campaign decks, the
cost of `ccx` having no hybrid element, mesh convergence on the layered cells,
iterative-solver tolerance, and scaling.

These scripts drive [SpaX](https://github.com/geck-000/SpaX) --- they generate
cells with `SpaX_Standalone.py`, solve them through `SpaX_CalculiX.py`, and
compare against the stored Abaqus tables in `validation_ccx/results/`. SpaX is a
separate repository, so point `SPAX_ROOT` at a checkout of it:

```bash
SPAX_ROOT=/path/to/SpaX  SPAX_CCX=/path/to/ccx  bash validation_ccx/validate_ccx.sh
```

Every script fails immediately with a clear message if `SPAX_ROOT` is unset or
does not look like a SpaX checkout. The decks and Abaqus result tables the
comparisons read travel with the folder, in `validation_ccx/params/` and
`validation_ccx/results/`, so nothing outside the two repositories is needed.
