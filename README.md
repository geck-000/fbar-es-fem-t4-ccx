# F-barES-FEM-T4 for CalculiX

A locking-free 4-node tetrahedral element for [CalculiX](http://www.calculix.de/),
implementing the F-barES-FEM-T4 formulation of Onishi, Iida and Amaya
([Int. J. Comput. Methods **15**(7) 1845003, 2018](https://doi.org/10.1142/S0219876218450032)).

The formulation is single-field: no pressure degree of freedom and no
saddle-point solve. Edge-based strain smoothing (ES-FEM) supplies the
deviatoric response, and a cyclic F-bar smoothing of the dilatation supplies
the volumetric one, with the cycle count `c` as the method's free parameter.
This repository delivers `c = 1`, the smallest complete setting. Against Abaqus
`C3D4H` on identical meshes the two formulations share a Richardson-extrapolated
modulus limit to 0.06%, while plain `C3D4` remains 4% stiff at the finest mesh.

Delivery is three `*USER ELEMENT` types over the underlying T4 connectivity:

- `U2` (edge ring): deviatoric edge-smoothed stiffness.
- `U3` (`E A^c` support): volumetric chain `(K V_h) tbar^T sbar`.
- `U4` (4 nodes): null base tetrahedron, geometry only.

## Layout

```text
fbar-es-fem-t4-ccx/
├── elements_ccx/          element sources and deck generator
│   ├── *.f, fbar_lock.c   U2/U3/U4 user elements
│   ├── fbares.py          C3D4 deck -> F-barES deck
│   └── tests/             verification and mesh-convergence checks
└── patches_ccx/           12 patches for a stock CalculiX 2.23 tree
```

## Dependencies

Building the element needs a CalculiX 2.23 source tree with `gcc`, `gfortran`
and `make`. The non-symmetric tangent is factored through PARDISO's
general-matrix path (`mtype = 11`), so CalculiX must be built with its PARDISO
solver (Intel MKL). The Python checks need Python 3 with numpy and scipy,
`verify_nltan.py` compiles a driver with `gfortran`, and the shell drivers need
`bash` and `git` (the patches apply with `git apply`).

## Build and run

Copy the sources from `elements_ccx/` into the CalculiX `src/` directory, add
them to `SCCXF`/`SCCXC` in `Makefile.inc`, apply the patches in `patches_ccx/`
(see its README), and rebuild; patches `0001` to `0005` carry the element
itself. Then `elements_ccx/fbares.py` rewrites a plain `C3D4` deck as
F-barES-FEM-T4, with `elements_ccx/tests/meshconv.py` as the mesh-convergence
driver.
The `*USER ELEMENT` interface caps connectivity at 255 nodes, so `c = 1` is
deliverable and `c >= 2` needs a direct global-assembly pass. Read results from
displacements and reactions, since no element in an F-bar deck carries stress.

## Tests

Run the checks with `python3 elements_ccx/tests/run_tests.py`; adding
`--ccx <binary>` also runs the CalculiX-side checks. `elements_ccx/tests/README.md`
lists what each script covers: the operator and finite-strain checks, the
Fortran chain against `S = E A^c`, the deck generator, and the
mesh-convergence cell.
