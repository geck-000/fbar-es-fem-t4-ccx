# F-barES-FEM-T4 tests

| Test | What it verifies |
|---|---|
| `verify_fbar.py` | operator checks V1--V4: unit row sums on `Q,P,E,R` and the chain `S = E A^c`; the patch test recovering `C1111 = K + 4G/3`; the `c = 0` collapse to selective ES-FEM-T4; the volumetric-constraint rank |
| `verify_convergence.py` | manufactured-solution convergence: divergence-free `u*` on the unit cube, single- and two-phase layouts, structured or Delaunay meshes, interior-rate option, `K/G` sweep, relative L2/H1 errors and rates; also the shared manufactured-field and mesh module for the CalculiX checks below |
| `mms_ccx.py` / `mms_ccx.sh` / `mms_report.py` | the same manufactured problem through the real CalculiX element (`fbares.py` conversion, `ccx_fbar`, `.dat` parsing), on structured and Delaunay meshes, with `--symmetric` for the Galerkin deck; the rates are collected in `mms.csv` |
| `verify_pairing.py` | coercivity of the two pairings: the `K/G` threshold from the generalized eigenproblem, the negative-eigenvalue counts, the fluctuation comparison and the error growth across the transition; `--metrics --ccx <binary>` adds the CalculiX cross-check of both pairings on the two-phase cube-inclusion MMS cell |
| `verify_fbar_nl.py` | finite-strain checks N1--N7: `f(0) = 0`, rigid translation, a 20% stretch patch, convergence to the small-strain reduction, frame indifference, the tangent `df/du`, and the rigid-body zero modes |
| `verify_tangent.py` | the consistent tangent of the finite-strain element: `tangent(0)` against the small-strain operator, central differences at 20/50% stretch and shear, and quadratic Newton convergence |
| `verify_nlccx.py` | the finite-strain element through CalculiX `NLGEOM`: reactions on a prescribed boundary against `FbarNL.force(u_ccx)` (1.4e-7, and 3.9e-7 at 1% strain with increments); `--uniaxial` runs the 20 % and 50 % stretch path cutoffs (free residual 4e-5, reactions 6e-5) |
| `verify_nltan.py` | the Fortran consistent tangent `u3nltan` against the prototype: the small-strain operator at `u = 0` (5e-15) and `FbarNL.tangent` at a finite field (6e-15), via a per-edge driver that assembles the global matrix |
| `verify_u3_chain.py` | the Fortran `u3vol.f` walk against the reference Python operator `S = E A^c` |
| `verify_fbares_deck.py` | the `fbares.py` deck generator, walked independently |
| `smoothing_proto.py` | the Python prototype of the smoothing operators the Fortran was checked against |
| `stability_modes.py` | spurious-mode census: distinguishes locking, stable, and unstable elements |
| `make_slabconv.py` / `slabconv.sh` / `report_slabconv.py` / `slabconv_extract.py` | the mesh-convergence cell and the Abaqus `C3D4H` comparison |
| `meshconv.sh` | same-geometry mesh sweep for the F-bar arm |

The two-element patch tests establish *consistency* (a uniform strain state is
reproduced to roundoff). Locking behaviour is covered by the stability census
and the mesh-convergence cell, because an unstable element passes a patch test.
