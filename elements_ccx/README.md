# The F-barES-FEM-T4 element (U2 + U3 + U4)

Fortran sources added to `ccx` to give it the locking-free 4-node tetrahedron
that Abaqus calls `C3D4H`. Three `*USER ELEMENT` types run over the underlying
T4 connectivity:

| Type | Nodes | Role |
|---|---|---|
| `U2` | edge ring | deviatoric edge-smoothed stiffness `V_h B~^T D_dev B~` (ES-FEM) |
| `U3` | `E A^c` support | volumetric chain `(K V_h) tbar^T sbar`, the F-bar term |
| `U4` | 4 | null base tetrahedron, geometry only |

`U3` has two pairings, selected by the deck alone. A letter in `lakon(3:3)` is
the ring size of the Petrov-Galerkin form; the `_` marker written by
`fbares.py --symmetric` selects the Galerkin form `(K V_h) sbar^T sbar`, a Gram
matrix that is positive semidefinite for every `K` and has no coercivity
threshold. `mastruct.c`, `mafillsmas.f` and `resultsmech_u3.f` must all respect
the marker: the support stays full, and the internal force spreads by `sbar`,
not `tbar`, or the symmetric boundary coupling is lost.

The formulation is displacement-based, with no pressure degree of freedom, and
suppresses volumetric locking by combining ES-FEM with the F-bar dilatation
treatment and a cyclic smoothing whose count `c` is the free parameter. The
companion papers carry the formulation and its operator analysis.

| File | Role |
|---|---|
| `u2edge.f` | edge-ring geometry: `V_h` and the smoothed gradients |
| `u3vol.f` | the chain as a row walk; returns `V_h`, `sbar`, `tbar` |
| `u3nl.f` | finite-strain force (NLGEOM only) |
| `u3nltan.f` | consistent tangent of that force (NLGEOM only) |
| `e_c3d_u2.f` | `U2` stiffness |
| `e_c3d_u3.f` | `U3` stiffness, either pairing; raises `nasym` |
| `e_c3d_u4.f` | `U4` null base tet |
| `resultsmech_u2.f` / `resultsmech_u3.f` | internal force through the same calls and pairing; under NLGEOM `U3` returns the finite-strain force and `U2` stays silent |
| `resultsmech_u4.f` | `U4` null base, no internal force |
| `fbar_lock.c` | mutex for the one-time node-to-element map |
| `fbares.py` | deck generator: base tets to `U4`, then `U2`/`U3` |

## Build

1. Copy the `*.f` and `*.c` sources into `<ccx_2.23>/src` and add them to
   `SCCXF`/`SCCXC` in `Makefile.inc`.
2. Apply the patches in `../patches_ccx/` (see its `README.md`).
3. Rebuild.

Two limits follow from the stencil widths: `*USER ELEMENT` caps connectivity at
255 nodes, so `c = 1` is deliverable and `c >= 2` needs a direct
global-assembly pass, and the non-symmetric tangent needs PARDISO's
general-matrix path (`mtype = 11`). Both are in the patch set; the `_` Galerkin
pairing rides the same asymmetric storage path (`nasym` is still raised by
`c > 0`).

## Generate a deck

```text
*USER ELEMENT,TYPE=U2,NODES=<ring size>,INTEGRATIONPOINTS=1,MAXDOF=3
*USER ELEMENT,TYPE=U3,NODES=<support size>,INTEGRATIONPOINTS=1,MAXDOF=3
*USER ELEMENT,TYPE=U4,NODES=4,INTEGRATIONPOINTS=1,MAXDOF=3
```

`fbares.py` writes the whole deck from a `C3D4` deck; see its docstring and
`tests/meshconv.py` for an end-to-end driver. No element in an F-bar deck
carries stress, because the smoothing domains have no shape function and `U4`
is null, so read results from displacements and reactions.

## Verification

`tests/verify_element.py` checks the element against closed-form answers: the
operators and patch test, the `u3vol.f` walk against `S = E A^c`, the
finite-strain force and its consistent tangent, the spurious-mode census, and
the deck generator. `tests/run_tests.py` runs it, plus the Fortran and
CalculiX-side checks; `tests/README.md` lists the entry points.
