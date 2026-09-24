# Patches to CalculiX

Applied against a stock `ccx_2.23` source tree, in order:

```bash
cd <ccx_2.23>
git apply -p0 <patch>
cd src && make -j8
```

| Patch | What it does |
|---|---|
| `0001-wide-user-element-matrix.patch` | widens the user-element assembly path from 60 to 765 DOF so the U3 volumetric stencil (up to 173 nodes at `c = 1`) fits |
| `0002-fbar-es-fem-t4-elements.patch` | wires `U2`/`U3`/`U4` into the element and results dispatchers and the build; adds the user-element digit list and the EVOL branch |
| `0003-asymmetric-user-element-path.patch` | opens ccx's asymmetric assembly/solve path to user elements (`U3` raises `nasym`), and stops truncating the node count at 127 |
| `0004-ring-support-reduction.patch` | stops allocating the identically-zero `(support − ring) × (support − ring)` block of `U3` |
| `0005-pardiso-out-of-core.patch` | opt-in out-of-core PARDISO for factors that do not fit in RAM |
| `0006-dedup-matrix-structure.patch` | deduplicates the matrix structure as it is built, so peak memory tracks the final `irow`/`jq` (~24x smaller) instead of the transient `mast1`/`next` list -- the only way past the 2^31 insertion wall at `L_mesh` 0.0080 and 0.0060 |
| `0007-pardiso-memreport-and-single.patch` | `CCX_PARDISO_MEMREPORT` reports what a factorisation needs from the ANALYSIS phase alone, in seconds instead of hours; `CCX_PARDISO_SINGLE` factors in single precision, -39% peak and -37% time for 1.3 ppm on the reaction |
| `0008-alloc-trace-and-inpc.patch` | `CCX_ALLOC_TRACE=<MB>` names the allocation about to be attempted, so a kernel OOM stops being silent; and `inpc` is sized to the characters actually read instead of 132 bytes per line -- 3.110 GB to 1.592 GB on the `L_mesh` 0.0060 cell, held for the whole run |
| `0009-symmetric-pairing.patch` | the `_` marker in `lakon(3:3)` selects the Galerkin pairing `K_vol = (K V_h) sbar^T sbar` from the deck (`fbares.py --symmetric`); `mastruct.c`/`mafillsmas.f` skip the ring decode for it. The element sources ship with the matching `sbar` force spread, so stiffness and force use the same pairing |
| `0010-finite-strain-force.patch` | stage 1 of the NLGEOM port: builds `u3nl.f` and passes `iperturb` to the U2/U3 force dispatchers, so under NLGEOM the U3 element returns the finite-strain force of eqs. (1)-(17) and U2 returns nothing (its small-strain force would double-count the deviatoric half) |
| `0011-finite-strain-tangent.patch` | stage 2: builds `u3nltan.f`, passes `iperturb`/`vold` to the U3 element and skips U2 in `mafillsm`, so under NLGEOM the Newton matrix is the exact consistent tangent of the finite-strain force (verified against the prototype to 6e-15) |
| `0012-opt-in-line-search.patch` | opens ccx's existing secant line search to any NLGEOM run via `CCX_LINESEARCH=1`, with the factor bracketed to [0.05, 1.0]; opt-in robustness aid, observed active but not required for the 20–50 % uniaxial stretch and not sufficient alone for the multi-axial MMS field |

The element sources themselves (`u2edge.f`, `u3vol.f`, `e_c3d_u2.f`,
`e_c3d_u3.f`, `e_c3d_u4.f`, `resultsmech_u2.f`, `resultsmech_u3.f`,
`resultsmech_u4.f`, `fbar_lock.c`) live in `../elements_ccx/` and are copied
into `src/` before building; add them to `SCCXF`/`SCCXC` in `Makefile.inc`.
See `../elements_ccx/README.md`.
