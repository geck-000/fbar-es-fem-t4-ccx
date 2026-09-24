# Patches to CalculiX

Apply in order to a stock `ccx_2.23` tree, after copying the element sources
from `../elements_ccx/` into `src/`:

```bash
cd <ccx_2.23>/src
git apply -p0 <patch>
make -j8
```

| Patch | What it does |
|---|---|
| `0001-wide-user-element-matrix.patch` | widens the user-element assembly path from 60 to 765 DOF so the U3 stencil (173 nodes at `c = 1`) fits |
| `0002-fbar-es-fem-t4-elements.patch` | wires `U2`/`U3`/`U4` into the element and results dispatchers, the digit list and the build |
| `0003-asymmetric-user-element-path.patch` | opens the asymmetric assembly and solve path to user elements (`U3` raises `nasym`); stops truncating the node count at 127 |
| `0004-ring-support-reduction.patch` | stops allocating the identically zero `(support − ring) × (support − ring)` block of `U3` |
| `0005-pardiso-out-of-core.patch` | opt-in out-of-core PARDISO for factors that do not fit in RAM |
| `0006-dedup-matrix-structure.patch` | deduplicates the matrix structure as it is built, so peak memory tracks the final `irow`/`jq` instead of the transient `mast1`/`next` list; the way past the 2^31 insertion wall at the finest meshes |
| `0007-pardiso-memreport-and-single.patch` | `CCX_PARDISO_MEMREPORT` reports what a factorization needs before it starts; `CCX_PARDISO_SINGLE` factors in single precision (−39% peak, −37% time, 1.3 ppm on the reaction) |
| `0008-alloc-trace-and-inpc.patch` | `CCX_ALLOC_TRACE=<MB>` names the allocation about to be attempted; `inpc` sized to the characters actually read (3.110 to 1.592 GB on the finest cell) |
| `0009-symmetric-pairing.patch` | the `_` marker in `lakon(3:3)` selects the Galerkin pairing `K_vol = (K V_h) sbar^T sbar` from the deck (`fbares.py --symmetric`); `mastruct.c`/`mafillsmas.f` skip the ring decode. The element sources ship with the matching `sbar` force spread |
| `0010-finite-strain-force.patch` | stage 1 of the NLGEOM port: builds `u3nl.f` and passes `iperturb` through the U2/U3 force dispatchers; under NLGEOM `U3` returns the finite-strain force and `U2` stays silent |
| `0011-finite-strain-tangent.patch` | stage 2: builds `u3nltan.f`, the consistent tangent of that force, verified against the prototype to 6e-15 |
| `0012-opt-in-line-search.patch` | opens the secant line search to any NLGEOM run through `CCX_LINESEARCH=1`, factor bracketed to [0.05, 1.0]; a robustness aid, not required for the 20 to 50% uniaxial stretch |

The element sources live in `../elements_ccx/`; add them to `SCCXF`/`SCCXC` in
`Makefile.inc`. See `../elements_ccx/README.md`.
