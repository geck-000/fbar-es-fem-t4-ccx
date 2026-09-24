# Solving the SpaX decks with CalculiX

`SpaX_CalculiX.py` (top level) translates the generator's Abaqus `.inp` into a
CalculiX deck, runs `ccx`, and reads the results back. The post-processor then
produces the same `results.csv` it produces from an ODB — with **plain
`python3` and no Abaqus licence anywhere in the chain**.

> **A note on script names in this document.** Some scripts named in the
> narrative below were one-off measurements and were removed once they had
> served their purpose; they are recoverable from the git history. What they
> established is reported here, and the tests that are still run live in
> `elements_ccx/tests/` and `validation_ccx/`.


```bash
python3 SpaX_Standalone.py  params.csv out/      # unchanged
python3 SpaX_CalculiX.py    all        out/      # translate + solve
python3 SpaX_PostProcess.py params.csv out/ results.csv
```

Converted decks are written beside the originals as `Job-<id>-<mode>-ccx.inp`,
so an Abaqus deck and its CalculiX twin coexist in one directory. The
post-processor prefers an `.odb` when both are present, and records which
solver produced each row in a new `solver` column. `SPAX_SOLVER=calculix`
reverses that preference, which is how to read both out of one directory.

## What this folder is

Verification runs, not part of the toolkit.

| Script | Question it answers |
|---|---|
| `validate_ccx.sh` | Does CalculiX reproduce the stored Abaqus results? |
| `validate_gas.sh`, `validate_gas_order1.sh` | The same check on the real campaign deck, at order 2 and order 1 |
| `mesh_convergence.sh` | Is the residual a mesh-density difference? |
| `iterative_tolerance.sh` | Does tightening the CG tolerance recover the direct answer? |
| `scaling_bench.sh`, `scaling_iterative.sh` | How far does this ccx build scale? |
| `hybrid_locking_test.sh` | What does CalculiX's lack of hybrid elements cost? |
| `compare_ccx.py` | Row-by-row relative difference between two results tables. |

Every script resolves the SpaX tree through `SPAX_ROOT` (default
the checkout this script sits in), so this folder is identical whether it is run from the
SpaX checkout or from the standalone F-barES-FEM-T4 element repo that vendors
it.

## Validation: the homogeneous cube

`params/rve_homog_qxy.csv` is the right reference because it has no inclusions:
`E_eff`, `nu_eff` and `G_eff` are then closed-form and mesh-independent, so a
disagreement there is a defect rather than discretisation.

Against the stored `results/results_homog_qxy.csv` (Abaqus, `SPAX_MESH_ORDER=2`):

| Quantity | Agreement |
|---|---|
| `E_eff`, `nu_eff`, `G_eff` | **3e-6 %** — eight significant figures |
| `D_classical` | 3e-6 % (derived from the above) |
| `D_rve`, `E_bending` | CalculiX 1.0–1.4 % higher |
| `l` | 13–49 % higher |

The first-order agreement is conclusive: the periodic-BC translation, the
volume averaging and the modulus fit are all correct.

**The bending gap is not attributable to the solver.** Re-meshing the same cell
at `L_mesh` 0.028 / 0.033 / 0.040 moves `D_rve` by 1.27 % *within CalculiX
alone*, while `E_eff`/`nu_eff`/`G_eff` do not move at all. The gap against a
table built on a mesh this tree cannot reproduce is therefore the same size as
the mesh noise, and the two cannot be separated from the stored table.
Separating them needs one Abaqus run on a deck generated here — the `.inp` the
converter reads is the same file Abaqus would solve, so that comparison is
exact when someone has a licence.

`l` amplifies the gap because it comes from the small difference
`D_rve - D_classical`; for the homogeneous cube the true length scale is zero
and the whole excess (0.5–3.3 % in Abaqus, 1.9–4.3 % in CalculiX) is numerical.
`analyze lengthscale` divides it out through the homogeneous baseline
`f_quad(L)`, so this calibrates away **provided the baseline and the porous
cells come from the same solver**. Do not mix solvers inside one length-scale
analysis.

## Validation: a real campaign deck

`params/rve_gas.csv` is the campaign's dominant configuration — L=0.50 at
L_mesh=0.033, the size shared by 859 of the deck rows — and its Abaqus results
are stored in `results/results_gas.csv`. `validate_gas.sh` reruns it through
CalculiX on the shipped default solver. Six cells, 46k–187k elements, up to
~1.16M equations.

Unlike the homogeneous cube this cell is randomly packed, and the campaign's
generation seed is not recorded in the tree, so a fresh run cannot reproduce
its exact microstructure. The comparison therefore bundles packing, meshing and
solver, and the achieved phase fractions are what say how much of each:

| run | total soft fraction, Abaqus / CalculiX | E_x [GPa] | Δ | E_z [GPa] | Δ | ν_x |
|---|---|---|---|---|---|---|
| GAS_v00 | 0.0202 / 0.0198 | 9.157 / 9.109 | −0.53 % | 9.179 / 9.157 | −0.24 % | 0.3303 / 0.3308 |
| GAS_v02 | 0.0397 / 0.0395 | 8.840 / 8.734 | −1.20 % | 8.890 / 8.845 | −0.51 % | 0.3303 / 0.3310 |
| GAS_v04 | 0.0593 / 0.0592 | 8.533 / 8.378 | −1.82 % | 8.626 / 8.530 | −1.12 % | 0.3304 / 0.3310 |
| GAS_v06 | 0.0793 / 0.0789 | 8.220 / 8.055 | −2.00 % | 8.353 / 8.222 | −1.57 % | 0.3303 / 0.3309 |
| GAS_v08 | 0.0990 / 0.0985 | 7.919 / 7.705 | −2.71 % | 8.065 / 7.950 | −1.42 % | 0.3304 / 0.3309 |
| GAS_v10 | 0.1183 / 0.1185 | 7.638 / 7.366 | −3.57 % | 7.805 / 7.632 | −2.22 % | 0.3304 / 0.3310 |

The cells hold the same amount of each phase — total soft fraction agrees to
0.1–2 % — but not in the same arrangement. Every solve converged
(`equilibrium_gap` 1.9e-7 to 1.7e-6).

ν agrees to 0.2 % and is flat across the series, while E is 0.5–3.6 % low with
the gap **growing monotonically with void content**. A solver discrepancy would
not be porosity-dependent; that shape belongs to geometry or discretisation.

Two things are worth stating about the reference before reading anything into
that. The current `results_gas.csv` is the **2026-08-07 re-run**, made after the
2026-08-06 mesher-geometry corrections (`98226e6`, `4aeb689`, `2f689f8`), so it
is a modern Abaqus reference and not the legacy table — the numbers it replaced
were 1–2 % different and carried no achieved-phase columns at all. But **the
element count behind it is not recorded anywhere in the tree.** (`results_scf.csv`
lists ~20 700 matrix elements for a run named GAS_v00, but that file dates from
2026-07-20, before those same fixes, so it describes a different mesh and cannot
be used to characterise this one.)

### It was quadratic against linear

The table above compares the wrong things. It was generated with
`SPAX_MESH_ORDER=2`, and the reference was solved with **linear** tets.

The campaign's own submitter records this. Every first-order campaign in the
table in `Paper/hpc/submit_campaign.sh` carries mesh order 1 — `weibull`,
`weibull_layer`, `nlgeom_layer`, `layercol`, all `-utx` — and only the bending
and torsion campaigns carry 2; `Paper/hpc/submit_campaign.sh` documents the variable
as "`SPAX_MESH_ORDER` (2 for bending)". `rve_gas.csv` is a first-order deck
(`Kappa=0`, `Mode` utx, `Mode2` utz), so it was solved on C3D4.

And the size was already measured here: `hybrid_locking_test.sh` found order 1
reads **+4.03 %** stiffer than order 2 on one frozen geometry. That is the
direction and magnitude of the "gap", and its growth with void content is what
linear-tet over-stiffening does when there is more geometry to resolve.

Rerunning the same deck at order 1 (`validate_gas_order1.sh`):

| run | E_x: Abaqus | ccx order 2 | ccx order 1 |
|---|---|---|---|
| GAS_v00 | 9.157e9 | −0.53 % | **+0.01 %** |
| GAS_v02 | 8.840e9 | −1.20 % | **−0.06 %** |
| GAS_v04 | 8.533e9 | −1.82 % | **−0.24 %** |
| GAS_v06 | 8.220e9 | −2.00 % | **−0.11 %** |
| GAS_v08 | 7.919e9 | −2.71 % | **+0.05 %** |
| GAS_v10 | 7.638e9 | −3.57 % | **−0.54 %** |

E_z behaves the same way: −0.24…−2.22 % at order 2, +0.44…−0.06 % at order 1.
ν agrees to 0.04 %.

The deficit collapses to **±0.54 %**, and — the part that matters — the
monotonic trend with void content is gone. What remains scatters around zero.

A re-packing study of GAS_v10 at four seeds, everything else fixed, gives a
population s.d. of **0.29 %** and a full spread of **0.77 %**: the residual
above sits inside its own packing noise. The same study separates the order
effect from the scatter — re-packing moves the answer by a third of a percent,
changing the element order by three and a half. The scatter was measured at
order 2 while the matched comparison is at order 1; it is a property of the
geometry rather than the element, so it should carry over.

**Nothing is left over for the solver.** Which is the same answer the
homogeneous cube gave at 6–8 significant figures, now confirmed on a real
microstructure at campaign scale.

The lesson is about the comparison, not the code: match the element order
before reading anything into a modulus difference. A 4 % discrepancy that grows
with porosity looked like physics and was a deck setting.

### The bug this comparison found

`phi_inclusion` initially came back at exactly **2×** the campaign's value while
the porosity printed beside it was right. ccx emits the EVOL block with every
element listed twice under one header; anything reading volumes through the
label→volume dict was immune — the volume-averaged stress, `V_solid`,
`porosity` and hence `E_eff` were all correct — so only the raw sum behind
`phi_inclusion` doubled. Total meshed volume now sums to L³ exactly on a
void-free cell.

The defect was in the reader, not the generator, and nothing about the moduli
would have revealed it. That is the argument for comparing achieved geometry
alongside results rather than results alone.

## Hybrid elements

CalculiX has none, and says so rather than guessing:

```
*ERROR reading *ELEMENT:
C3D4H    is an unknown element type
```

There is no route to one *from the deck*. So the converter strips the `H`, and
the phase the generator marked hybrid — the brine, at `nu >= SPAX_HYBRID_NU` —
runs on the plain displacement element that is supposed to volumetrically lock
as `nu -> 0.5`.

`*USER ELEMENT` is not the escape hatch it is in Abaqus, but it is also not
nothing — see *What `*USER ELEMENT` actually is* below, which corrects an
earlier claim in this file that it was merely the substructure interface.

The CalculiX community thread on this
([discourse 1509](https://calculix.discourse.group/t/incompressible-material-models/1509))
confirms the limitation from the manual — *"Perfectly incompressible materials
require hybrid finite elements. CalculiX does not provide such elements"* — and
offers two remedies. **Neither transfers to this problem**, and it is worth
saying why, because both look applicable at first glance:

* **Reduced integration (`C3D20R`), the manual's advice for isochoric
  behaviour.** CalculiX has reduced integration only for hexahedra — `C3D8R`
  and `C3D20R`. There is no reduced-integration tetrahedron (`C3D4`, `C3D10`,
  `C3D15` only), and the periodic mesher produces tets for these geometries.
  Not available.

* **Accept `nu ≈ 0.475`** (ccx's own default when a hyperelastic `D1` is zero).
  That advice is for `*HYPERELASTIC`, where incompressibility is a modelling
  choice. Here the brine is `*ELASTIC` and its `nu` is not a knob — it encodes
  a measured bulk modulus. Because `G` is only 0.44 MPa, `K` is violently
  sensitive to `nu` near 0.5:

  | ν | K | vs the specified 2.2 GPa |
  |---|---|---|
  | 0.49990 | 2.2 GPa | as specified |
  | 0.499 | 220 MPa | 0.10× |
  | 0.490 | 21.9 MPa | 0.010× |
  | 0.475 | 8.65 MPa | **0.0039×** |

  Capping at 0.475 cuts the brine's bulk modulus by **254×**, turning the
  undrained cell into something close to its drained twin. Since K = 2.2 GPa
  against the ice's 9.25 GPa is exactly what makes the confined brine a load
  path, that does not work around the numerics — it deletes the physics under
  test.

### Fluid elements (`F3D*`) are not a third remedy

The brine *is* a fluid, so mapping it to CalculiX's fluid elements looks like
the natural way out of needing a hybrid solid element. It is not, and the
reason is worth recording so it is not re-proposed.

**There is no `F3D10`.** ccx 2.23 recognises exactly four fluid element names,
in `src/elements.f`:

```fortran
!     3D fluid element
      elseif((label.eq.'F3D8    ').or.
     &       (label.eq.'F3D8R   ').or.
     &       (label.eq.'F3D4    ').or.
     &       (label.eq.'F3D6    ')) then
```

Linear hex, reduced hex, linear tet, linear wedge. Anything else — `F3D10`
included — falls through to `*ERROR reading *ELEMENT: ... is an unknown element
type`, the same message `C3D4H` gets. There is no quadratic fluid element
because the CFD discretisation is not built on one.

**And they are not a material model — they are a different analysis.** `F3D*`
elements are the mesh for ccx's CFD solver, reached through a `*CFD` step.
`src/cfds.f` sets `nmethod=4` and comments *"ONLY CFD-CALCULATIONS WITH THE
CBS-METHOD ARE ALLOWED"*: a transient, characteristic-based-split
Navier–Stokes solve for velocity, pressure and temperature fields of a
**flowing** fluid, with its own boundary conditions, its own time integration
and its own output. It is not a constitutive law that a `*STATIC` step can
assign to part of a solid mesh, and there is no monolithic solid/fluid
stiffness assembly in ccx for a homogenisation to ride on.

What the RVE needs is the opposite of what CFD computes. The brine is trapped,
not flowing; its velocity field is not the unknown; the wanted answer is a
stiffness contribution to a static periodic cell. Running the brine through a
CBS Navier–Stokes solve would produce a transient flow field and no term in the
effective elasticity tensor.

**Acoustic elements do not exist either.** A pressure-only acoustic element is
the other standard way to represent a trapped inviscid fluid, and ccx has none
— `grep -ril acoustic src/*.c src/*.f` returns nothing in 2.23.

**The idea does have a correct cousin, and it is Abaqus's `F3D3`/`F3D4`.**
Those names collide confusingly with ccx's CFD elements but are a completely
different thing: *hydrostatic fluid (cavity)* elements — **surface** elements
lining a closed cavity, carrying a single cavity-pressure degree of freedom at
a reference node and a pressure–volume law `p = -K (V - V₀)/V₀`. That *is* the
right way to say "the brine is a fluid" inside a stress analysis, and it
sidesteps volumetric locking completely, because the incompressibility stops
living in a displacement element at all.

It is also emulable in ccx **without touching Fortran**, for the linear
kinematics these decks use. Under small strain the cavity volume change is a
*linear* form in the boundary displacements, `ΔV = Σᵢ Aᵢ·uᵢ` with `Aᵢ` the
nodal area vectors of the brine surface, so:

* delete the brine solid elements, leaving a cavity;
* add one node `R`, whose DOF 1 carries `ΔV`;
* one `*EQUATION` tying `ΔV` to the surface displacements — the converter
  already writes `*EQUATION` for the periodic BCs;
* one grounded `SPRING1` on that DOF with stiffness `K/V₀`, which is the
  fluid's energy `½ (K/V₀) ΔV²`.

The approximation it makes is dropping the brine's shear stiffness, and here
that is `G = 0.44 MPa` against the ice's `3550 MPa` — 0.012 % of the matrix.

This is a real route, and the measurements below say it is now worth building.
Against Abaqus's own hybrid results the undrained layered cells lose 8–12 % of
their across-layer modulus and the loss grows as the mesh improves, while the
two cheap alternatives — element-level B-bar and order 2 — are respectively a
no-op on C3D4 and unconvergent. See *Refine, and the gap opens*.

Which leaves the two routes below, and the cheap alternative of solving those
cells at order 2.

**Measured, it does not.** One frozen packing (`SPAX_SAVE_PACKING`), the same
geometry throughout, inclusion fraction 0.295:

| Case | `E_eff` order 1 → 2 | `G_eff` order 1 → 2 |
|---|---|---|
| brine, `nu` = 0.4900 | +3.62 % | +4.59 % |
| brine, `nu` = 0.49993 | +3.78 % | +4.82 % |
| **compressible twin, `nu` = 0.300** | **+4.03 %** | **+4.51 %** |

The twin has the *same* inclusion Young's modulus (`E = 9KG/(3K+G) = 1.320e8`)
and differs only in compressibility. Its order-1 → order-2 change is the same
as the brine's, so that change is ordinary linear-tet over-stiffness, not
locking. Had the missing hybrid formulation been biting, the near-incompressible
rows would have moved much further than the twin.

Refining the order-2 brine mesh (`L_mesh` 0.050 → 0.035) moves `E_eff` by
0.22 % at `nu` = 0.490 and 0.23 % at `nu` = 0.49993 — converged, and equally so
at both.

The reason is physical, and worth stating because it means the result is
structural rather than luck. The brine's ν = 0.49 comes from a large K/G ratio,
not from being volumetrically stiff:

| | K | G | E | ν |
|---|---|---|---|---|
| ice | 9.25 GPa | 3.55 GPa | 9.43 GPa | 0.330 |
| brine | 2.20 GPa | 0.044 GPa | 0.132 GPa | 0.490 |

The brine is **four times more compressible in bulk than the ice around it**,
and 80× softer in shear. Volumetric locking bites when a phase has to represent
isochoric deformation *while carrying load* — an incompressible phase acting as
a constraint on its surroundings. This one does the opposite: it is a soft,
nearly-void inclusion that the matrix flows around. It neither carries the load
nor constrains the matrix volumetrically, so the element's ability to represent
isochoric motion inside it barely enters the effective stiffness. The
incompressibility *is* being represented — `nu_eff` reads 0.336 with the brine
against 0.307 with the compressible twin — it simply is not what governs.

That also says exactly when the conclusion would stop holding: a soft phase
that **percolates** and carries load, or geometric **confinement**. Both are
present in the layered decks — see below.

### The layered cells: where it does show up

`layered_incompressible.sh` tests the case the spherical result does not cover.
The layered decks carry drained/undrained pairs:

| | K | G | ν | Abaqus element |
|---|---|---|---|---|
| drained | 2.2 MPa | 0.44 MPa | 0.406 | C3D4 (below the 0.45 threshold) |
| **undrained** | **2.2 GPa** | 0.44 MPa | **0.49993** | **C3D4H — hybrid** |

The undrained cell has everything the spherical one lacked: the brine is a
cell-spanning slab, confined between ice plates, with the ice bridges as
constrictions — and K = 2.2 GPa is only 4× below the ice, so it genuinely
carries load. Pressure transmission through that confined brine is the
mechanism the layered closure rests on. The campaigns also run it on linear
tets (`5196ff1`, "Run the layered sweeps with linear elements to fit the
face-constraint limit").

One geometry per element order; the drained twin made by rewriting the single
inclusion elastic card on the **identical mesh**, so drainage is the only
difference:

| order 1 vs order 2 | E_x (across layers) | E_z (in-plane) |
|---|---|---|
| drained, ν = 0.406 | **+9.08 %** | +0.75 % |
| undrained, ν = 0.49993 | **+11.19 %** | +1.06 % |
| **attributable to incompressibility** | **+2.11 points** | +0.32 points |
| *spherical brine vs its twin, for contrast* | *−0.41 points (none)* | — |

**Read the control before the headline.** The raw +11.19 % looks like locking
and mostly is not: the drained cell, which Abaqus also meshes without hybrid
elements, loses +9.08 % on the same geometry. That is ordinary linear-tet
stiffness amplified by thin slabs and narrow bridges, and **Abaqus carries it
too** — it is not a CalculiX-versus-Abaqus difference at all.

What *is* CalculiX-specific is the excess: **+2.1 points** in E_x, the part
Abaqus recovers with C3D4H and CalculiX cannot. Small next to the geometry
term, but real, and qualitatively unlike the spherical case where the
incompressibility penalty was zero (−0.41 points, i.e. the brine tracked its
compressible twin exactly). It is also directional — 2.1 points across the
layers against 0.3 in-plane — which is the signature of the confined phase
having to deform near-isochorically under pressure.

So the expected CalculiX-vs-Abaqus discrepancy on an undrained layered cell at
order 1 is of order **2 % in the across-layer modulus**, not 11 %.

**Exposure: 204 undrained cells across 14 decks** (`rve_bracket_*`,
`rve_layerb`, `rve_layercol*`, `rve_layermesh`, `rve_layerskel`,
`rve_weibull_layer`). The drained-only decks (`rve_eringen_layer`,
`rve_nlgeom_layer`, `rve_torsion_layer`) sit at ν = 0.406, below the hybrid
threshold, so Abaqus used no hybrid there either and CalculiX matches it.

### Do not rely on the numbers above: that cell was under-resolved

The table was measured on a cell whose brine slab is barely one element thick.
The generator refines to `lc_fine = 0.4 x L_mesh` near an inclusion, so:

| | slab thickness | `lc_fine` | elements through the slab |
|---|---|---|---|
| the test above, L_mesh=0.020 | 0.0075 | 0.0080 | **0.9** |
| the test above, L_mesh=0.011 | 0.0075 | 0.0044 | **1.7** |
| campaign `rve_layerb` | 0.0100 | 0.0020 | **5.0** |
| campaign `rve_layercol_p060` | 0.0075 | 0.0022 | **3.4** |

With one element across it, the confined brine layer — the whole feature under
test — is not represented at all. Both the +9.08 % geometric term and the
+2.11-point incompressibility excess were measured on a cell that cannot carry
the mechanism they are about, and the +7.40 % order-2 drift under refinement is
at least as likely to be the unresolved slab as locking. **These numbers do not
transfer to the campaign's cells.**

`layered_incompressible.sh` now takes `NSLABS` and `SLABVOF` so the local
resolution can be held at the campaign's while the cell is shrunk to something
this machine can solve — `L=0.15 NSLABS=2 LMESH=0.0045` gives 4.2 elements
through the same 0.0075 slab. That is the measurement to trust.

### Re-resolved: the incompressibility excess is not there

`L=0.15 NSLABS=2 SLABVOF=0.1000 LMESH=0.0045`, 4.2 elements through the slab,
inside the campaign's own 3.4–5.0. One tetrahedralisation per order — 360 435
elements at order 1 (178 217 equations) and the same tets carrying mid-side
nodes at order 2 (1 449 110 equations) — with the drained twin the identical
mesh and one elastic card rewritten. Every solve converged (equilibrium gap
8.4e-8 to 2.7e-6). `out_layerres/`:

| | E_x (across layers) | E_z (in-plane) |
|---|---|---|
| undrained, ν = 0.49993, order 1 → 2 | +2.073 % | +0.678 % |
| drained, ν = 0.406, order 1 → 2 | +3.916 % | +0.431 % |
| **excess attributable to incompressibility** | **−1.84 points** | +0.25 points |

Against the under-resolved cell, which gave +11.19 / +9.08 / **+2.11 points**.
Two things changed, and both matter:

* **The raw order-1 error shrank 5×**, from +11.19 % to +2.07 %. Most of what
  looked like a locking signal was the unresolved slab, exactly as suspected.
* **The excess changed sign.** At campaign resolution the undrained cell loses
  *less* going from order 1 to order 2 than its drained twin does. There is no
  order-1 stiffness penalty left for a hybrid element to recover.

**Do not read the −1.84 as incompressibility helping.** It says the paired
design's control is imperfect, and the reason is visible in the moduli: E_x is
3.36 GPa drained against 6.35 GPa undrained. Filling the slab with a phase four
times more compressible than the ice does not just remove a locking mechanism,
it moves the load path off the bridges — and the bridges are what the linear
tet resolves badly. So the drained twin does not carry the *same* geometric
discretisation error, and the subtraction that isolates locking rests on
assuming it does.

What survives that objection is the direction. Locking makes order 1 too stiff;
if it were biting, the undrained cell's order-1 error would exceed the drained
one's. It is smaller, by more than the in-plane term is large.

**And that inference is wrong, for a reason this comparison structurally cannot
see.** An order-1-vs-order-2 test can only detect an error the two orders do
*not* share. C3D10 locks on this cell too — the order-2 convergence sweep
already refused to certify it — so both orders are stiff, the difference
between them says little, and the whole method is blind here. The direct
measurement against Abaqus's hybrid element, two sections down, puts the cost
at 8–12 %. Read *Refine, and the gap opens* before using anything in this
subsection.

One caveat stands either way: this is one morphology, `n_bridges=2`,
`bridge_fraction=0.29`. Narrower bridges mean tighter constrictions, so the
penalty has to be re-checked at the low-`bridge_fraction` end of
`rve_bracket_bridge.csv` before any number is relied on.

The other caveat — that this infers the hybrid benefit from an order-1-vs-order-2
comparison rather than measuring C3D4 against C3D4H directly — no longer stands.
See below.

### Measuring C3D4 against C3D4H, with no Abaqus licence

Everything above compares element *orders* inside CalculiX. That can only say
the two differ. The question is what Abaqus gets from C3D4H that CalculiX
cannot, and the tree already holds Abaqus's answer: the layered campaigns
stored drained *and* undrained `E_x` for hundreds of cells.

The obstacle is the one `validate_gas.sh` hit — the campaign's packing seed is
not recorded, so a fresh cell is not the same cell. The way through is that the
decks come in **drained/undrained pairs**:

```
R = E_x(undrained) / E_x(drained)
```

Both codes mesh the drained cell (ν = 0.406) with the plain C3D4. Only the
undrained cell (ν = 0.49993) differs — C3D4H in Abaqus, C3D4 here. Geometry
enters `R` through a ratio, where it largely cancels, and the stored tables
carry two or three **seeds** per condition, so **Abaqus's own seed spread in
`R`** is the noise floor the difference has to clear. Locking makes a
displacement element too stiff and only the undrained cell can lock, so it
inflates `R`.

`layered_abaqus_ratio.sh` runs it; `report_abaqus_ratio.py` reports it without
re-solving.

**`rve_layermesh` at the campaign's coarse level** (`L_mesh` = 0.0240, L = 0.50,
345 029 elements, 181 850 equations, both solves converged at ~1.3e-6):

| | `E_x` undrained | `E_x` drained | porosity | `phi_soft` |
|---|---|---|---|---|
| **CalculiX C3D4 / C3D4** | 6.3161e9 | 2.5170e9 | 0.00991 | 0.13021 |
| Abaqus C3D4H / C3D4, s1 | 6.2534e9 | 2.5394e9 | 0.01034 / 0.01014 | 0.12965 / 0.12893 |
| Abaqus C3D4H / C3D4, s2 | 6.1903e9 | 2.4985e9 | 0.01013 / 0.00985 | 0.12944 / 0.12852 |

Read the geometry columns first: the achieved phase fractions agree to
0.5–1.3 %, so unlike the gas comparison these really are equivalent packings
and even the **absolute** moduli can be read. They agree to +1.0…+2.0 % on the
undrained cell and −0.9…+0.7 % on the drained one, against an Abaqus seed
spread of 1.0 % and 1.6 % respectively. **A CalculiX C3D4 undrained layered
cell reproduces an Abaqus C3D4H one to within its own packing noise.**

The ratio sharpens that:

| | `R` | Abaqus seed spread | excess |
|---|---|---|---|
| `LMESH_m0p0240` | 2.5094 vs 2.4701 | 0.61 % | **+1.59 %** |

So there *is* a signal, and it has the sign locking predicts — CalculiX reads
1.6 % stiffer in the undrained-to-drained ratio than Abaqus's hybrid element
does, against a 0.61 % noise floor. **1.6 % is the number to carry**, not the
+2.1 points the order comparison suggested and not the 11 % the under-resolved
cell suggested.

### The campaign already answered the resolution question — in Abaqus, with C3D4H

Before reading anything into 1.6 %, note what `9104261` (2026-08-13, *"Mesh
gate result: size every layered cell per its own layer thickness"*) measured.
It is the same convergence study, run **in Abaqus with the hybrid element**, at
φ = 0.10, b = 0.293:

| elements across the layer | drained error | undrained error |
|---|---|---|
| 3.0 | 0.0 % | 0.0 % |
| 2.2 | 0.3 % | 0.3 % |
| 1.5 | 1.0 % | **20.7 %** |
| 0.7 | 8.7 % | **35.0 %** |

**The undrained response needs about twice the resolution the drained one
does, and C3D4H does not rescue it.** A mesh that resolves the drained modulus
to 1 % is 21 % out on the undrained one — in Abaqus, with the mixed
formulation, on the very cells this section is about. Whatever makes a confined
near-incompressible layer hard to discretise, the hybrid element is not the
cure for it.

That gate is also the cross-check that `R` is the right metric, because the
campaign reached for the same quantity without calling it that: *"the
undrained/drained ratio there is 2.47 against 1.99 converged."* Those are `R`
at `L_mesh` 0.0240 and 0.0060 in `results_layermesh.csv` — 2.4701 and 1.9897.
The ratio this comparison uses to cancel packing is the campaign's own
convergence measure.

So the coarse level has to be read in that light. `L_mesh` = 0.0240 is 0.7–1.3
elements across the layer, where the gate says **Abaqus is 35 % out**. The
CalculiX-vs-Abaqus excess measured there is 1.6 %. At that mesh the
discretisation error swamps everything and the two codes agree because they are
both wrong in the same way — which is exactly why the comparison has to be
repeated where Abaqus starts getting it right.

### Refine, and the gap opens: it is locking after all

Repeating at `L_mesh` = 0.0120 (1 211 410 elements) does not shrink the excess.
It multiplies it by five.

| | `R_ccx` (C3D4/C3D4) | `R_abq` (C3D4H/C3D4) | Abaqus seed spread | excess |
|---|---|---|---|---|
| `LMESH_m0p0240` | 2.5094 | 2.4701 | 0.61 % | +1.59 % |
| `LMESH_m0p0120` | 2.5852 | 2.3786 | 0.75 % | **+8.69 %** |

The mechanism is visible in how each code responds to the refinement. Mean over
the Abaqus seeds, against the single CalculiX cell:

| refining 0.0240 → 0.0120 | Abaqus | CalculiX |
|---|---|---|
| drained (both C3D4) | −7.11 % | −3.97 % |
| **undrained** (C3D4H vs C3D4) | **−10.55 %** | **−1.06 %** |

On the drained cell both codes refine downward at comparable rates — ordinary
convergence, and the control that says the meshes are comparable. On the
undrained cell Abaqus moves 10.6 % toward its converged answer and **CalculiX
moves 1.1 %.** The displacement element is stuck. Refining the mesh does not
relieve the volumetric constraint, which is what volumetric locking *is*.

Read as an absolute disagreement on `E_x`:

| | undrained | drained |
|---|---|---|
| `L_mesh` = 0.0240 | +1.52 % | −0.08 % |
| `L_mesh` = 0.0120 | **+12.28 %** | +3.30 % |

The achieved phase fractions match throughout (`phi_soft` 0.12952 against
0.12834–0.12962), and every solve converged (equilibrium gap 6.6e-7 to 6.3e-6).

**This corrects what the sections above concluded.** "The missing hybrid
element costs nothing measurable" was measured two ways that both hide it: an
order-1-vs-order-2 comparison inside CalculiX, which cannot see an error both
orders share, and a single coarse mesh, where discretisation error dominates.
Against Abaqus's own hybrid results, on the campaign's own decks, at the mesh
density the campaign now specifies, **it costs 8–12 % of the across-layer
undrained modulus, and the cost grows as the mesh improves.**

Two more morphologies, both at `L_mesh` = 0.0240 (`rve_bracket_bridge`, three
seeds):

| | `R_ccx` | `R_abq` | spread | excess |
|---|---|---|---|---|
| `BRKB_b020` (`bridge_fraction` = 0.02) | 19.5023 | 19.3558 | 7.99 % | +0.76 %, inside |
| `BRKB_b280` (`bridge_fraction` = 0.28) | 3.0425 | 2.9150 | 2.07 % | +4.37 % |

`b020` is the tightest constriction in the tree and the README's flagged worst
case; it comes back inside the noise, but its noise floor is 8 % because the
drained cell there is nearly a stack of disconnected plates (`E_x` = 0.28 GPa)
and its packing scatter is large. It is a weak "inside", not a clean one. Both
sit at the retired coarse setting, so by the trend above both understate.

**Consequence for the 204 undrained cells.** Production sizing is now
`L_mesh = t/2.5` clipped to [0.005, 0.012] — at or finer than the level where
the gap is 8.7 %. Those cells cannot be solved in CalculiX at production
resolution without a mixed formulation, and the drained cells are unaffected.
This is the one place in the repository where the missing hybrid element is
disqualifying rather than merely measurable.

**And it kills the cheap fixes.** B-bar is a no-op on C3D4 — one integration
point, so the element mean of the divergence is the pointwise value
(`../patches/0002-bbar-mean-dilatation.patch`, verified bit-identical). C3D10 has four
integration points, so B-bar is *not* trivial there — but see below, it does
not help either. What is left is nodal-averaged B-bar (needs `mastruct.c` and a
new assembly path), a true mixed element (38 files), or the hydrostatic fluid
cavity sketched under *Fluid elements* above.

### What `*USER ELEMENT` actually is

CalculiX has **no Abaqus-style UEL**: there is no separately compiled
subroutine you supply and link, no `AMATRX`, no `RHS` callback
(`grep -rn amatrx src/*.f src/*.c` returns nothing). Any new element means
editing ccx's Fortran and rebuilding.

But this file previously called `*USER ELEMENT` "the substructure/superelement
interface", and that is wrong — substructures are only its fall-through branch.
It is a genuine user-element slot, and **the plumbing a mixed element needs is
already in place**:

| what a mixed u/p element needs | where ccx already has it |
|---|---|
| declare a new element type | `*USER ELEMENT, TYPE=Uxxxx, NODES=n, INTEGRATIONPOINTS=p, MAXDOF=d` — same syntax as Abaqus, parsed in `userelements.f` |
| raise the per-node DOF ceiling | `allocation.f:1109-1110`, `mi(1)=max(mi(1),iuel(2,id))` and **`mi(2)=max(mi(2),iuel(3,id))`** — `MAXDOF` raises the global limit, so `mt=mi(2)+1` grows |
| room in the matrix structure | `mastruct.c:137` reads `ndof` per node from the element label, `:186` sets `nactdof[mt*node+k]=1` for `k=1..ndof` — **extra nodal DOFs beyond 3 are already supported** |
| get called during assembly | `mafillsm.f:167` recognises `lakon(1:1)=='U'` and reads `ndof`/`nope` from the label; `:302` calls `e_c3d_u` |
| get called for stress recovery | `resultsmech.f:115` calls `resultsmech_u` |

`e_c3d_u.f` is a hard-coded dispatcher on the type character — `U1` is a
Timoshenko beam, `US45`/`US3` are flat shells, anything else falls through to
the substructure path. So **adding a mixed element means adding a branch to
`e_c3d_u.f` and `resultsmech_u.f`**, not touching `lakon` handling across 38
files. The global-DOF work the earlier estimate feared is done.

**The one real awkwardness** is that `ndof` is uniform across every node of the
element type, while Taylor-Hood P2/P1 — which is what `C3D10H` is — wants
pressure on the four corner nodes only. The practical route is to declare
`MAXDOF=4` on all ten nodes and constrain the six midside pressure DOFs to
zero, which gives genuine P2/P1 at the cost of six dead DOFs per element.
P2/P0 is the alternative that needs no constraining, and the section below
shows it is too weak to be worth it.

### B-bar on C3D10: it only acts where the mesh is already wrong

`0002-bbar-mean-dilatation.patch` advertises that it *"acts only at order 2
(C3D10, 4 integration points), where it moves E_x on the confined undrained
layered cell by −11.5 %"*. That is true and it is misleading, because of which
cell it was measured on.

| C3D10 undrained, `CCX_BBAR_NU` = 0.45 | elements through the slab | `E_x` change |
|---|---|---|
| `out_layerinc` (L = 0.30, `L_mesh` = 0.015) | 0.9–1.7 | **−11.52 %** |
| `out_layerres` (L = 0.15, `L_mesh` = 0.0045) | **4.2**, the campaign's own | **−0.003 %** |

Same patch, same threshold, same morphology, both confirmed active (ccx echoes
`SPAX B-bar active for nu >= 0.45` and the equilibrium gap flips from 3.96e-7
to exactly 1.00, which is the documented signature of the switch engaging).
`E_z` moves −0.069 % on the resolved cell.

**Element-level B-bar does something on C3D10 only where the slab is
unresolved — which is precisely where the answer is wrong for an unrelated
reason.** Resolve the slab and C3D10 is not locking, so there is no volumetric
penalty for the mean-dilatation term to remove. The −11.5 % is not a recovery;
it is one under-resolved wrong answer moving toward a different wrong answer.

The reason is structural, not a defect in the patch. Element-mean B-bar on a
quadratic tet is a P2/P0 pressure approximation: it imposes one dilatation
constraint per element against 30 displacement DOFs. On a linear tet that one
constraint is the whole element (and hence trivially satisfied); on a quadratic
one it is too weak to change a field that was not badly constrained to begin
with. Making it bite needs the pressure space refined with the displacement
space — P2/P1 Taylor-Hood — which is the true mixed element, not a B-bar.

**Caution when using the patch at all.** With B-bar active only the stiffness is
patched; `resultsmech.f` still recovers internal forces with the standard B, so
`equilibrium_gap` reads ~1.0 and the repository's only reference-free
convergence check is unavailable on exactly the runs that most need it. The
volume-averaged stress stays consistent — the element mean of λ·div(u) *is* the
B-bar pressure — so `E_eff` is readable, but always solve an unpatched twin
beside it as the control.

**One limitation to keep in view.** This validated `E_eff`, `G_eff` and
`nu_eff` — homogenised quantities, which average over the inclusion interiors.
Volumetric locking distorts the *local* stress field well before it moves a
homogenised modulus. The per-slice principal-stress and SCF extraction is not
ported to CalculiX (see *Not ported*), and if it ever is, this measurement does
not carry over to it: the fields inside a near-incompressible phase would need
checking on their own terms.

**For the spherical decks: use `SPAX_MESH_ORDER=2` and the missing hybrid
element costs nothing measurable.** Order 1 is a separate 4 % error that Abaqus
has too. This was measured at one inclusion fraction on one packing.

**For the undrained layered decks that conclusion does not hold** — see
*Refine, and the gap opens* above, where the direct comparison against Abaqus's
hybrid results puts the cost at 8–12 % and rising with mesh quality. The
condition named right here as the one that would break the spherical result —
a soft phase that percolates and carries load — is exactly the condition those
cells satisfy.

### If it ever does need building

Two routes, and they are not the same size of job.

**A true mixed displacement/pressure element** — an extra pressure unknown
entering the global system — was estimated here as touching 38 source files,
because `lakon` is consumed in that many. **That estimate was wrong**, and the
reason is the user-element slot: see below.

**B-bar (mean dilatation) is the tractable one, and it is the standard fix.**
`e_c3d.f` assembles the isotropic stiffness directly from the Lamé constants,
with the volumetric and deviatoric parts already written as separate terms:

```fortran
s(ii1,jj1) = s(ii1,jj1) + (al*w(1,1) +
&              um*(2.d0*w(1,1)+w(2,2)+w(3,3)))*weight
```

Locking lives entirely in the `al` (λ) term, which blows up as ν → 0.5. B-bar
replaces the discrete divergence operator in that term by its element average —
accumulate it through the integration loop, then add the volumetric
contribution once afterwards instead of at every point. No extra DOF, no change
to the matrix structure, no solver changes; roughly 50–100 lines in `e_c3d.f`
plus registering a type name, and a matching change in `resultsmech.f` so the
recovered volumetric stress uses the same operator.

**What B-bar would and would not buy here.** For the spherical decks, nothing:
the 4 % order-1 error is not volumetric — the compressible twin shows the same
4 % — so B-bar would not let those cells drop from quadratic to linear
elements, and there is no incompressibility penalty to recover.

For the **undrained layered** decks it would buy the +2.1 points measured
above, which is the part Abaqus gets from C3D4H. That is the one place in this
repository where the missing hybrid element costs something measurable. It is
worth weighing against the alternative, which costs nothing to try: solve those
cells at `SPAX_MESH_ORDER=2`, where the penalty is by construction absent.
The campaigns moved to linear elements for a reason (`5196ff1`, the
face-constraint limit), so that alternative may not be free in practice — but
it should be priced before writing Fortran.

## Two ccx incompatibilities worth knowing

**Node sets in `*EQUATION`.** CalculiX resolves set names everywhere else but
wants a node *number* inside `*EQUATION`. Untranslated it stops with one
`*ERROR reading *EQUATION` per periodic-BC equation — loud, at least.

**The 20-character field limit.** CalculiX copies each comma-separated field
into a fixed 20-character buffer and **truncates** anything longer. Python's
`repr` overruns it routinely: a coordinate of `3.8163916471489756E-17` is 22
characters and truncates to `3.8163916471489756E-`, which fails to parse.

That failure is the lucky case. Truncation only errors when it lands
mid-exponent — a field that truncates to a still-valid number is accepted
silently at the wrong value. An `*Equation` coefficient of
`3.469446951953614e-18` (21 characters) truncates to `3.469446951953614e-1`, a
periodicity constraint off by seventeen orders of magnitude in a deck that
solves without complaint. The converter therefore re-renders **every** numeric
field at 12 significant digits, not only the ones observed to break.

## Scale, and why the default is the iterative solver

A stock `ccx` is usually linked against SPOOLES only — check with
`ldd $(which ccx)`. SPOOLES is a *direct* solver, and at production size that
is not a trade-off but a wall:

| equations | SPOOLES | ITERATIVE CHOLESKY |
|---|---|---|
| 167 306 | 23.8 s, 1276 MB | 20.6 s, **370 MB**, 338 iters |
| 488 390 | 149.8 s, 4402 MB | 64.8 s, **1072 MB**, 306 iters |
| ~1.4 M | still factorising at 19 GB after ~25 min, stopped | — |

Two things matter in that table. The iterative memory is **linear** in the
model (2.90× for 2.92× the equations) while SPOOLES' is super-linear and its
time worse still (6.3× for 2.92×). And the **iteration count did not grow** —
338 → 306 — so the incomplete-Cholesky preconditioner is holding up as the cell
grows, which is what makes the method usable at all here.

So `SPAX_CCX_SOLVER` now defaults to `ITERATIVE CHOLESKY`. Set it to `SPOOLES`
for a small-cell reference, or `DEFAULT` to write no `SOLVER` parameter and let
ccx choose.

### The tolerance, which is not optional

Out of the box the iterative solvers return the **wrong answer, quietly**. The
convergence test in `pcgsolver.c` stops once the residual max-norm falls to
0.5 % of the mean load (`c1 = 0.005`, a local constant no input deck reaches),
which leaves `E_eff` 0.15 % below the direct answer on the same deck.

`../patches/0001-iterative-tolerance.patch` exposes it as `CCX_ITER_TOL`:

| `CCX_ITER_TOL` | iters | wall | peak | `E_eff` error | equilibrium gap |
|---|---|---|---|---|---|
| 0.005 (stock) | 70 | 14.0 s | 370 MB | **−0.1477 %** | 1.8e-3 |
| 1e-3 | 94 | 14.6 s | 370 MB | +0.0133 % | 3.3e-6 |
| 1e-5 | 338 | 20.9 s | 370 MB | −0.0000 % | 1.8e-7 |
| SPOOLES | — | 23.2 s | 1275 MB | reference | 7.4e-8 |

`1e-5` is the default `SpaX_CalculiX.solve()` exports. Memory is flat in the
tolerance — tightening costs iterations, not storage.

Beware the decoy: the `eps` argument threaded down from `preiter.c` looks like
the tolerance and is documented as "required accuracy", but is never read — PCG
and CG overwrite it with the final residual and return it. Patching `eps`
changes nothing.

An unpatched ccx silently ignores `CCX_ITER_TOL`, so `solve()` checks the log
for the tolerance the patched binary echoes and warns plainly when it is
absent.

### The check that survives into production

`equilibrium_gap`, recorded per RVE in `results.csv`, is the relative
disagreement between the volume-averaged stress and the reference-point
reaction — two independent measurements of the same macroscopic stress, linked
by nothing but the model being right and the system being solved. It tracked
the true error across the whole sweep above (1.8e-3 when the solve was bad,
1.8e-7 when good) and it needs **no reference solve**, which is the point: on a
cell too large to solve directly even once, it is the only convergence evidence
available. Treat anything above ~1e-4 as a solve not to be believed.

## Not ported

`extract_principals` — the per-slice principal-stress and SCF field CSVs behind
`results_scf.csv` and the stress-field figures. `SpaX_CalculiX.extract_averages`
supplies the volume-averaged stress and strain tensors that the 6x6 elasticity
route needs, but not the slice fields.
