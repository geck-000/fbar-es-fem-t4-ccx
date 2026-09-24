# F-barES-FEM-T4 tests

Two commands cover the suite:

```bash
python3 run_tests.py                  # element checks, no ccx needed
python3 run_tests.py --quick          # operators, walk and deck only
python3 run_tests.py --ccx /path/to/ccx_fbar
python3 meshconv.py --ccx /path/to/ccx_fbar --ns 10 20 30 40 --kg 500
```

| File | Role |
|---|---|
| `run_tests.py` | suite entry point: element checks, the Fortran tangent, and with `--ccx` the CalculiX checks |
| `verify_element.py` | every check that needs no ccx: operators and patch test (V), the `u3vol.f` walk against `S = E A^c` (W), finite-strain force (N), consistent tangent (T), spurious-mode census (S), deck generator walk (D), assembled-structure audit (A) |
| `verify_nltan.py` | the Fortran `u3nltan` against the prototype, through a gfortran-built per-edge driver |
| `verify_nlccx.py` | the finite-strain element through CalculiX `NLGEOM`; `--uniaxial` for the 20% and 50% stretch paths |
| `verify_pairing.py` | coercivity of the two pairings; `--metrics --ccx` adds the cross-check |
| `verify_convergence.py` | manufactured fields, meshes and error norms shared by the CalculiX checks |
| `smoothing_proto.py` | Python prototype of the smoothing operators |
| `mms_ccx.py` | the manufactured problem through CalculiX; `--sweep` runs n = 8 to 24 over K/G = 10 to 5000 and prints the rates |
| `meshconv.py` | mesh convergence on the resolved slab cell, three arms, with report |
| `make_slabconv.py`, `report_slabconv.py` | cell generator and report used by `meshconv.py` |
| `nl_tan_driver.f` | per-edge driver compiled by `verify_nltan.py` |

`verify_element.py` takes `--quick` (skip N, T and S), `--skip-tangent`,
`--skip-stability`, `--n` for the census size, and `--deck ORIG GEN` to run the
deck checks on existing decks instead of generating a small cell.

The two-element patch tests establish *consistency* (a uniform strain state is
reproduced to roundoff). Locking behaviour is covered by the stability census
and the mesh-convergence cell, because an unstable element passes a patch test.
