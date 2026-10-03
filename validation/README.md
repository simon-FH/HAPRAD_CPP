# Validation harness

Compares this C++ implementation against the **original HAPRAD 2.0 Fortran** in
[`../haprad2`](../haprad2) -- the code it was ported from. Same formalism, same
models, so any disagreement is a defect in the translation rather than a
modelling choice.

## Build and run

```bash
export CERN_LIB=$HOME/cernlib/lib     # needed for the Fortran reference only
make
./compare.py --grid rge               # or --grid eg2
```

`rc_point` (C++) is always built. `haprad2_point` (Fortran) needs PDFLIB for its
PDF x FF semi-inclusive model, so it is built only when `CERN_LIB` is set; see
the CERNLIB section of the top-level README. Without it you still get Tier 0.

Both drivers take `E x Q2 z pt phi_deg` and print `@name value` lines, so they
can be diffed directly:

```bash
./bin/rc_point      10.5473 0.234 2.5 0.34 0.424 60
./bin/haprad2_point 10.5473 0.234 2.5 0.34 0.424 60
```

Run them from `data/`, which holds symlinks to the inputs both codes open by
relative path (`pi_n_maid.dat`, the `*.grid` fragmentation functions).

`bin/phi_convention "<file glob>"` checks that RG-E's `phi_PQ` is the `phi_h`
HAPRAD expects, by comparing `V1 = 2 k1.p_h` from the four-vectors against
`THadronKinematics::V1()`. It also classifies every event HAPRAD rejects.

### Structure-function tables

The C++ needs a table of azimuthal amplitudes (format in
`../TSemiInclusiveModel.h`); without one sigma_Born is exactly zero.

```bash
./bin/make_toy_table /tmp/toy.root reduced          # known analytic amplitudes
./compare.py --grid rge --table /tmp/toy.root       # also prints table coverage
HAPRAD_SI_TABLE=/tmp/toy.root ./bin/rc_point 10.5473 0.234 2.5 0.34 0.424 60
```

With a table, `rc_point` also reports `sf_lookups`, `sf_oor`, `sf_empty` and
`sf_unphys`: how many structure-function lookups the integrals made in that call,
and how many fell outside the table, in an unfitted cell, or below threshold.
`HAPRAD_SI_NEAREST=1` switches the reader from interpolation to nearest-cell
lookup, for comparison only -- see PLAN.md 0.2 for why that mode breaks the
inner R integral.

### Tables from HAPRAD 2.0's own model (PLAN.md 2.3)

Tier 3 can only agree if both codes see the same physics. These tools build a
table *from* the FORTRAN's PDF x FF model, so it can be fed to the C++:

```bash
B="18:0.5:10,16:1:10.5,20:0:1,24:0:1.5"            # nQ2:lo:hi,nnu:..,nz:..,npt2:..
./bin/make_model_table points $B 10.5473            > cells.txt
(cd data && ../bin/born_harmonics < ../cells.txt)  > harm.txt
./bin/make_model_table build  $B 10.5473 harm.txt table.root
MODEL_TABLE_INTERP=log ./table_convergence.py       # Tier 3 vs grid resolution
```

* `born_harmonics` is HAPRAD 2.0's Born cross section alone -- the setup of
  `ihaprad()` and the vv10/vv20 part of `sphih()`, copied without the tail
  integrals. It matches the full `haprad2_point` to 6e-13. It links a copy of
  `semi_inclusive_model.f` with the unconditional `stop` at line 178 removed
  (the build fails if the substitution does not happen); the reference used by
  `haprad2_point` is untouched.
* `make_model_table build` divides the FORTRAN harmonics by the C++'s own
  kinematic factor K (sigma_Born with a unit table, via
  `TRadCor::CalculateBorn`), then checks the result through the real reader at
  every cell centre. That check is not exact: the C++ inversion drops the
  `-m_e^2 lambda_q` term of lambda, which reaches 1e-5 at Q2 = 0.4, y = 0.98.
* `MODEL_TABLE_INTERP=log` writes `interpolation = "log"` into the table, so the
  reader interpolates log(A) and the ratios Ac/A, Acc/A instead of the raw
  values.

### Producer closure (PLAN.md 2.2)

```bash
./closure_producer.py [n_events] [seed]
```

Toy pions in the RG-E ntuple format with a known azimuthal modulation
(`bin/make_toy_events`) go through `MakePhiTable` fill and fit; the fitted
modulation in each cell is compared with the truth as pulls. Expected: mean 0,
width 1. Average a few seeds -- one seed's 187 cells only pin the mean to
+-0.07.

### Normalisation closure (PLAN.md Phase 1)

```bash
./closure_normalisation.py [n_events]
```

Toy pions uniform in (Q2, nu, z, pt2, phi), weighted by HAPRAD 2.0's Born
cross section (`bin/make_weighted_events` + `bin/born_harmonics`), go through
`MakePhiTable` with the generator's luminosity, on 4^4, 6^4 and 8^4 grids.
Per cell it reports producer / model (`make_model_table`, end to end), true
cell average / model (bin-centring), and producer / true cell average
(`MakePhiTable` alone). Expected: the last is 1 for A, and 1 within errors for
Ac/A and Acc/A, flat in every variable. 3M events take about 4 minutes.

### Iteration closure (PLAN.md Phase 4)

```bash
./closure_iteration.py [n_events] [--iterations 4] [--jobs N] [--keep DIR]
```

Toy pions weighted by HAPRAD 2.0's Born cross section give a true table T_B
(`MakePhiTable`) and its RC grid; the same events weighted by Born x RC are the
"observed" data, and `Utilities/rc_iterate.py` runs on them. T_B is an exact
fixed point of the loop and the events are the same throughout, so the
distance of table_k from T_B and of grid_k from the true grid measures the
convergence alone. About 12 minutes for 2M events and 4 iterations on 14
cores (six RC grids of 5,250 nodes).

## The tiers

They are ordered so that a failure in one invalidates the ones below it.

| Tier | What | Depends on the SF model? |
|---|---|---|
| **0** | C++ self-consistency: phi symmetry, `sigma_Born` varies with z and pt, everything finite | no |
| **1** | Kinematics vs Fortran: invariants, hadron kinematics, tau limits | no |
| **2** | Exclusive radiative tail vs Fortran | **no** |
| **3** | Cross sections and RC factors vs Fortran | yes |

Tier 2 deserves the emphasis: `tai_ex` is built from the MAID exclusive model
(`TSffun`) and the theta matrix alone and never touches
`TStructFunctionArray`, so it is fully comparable **today**, while the
semi-inclusive model is still missing.

## Current status

```
                no table                      --table toy_reduced.root
Tier 0 : pass   (4 of 5 checks SKIP)          pass, all 5 checks run
Tier 1 : pass   ~4e-14 at 6 and 10.5 GeV      unchanged
Tier 2 : pass   3e-5 .. 1.5e-4                unchanged
Tier 3 : FAIL   sigma_Born = 0                FAIL by construction: a toy table and
                                              HAPRAD 2.0's PDF x FF model are
                                              different physics (PLAN.md 2.3)
```

Tier 1 passing at RG-E energies is the substantive good news: the whole
kinematics layer is a faithful port, so the remaining defects are localised to
the models and the integration.

### What Tier 2 found

Tier 2 failed on first run, with the C++ exclusive tail a factor 2-3.5 below the
Fortran. Dumping the integrand `rv2tr(tau, phi_k)` pointwise from both codes
(`HAPRAD_DUMP_EXC=1`, the `EXC` block) showed it agreeing to **2e-7** over 45
points, and swapping the integrator changed nothing -- so neither the integrand
nor the quadrature was at fault.

The cause was in `TRadCor::Initialization()`. It computed the normalisation as

```cpp
fHadKin->Evaluate();
...
if (fKin->T() >= 0) { N = N * sqrt(lambda_q) / (2 M p_l); }
```

but `THadronKinematics::SetMomentum()` overwrites `fKin`'s `T` with the computed
invariant `t`, which is **negative**, whenever `p_t` was the quantity supplied.
The test on the next line was therefore always false and the p_t-differential
Jacobian was never applied. The FORTRAN applies it inside its
`IF (tdif .GE. 0)` block, before `tdif` is overwritten (`ihaprad.f`).

The predicted deficit `sqrt(lambda_q) / (2 M p_l)` reproduced the observed ratios
to four decimals (2.1753 vs 2.1749; 3.1370 vs 3.1368).

**Scope of the bug.** It scaled every absolute cross section -- `GetSigBorn`,
`GetSigObs`, `GetTail` -- by a factor of 2 to 3.5 over the RG-E range. It did
**not** affect the radiative correction factors, because `N` multiplies
sigma_Born and both tails alike and therefore cancels in `GetFactor1/2/3`. This
was checked directly, not just argued: with the structure functions stubbed to
constants so that sigma_Born is real, the fix moves `sib` from 8.5670e4 to
2.6875e5 (ratio 3.1370, the Jacobian) while `f1` stays at 1.048191. So the
numbers `GetRC` writes out were never wrong on this account.

A second, smaller change came out of the same investigation: the exclusive tail
now uses deterministic adaptive cubature (Genz-Malik, the same family as the
FORTRAN's NAG D01FCE) rather than MISER Monte Carlo capped at 20000 calls. That
was *not* the bug -- both integrators gave the same answer -- but on this
sharply peaked integrand MISER agreed with the reference only to ~1.4e-3 where
adaptive reaches ~5e-5.

### What Tier 0 found (hurdle H4)

`tai_in` came back NaN at exactly phi = 180 deg, while 179 and 181 agreed to ten
digits. Instrumenting the integrand showed the theta matrix finite but the
structure functions returning `H(2) = inf` with the rest NaN -- the signature of
`0 * inf`.

The cause is in `TStructFunctionArray`. Because this model inverts the *measured*
azimuthal amplitudes, and those enter the Born cross section multiplied by
`b ~ p_t` and `b^2 ~ p_t^2`, the inversion carries the reciprocals:

```cpp
H3z = ... / Sqrt(tldPt2);
H4z = ... / tldPt2;
```

The guard above admits `tldPt2 == 0` (it tests `>= 0`), and the h_i then combine
an infinite `H4z` as `tldPt2 * H4z`. At phi_h = 180 deg the shifted p_t passes
exactly through zero inside the integration domain. The FORTRAN has no such
inversion -- `semi_inclusive_model.f` builds H3/H4 directly from PDFs x FFs -- so
it stays regular.

`Evaluate()` now zeroes `fArray` on entry, as `strf()` does with `sfm(1..8)`, and
refuses to publish non-finite h_i. That removes the NaN and the uninitialised
read (`TPODINL` constructs a fresh `TStructFunctionArray` per evaluation, so
anything unassigned is stack garbage).

**But that is a robustness fix, not a cure.** With the NaN gone, phi = 180 reads
5.3e7 against 1.98e4 at 179 -- the NaN was the tip of a real pole, growing like
1/(180 - phi):

```
phi=179      1.978841e+04        phi=179.999  8.555769e+05
phi=179.9    2.679992e+04        phi=180      5.327075e+07
```

The pole is an artifact of the *flat* amplitudes this code is currently fed.
Giving them their physical dependence -- `A_c ~ p_t`, `A_cc ~ p_t^2`, which is
what the Born cross section requires -- removes it completely:

```
phi=179      1.654784e+04        phi=179.999  1.654919e+04
phi=179.9    1.654917e+04        phi=180      1.654919e+04
```

and takes the GSL integrator warnings from 43586 to zero. So the integrator was
never at fault, and no integrator change was made. **H4 is a symptom of H1**:
it closes for real only when the structure-function model supplies azimuthal
amplitudes with the right p_t behaviour.

The Tier 0 "no spikes in phi" check exists to catch exactly this, and reports
SKIP until a live model makes it meaningful.

## Notes on comparability

Deliberate differences that will keep Tier 3 from ever reaching machine
precision, even once H1 is fixed:

* **H11** -- the C++ forms `sigma_B * exp(delta_inf) * (1 + delta_VR + delta_vac)`
  while the Fortran forces `delta_inf = 0` and uses `sigma_B * (1 + alpha/pi * delta)`.
  Equivalent to O(alpha^2), but not digit-for-digit. Measured (PLAN.md 2.3):
  it moves sigma_obs / sigma_Born by 0.2-0.4% at the RG-E test points.
* The Fortran `sig` already includes both tails; `haprad2_point` subtracts them
  so that `sig_obs` means the same thing on both sides.
* `rc_point` defaults its missing-mass threshold to 0, disabling the gate in
  `TRadCor::CalculateRCFactor`, because the Fortran driver calls `ihaprad`
  directly. Kinematic rejection still happens inside both codes, identically.
  Pass `(m_n+m_pi)^2 = 1.16168` as a 7th argument to reproduce production `GetRC`.
