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
Tier 0 : FAIL   f1/f2/f3 = inf at phi=180 (hurdle H4: GSL QNG fails to converge)
                the other three checks SKIP -- they need a working model
Tier 1 : pass   agreement to ~4e-14 at both 6 GeV and 10.5 GeV
Tier 2 : pass   agreement to 3e-5 .. 1.5e-4  (was failing by a factor of 2-3.5)
Tier 3 : FAIL   expected: sigma_Born is uninitialised memory (hurdle H1)
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

## Notes on comparability

Deliberate differences that will keep Tier 3 from ever reaching machine
precision, even once H1 is fixed:

* **H11** -- the C++ forms `sigma_B * exp(delta_inf) * (1 + delta_VR + delta_vac)`
  while the Fortran forces `delta_inf = 0` and uses `sigma_B * (1 + alpha/pi * delta)`.
  Equivalent to O(alpha^2), but not digit-for-digit.
* The Fortran `sig` already includes both tails; `haprad2_point` subtracts them
  so that `sig_obs` means the same thing on both sides.
* `rc_point` defaults its missing-mass threshold to 0, disabling the gate in
  `TRadCor::CalculateRCFactor`, because the Fortran driver calls `ihaprad`
  directly. Kinematic rejection still happens inside both codes, identically.
  Pass `(m_n+m_pi)^2 = 1.16168` as a 7th argument to reproduce production `GetRC`.
