# HAPRAD_CPP for CLAS12 / RG-E — game plan

Companion to [H1_NOTES.md](H1_NOTES.md) and [validation/README.md](validation/README.md).
Status as of commit `3ad1a7b`.

---

## Where we are

```
Tier 0 : pass     self-consistency; every model-independent quantity finite
Tier 1 : pass     kinematics vs HAPRAD 2.0 FORTRAN, ~4e-14, at 6 and 10.5 GeV
Tier 2 : pass     exclusive radiative tail vs FORTRAN, 3e-5 .. 1.5e-4
Tier 3 : FAIL     blocked on H1 -- sigma_Born reads uninitialised memory
```

Everything that can be validated without a structure-function model now agrees
with the reference implementation. The QED machinery is not in question; the
physics *input* is missing.

### Settled

| | |
|---|---|
| CERNLIB | no longer required; PDFLIB 8.04 built at `~/cernlib` as insurance, and it is what makes `haprad2` runnable as the reference |
| Validation | 4-tier harness against the original FORTRAN, `validation/compare.py` |
| p_t Jacobian | `Initialization()` tested `fKin->T() >= 0` *after* `Evaluate()` had overwritten T; every absolute cross section was low by `sqrt(lambda_q)/(2 M p_l)`, a factor 2-3.5. RC factors were unaffected (N cancels) |
| H4 | NaN at phi=180 contained, and root-caused to H1 rather than to the integrator |
| Latent bugs | `fLepton` uninitialised, swapped setters, stale results between calls, non-finite RC factors, dead members, unchecked file open, empty build dependencies -- all fixed, all verified numerically inert |
| H1 design | confirmed with the port's author: deliberately data-driven, no external model. The A/Ac/Acc -> H1..H4 inversion **round-trips exactly**; nothing there needs porting |

### The one open blocker

`TStructFunctionArray` needs `A`, `Ac`, `Acc` per `(Q2, x_B, z_h, p_t)` cell.
Nothing currently produces them. The chain instead fits the *integrated* sample
and seds three constants into the source, which makes `sigma_Born` independent
of z and p_t by construction -- and is the direct cause of H4.

---

## The sequence

### 1. Binning study — hours, no dependencies

Choose the grid the producer writes on. Preliminary numbers, run 020026, 275737
pi+ DIS events:

```
grid                  cells  filled   N>=1000  median N  sig(Ac/A)
EG2-like 6x5x10x5      1500     959        55        22      0.302
moderate 4x4x5x4        320     254        36        79      0.159
coarse   3x3x4x3        108      96        27       194      0.102
```

`sigma(Ac/A) ~ sqrt(2/N)`; physical modulations are 5-20%, so the bar is
`sigma <~ 0.03`, i.e. **N >~ 2200 per cell**.

**`runs_all.txt` lists 384 runs and exactly one is processed.** Scaled to the
full set, the EG2-like grid reaches `sigma ~ 0.015` and a finer 8x6x10x6 reaches
`~0.021`. So: **choose the grid for the final dataset, not for what exists
today**, and make the binning a parameter so it can be re-run. Treat
`sqrt(2/N)` as a floor -- acceptance correction inflates it.

### 2. The producer — data -> `newphihist.root`

Start from `git show b303487^:PhiHist/phihist.cpp`, which is the original and
was deleted in `b303487`. It reads a 5-D table `(Q2, Xb, Zh, Pt, Phi, Val, Err)`
and fits `A + Ac cos(phi) + Acc cos(2 phi)` per cell, writing NTuple
`AAcAcc_data` with `Q2:Xb:Zh:Pt:A:AErr:Ac:AcErr:Acc:AccErr:ChiSQ`.

`GetCentroids` already builds acceptance-corrected per-bin `phi_PQ` histograms
and keeps only their mean; fit them instead of averaging them.

Watch for:

* the model grid is in **x_B and p_t**, while the analysis binning is in
  **nu and p_T^2**. Different variables, not just different edges
  (`bin[3] = sqrt(pt2)` in `TSemiInclusiveModel.cxx`).
* `MIN_BIN = 4` in the original: cells with fewer than 4 filled phi bins get no
  fit and vanish from the table, which `THnD` then reads back as zero.
* `SetBinError(j, Err*1.04)` -- a 4% inflation worth understanding before copying.
* this step absorbs the rest of H3 (RG-E branch names, targets, vertex cuts).

### 3. Re-enable the consumer

* uncomment `SemiInclusiveModel(...)` at `TStructFunctionArray.cxx:87`
* make the `THnD` axes match the producer's grid -- they are hardcoded to the
  EG2 values and **must** agree or every lookup lands in the wrong cell
* make the input filename a parameter and allow a reload. It is currently
  hardcoded to `newphihist.root` in the cwd with `static` one-shot init, so EG2
  had to `cp`/`mv`/`rm` files around per target. An iteration loop doing that
  per pass per target is where mistakes will hide
* **decide the out-of-range policy** (below)

### 4. Closure test

Toy 5-D distribution with known A/Ac/Acc -> producer -> table -> HAPRAD.
Tier 0's four skips and Tier 3 should all go green. This validates the whole
chain *before* acceptance correction lands, which is exactly the "get this ready
before that" the professor asked for.

### 5. Blocked externally

Acceptance correction -> real amplitudes -> Tier 3 against real data ->
iteration.

### 6. Iteration loop

The existing chain does **zero** iterations; its only loop is over targets. The
papers require iteration but never give a count. From our own measured
`delta ~ 1.02-1.11`, model error is suppressed by `(delta - 1) ~ 0.05` per pass,
so 2 passes reach sub-percent and 3 are comfortably converged -- worse at high z
where we measured `delta ~ 0.86`. **Iterate to a tolerance, not a fixed count**,
and report how many it took; failure to converge is itself a diagnostic.

---

## Decisions needed

1. **Table or functional form?** A per-cell table cannot be evaluated outside
   the measured region, and a single Born point at `x = 0.234` has its tails
   sampling out to `x ~ 0.83, Q2 ~ 8.5` -- beyond RG-E statistics. Today that is
   silent edge-bin clamping, which is the failure mode documented in
   `~/externals/CLAS12_MIGRATION_PLAN.md`. A smooth fit in `(Q2, x, z, p_t)`
   extrapolates; precedent exists in `h3.f`/`h4.f` (6 parameters each) and in
   the dead `ConfigFile` hook asking for `par0..par4, A1..E1`.
   **Suggested compromise:** build the table first, then fit a smooth function
   *to the table* rather than to the raw 5-D data -- ordinary least squares on
   ~1500 points with error bars, not a 5-D likelihood.
2. **Event weights or per-bin centroids?** Centroids assume the RC is linear
   across a bin and couple the analysis binning to the RC evaluation. Weights
   decouple them, allow free rebinning, and are what makes iteration work --
   corrected events can be re-histogrammed onto whatever grid the next pass
   wants. HAPRAD is far too slow per event, so this means a precomputed grid in
   `(E, x, Q2, z, p_t, phi)` plus interpolation.
3. **Hadron species.** Everything is pi+ on a proton: `kMassDetectedHadron`, the
   fragmentation set, and the MAID grid (pi+ n). pi- needs its own exclusive
   grid; `haprad3` ships one.
4. **Nuclear targets (H15).** This is a free-proton calculation.
   `GetFactor3`'s `Z/A` on the exclusive tail is ad hoc. Does the RC cancel in
   your multiplicity ratios?

## Deliberately deferred

* **MAID2003 -> 2007.** Verified drop-in (same 18x47x61 layout), but swapping it
  now breaks Tier 2, which agrees with `haprad2` *because both use MAID2003*.
  Do it on both sides at once, after H1.
* **H11**, the `sigma_obs` form: C++ uses `sigma_B e^{delta_inf}(1 + delta_VR +
  delta_vac)`, the FORTRAN forces `delta_inf = 0` and uses
  `sigma_B (1 + alpha/pi delta)`. Equivalent to O(alpha^2). This sets the floor
  on Tier 3 agreement at ~1e-3; do not expect machine precision.
* `TQQTPhi`'s `+ep` on the upper tau limit where `qqtphi` uses `-ep`, and the
  missing-mass gate differing from `fhaprad`'s. Both are core calculation;
  recorded, not touched.

## Definition of done

```
./validation/compare.py --grid rge
```

Tier 0 all five PASS with no SKIPs; Tier 3 agreeing with HAPRAD 2.0 at ~1e-3.
H4 closes on its own, since real fitted amplitudes vanish as p_t -> 0.
