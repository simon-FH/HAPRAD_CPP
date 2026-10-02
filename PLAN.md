# HAPRAD_CPP for CLAS12 / RG-E — work plan

Companion to [H1_NOTES.md](H1_NOTES.md) and [validation/README.md](validation/README.md).

## Decisions taken

| | |
|---|---|
| Structure-function input | **iterative table** (`newphihist.root`), as originally designed. Fitting an equation for the amplitudes (A) or for the RC factor (B) comes *after* this works, using its outputs as training data |
| Applying the correction | **per-event weights**, interpolated from a precomputed RC grid — not one HAPRAD call per analysis-bin centroid |
| Hadrons | **pions**. pi+ first: every input exists. pi- once an exclusive grid exists (Phase 6) |
| Iterations | to a **tolerance**, not a fixed count. Expect 2-3: model error is suppressed by roughly `delta - 1 ~ 0.05` per pass at our measured `delta ~ 1.02-1.11` |

## Starting point

```
Tier 0 : pass     self-consistency; every model-independent quantity finite
Tier 1 : pass     kinematics vs HAPRAD 2.0 FORTRAN, ~4e-14, at 6 and 10.5 GeV
Tier 2 : pass     exclusive radiative tail vs FORTRAN, 3e-5 .. 1.5e-4
Tier 3 : FAIL     no structure-function table exists yet
```

The QED machinery is validated, and the A/Ac/Acc -> H1..H4 inversion round-trips
exactly. What is missing is the table, and the code around it to build, consume,
apply and iterate it.

**Measured cost:** one HAPRAD call at RG-E kinematics takes **0.45-1.5 s**. A
grid of 18,000 points is roughly 4-8 core-hours — overnight on a laptop with
`xargs -P`, or minutes as a cluster job array.

---

## Phase 0 — Foundations (local, no data needed)

### 0.1 Check that RG-E's `phi_PQ` is HAPRAD's `phi_h`

A convention mismatch here is silent: shifting phi by pi flips the sign of A_c
and nothing complains.

`phi_pq()` in `clas12-rge-analysis/src/rge_particle.c` rotates the virtual photon
onto +z with the scattered electron in the xz-plane, then takes
`atan2(p_y, p_x)` — the Trento construction, with phi = 0 on the lepton side.
HAPRAD writes `V_{1,2} = 2(a_{1,2} + b cos phi_h)` with `b = -p_t sqrt(lambda/lambda_q) < 0`,
which also puts phi_h = 0 on the lepton side. **So they very likely agree** —
but that rests on the sign convention of `rge_rotate_y`, which I have not
checked.

**Test:** for real events, pair each pion with its trigger electron by
`event_num`, compute `V1 = 2 k1·p_h` directly from the four-vectors, and compare
with `THadronKinematics::V1()` evaluated from `(x, Q2, z, p_T, phi_PQ)`.
Agreement settles it. Agreement only under `phi -> pi - phi` or `phi + pi` means
a conversion is needed.

**Done when:** V1 agrees to rounding on run 020026.

**DONE.** `validation/bin/phi_convention`, over 268,258 pi+ from run 020026:

```
hypothesis                   median          95%     frac < 1e-4
phi_h = phi_PQ            2.845e-06    1.215e-05         1.0000
phi_h = phi_PQ + pi       7.085e-01    4.223e+00         0.0001
```

The conventions agree; the residual is the ntuple's float storage. The test has
plenty of leverage -- shifting phi by pi moves V1 by 71% at the median -- so a
mismatch could not have hidden.

**Found along the way: the analysis needs a missing-mass cut.** 6,641 events
(2.4%) pass the DIS cuts but HAPRAD rejects them, *all* because
`p_x^2 < (M_p + m_pi)^2`. An independent M_x calculation from the ntuple flags
exactly the same 6,641. They form the exclusive peak at M_x ~ 0.94 GeV
(e p -> e' pi+ n) with its smeared tail, at mean z = 0.555. HAPRAD is right to
refuse them -- the semi-inclusive RC does not apply below the two-body threshold
-- and the analysis should cut them anyway. The same cut also disposes of the
838 backward-going pions (0.30%), which HAPRAD's p_t input would otherwise
silently treat as forward:

```
M_x > 1.2 GeV : keeps 96.1%, backward left 3
M_x > 1.4 GeV : keeps 93.1%, backward left 0
M_x > 1.5 GeV : keeps 91.0%, backward left 0
M_x > 1.6 GeV : keeps 88.2%, backward left 0
```

The exact value (1.4-1.6 GeV is typical, to clear the resonance region too) is
an analysis decision. Whatever it is, the producer and the RC grid must use the
same one, since it defines what "semi-inclusive" means for the exclusive tail.

### 0.2 Rework the consumer (`TSemiInclusiveModel`)

* **Read the binning from the file**, not from hardcoded `THnD` axes. The
  producer writes axis edges and variable names as metadata into the ROOT
  file, and the consumer builds its histograms from that. Producer and
  consumer then *cannot* disagree on the grid — which is otherwise the easiest
  way to make every lookup land in the wrong cell.
* **Filename as a parameter, with a reload.** Drop `static` one-shot init and
  the hardcoded `newphihist.root` in the cwd. This is what forced EG2's
  `cp`/`mv`/`rm` shuffle per target, and an iteration loop would repeat it
  every pass.
* **Explicit out-of-range and empty-cell handling, with counters.** Today it
  clamps to the edge bin silently. For v1 keep clamping, but *count* the
  fraction of integrand evaluations that fall outside the table or in an
  unfitted cell, and report it per call. That turns H7 from a guess into a
  measurement, and the decision about what to do there can wait for the number.
* **Un-comment the call** at `TStructFunctionArray.cxx:87`.

**Done when:** a synthetic table loads from an arbitrary path, the harness
reports the out-of-range fraction, and Tier 1/2 values are bit-identical.

**DONE.** `TSemiInclusiveModel` is now a class owned by `TRadCor`
(`LoadSemiInclusiveTable(path)`), and its file format is documented in
`TSemiInclusiveModel.h`: THnDs `A`, `Ac`, `Acc`, optional `fitted`, with each
axis *named* for its variable, so any four of `Q2, x, nu, z, pt, pt2` work in any
order. Malformed tables are refused with a specific message. `GetRC` takes `-m`
and refuses to run without a table; the `sed` hack is gone from
`exec_rad-corr_chain.sh`. With no table the amplitudes are exactly zero, and
all 108 model-independent values stayed bit-identical.

Two things were found that were not in the plan, and both shaped the result.

**Interpolation is required, not optional.** A nearest-cell lookup makes the
integrand piecewise constant, and the non-adaptive rule used for the inner R
integral (`TRV2LN`) cannot converge across the steps. With a toy table this gave
thousands of GSL tolerance failures per call and, at one point, an RC factor of
**2.24 that was pure integrator garbage**. A table with no steps gave zero
failures, which pins the cause. The reader now interpolates multilinearly
between cell centres:

```
toy table, x=0.3 Q2=4 z=0.5 pt=0.5   GSL failures   f1        lookups
  phi=150  nearest cell                     1533    2.2355    133k
  phi=150  interpolated                       14    1.0430     63k
  phi=30   nearest cell                    15096    0.9939    1.3M
  phi=30   interpolated                      446    0.9326     64k
```

For a constant table the two modes agree to 4e-15 (rounding). The remaining
failures come from the kinks multilinear interpolation leaves at cell centres;
their effect on the result is to be measured in Phase 2.3, where there is a
correct answer to compare against.

**Store the reduced amplitudes.** A raw table holds `Ac` and `Acc`, which stay
non-zero across the lowest p_t cell; the inversion then divides by p_t and
p_t^2 and the phi = 180 pole of H4 returns. The format therefore supports
`pt_scaling = "reduced"`: the table holds `Ac/p_t` and `Acc/p_t^2` and the reader
multiplies back by the p_t of each lookup, so the amplitudes vanish with the
right power by construction:

```
phi        f1 (raw)     f1 (reduced)
170        1.0362       1.0292
179.99     1.2267       1.0292
180      105.6460       1.0292
```

**The producer should write reduced tables.** With one, all five Tier 0 checks
pass -- none skipped -- and the table-coverage report shows 0.30% of lookups
outside the toy table over the RG-E grid.

### 0.3 Choose the model-grid variables

**Recommendation: bin the table in `(Q2, nu, z_h, Pt2)`, not `(Q2, x_B, z_h, p_t)`.**

RG-E's acceptance tool, `acc_corr`, already produces a 5-D acceptance ordered
`[Q2][nu][z_h][Pt2][phi_PQ]`. A table in the same variables and edges can use
that acceptance map directly. In `(x_B, p_t)` we would need a second
acceptance in different variables. The consumer converts trivially:
`nu = Q2 / (2 M x)`, `Pt2 = p_t^2`.

**Supported as of 0.2** -- the reader takes the variables from the axis names,
so this is now purely the producer's choice. `make_toy_table` already writes
`(Q2, nu, z, pt2)`.

---

## Phase 1 — The producer (develop on run 020026, run on the cluster)

`MakePhiTable`: RG-E ntuples -> `newphihist.root`. Starting point is
`git show b303487^:PhiHist/phihist.cpp`, deleted from this repo in `b303487`.

1. **Read the `DT` ntuples** — pi+, DIS cuts, vertex window, and a
   **missing-mass cut** (see 0.1 — without it 2.4% of events are exclusive
   and HAPRAD refuses them). Takes a file list,
   writes one small ROOT file. Built to run unchanged on the cluster.
2. **Binning from a config file**, never compiled in.
3. **Acceptance as a pluggable input** reading the `acc_corr` format. v0 runs
   with unit acceptance, *flagged in the output file*, so uncorrected tables
   cannot be mistaken for real ones.
4. **Per-cell fit** of `A + Ac cos(phi) + Acc cos(2 phi)`, keeping errors,
   chi^2 and a fit-status flag. The original's `MIN_BIN = 4` silently dropped
   sparse cells; record them as unfitted instead.
5. **Output** `AAcAcc_data` plus the binning metadata 0.2 reads.

`phi_PQ` is in radians in RG-E; the fit and HAPRAD both work in degrees.
Convert once, at the reader.

**Also from this phase:** the occupancy table for the binning study falls out
for free — the producer already fills every cell. Run it over the full dataset
on the cluster and choose the final grid from that, not from one run.

**Done when:** a table from run 020026 loads in HAPRAD and Tier 3 runs (it won't
*agree* yet — no acceptance, one run).

**DONE, except acceptance and per-target vertex windows (both blocked
externally).** `Utilities/bin/MakePhiTable`, configured by
`Utilities/config/rge_pip.cfg`:

```bash
MakePhiTable fill config.cfg counts_job17.root <ntuple files or globs>   # per cluster job
hadd counts.root counts_job*.root
MakePhiTable fit  config.cfg counts.root table.root
```

* **What A is, settled in practice.** K = sigma_Born / A does not follow any
  standard yield definition: scanning it, K goes roughly like 1/nu and 1/z with
  a non-power Q2 dependence. So the table *defines* A by it:
  `A = (fitted yield) x 2 M E nu / (cell volume) / K(cell centre)` -- the
  measured cross section in HAPRAD's variables, divided by K. Then K's formula
  cancels in the calculation. Tested: replacing line 99's `Sqrt(Q2 + y^2)` with
  `|q|` moved the 2.3 convergence results by at most 0.03 percentage points.
  This needs the conversion to use the real beam energy (2.3, Finding 2).
* **Bin-averaged fit.** A counted phi bin integrates the distribution, so
  fitting basis functions at bin centres attenuates the harmonics -- 1.1% for
  cos(phi), 4.5% for cos(2 phi) with 12 bins. The original `phihist.cpp` did
  this. The fit now uses bin-averaged basis functions.
* Writes `pt_scaling = reduced`, `interpolation = log` and
  `acceptance_corrected = no`; keeps fit errors, chi^2/ndf, event counts and the
  raw yield harmonics per cell; sparse cells are recorded as unfitted.
* **Merging is exact:** two fill jobs merged with `hadd` give a table identical
  to a single pass (0 of 14,406 cell values differ).

**On run 020026** (3.8 s to fill): 4,916,009 particles -> 275,737 pi+ after the
DIS cuts (the same count `phi_convention` found independently) -> 243,390 after
the vertex window and M_x > 1.5 GeV. On the full-dataset grid of the config
(138,240 cells) one run fits 1,120 cells; on a coarse 5x5x5x5 grid, 176.

With the coarse table loaded, **all five Tier 0 checks pass**, and the table
coverage gives the first real measurement of H7: over the RG-E test points the
tail integrals looked up structure functions **outside the table 4.9% of the
time, and in an unfitted cell 12.5%** (6.1% and 14.1% after the fit fix below:
the integrators adapt to the table, so the sampled points move). The full
dataset will fill empty cells; it will not reach outside the data's kinematic
range.

**Normalisation closure: DONE** (`validation/closure_normalisation.py`).
3M toy pions uniform in (Q2, nu, z, pt2, phi), each weighted by HAPRAD 2.0's
Born cross section times 1/(2 M E nu), through `fill` (new optional
`weight_branch`) and `fit`; compared cell by cell with `make_model_table`'s
table of the same model. The generator knows its luminosity, so the test is
**absolute**. Three comparisons per cell, medians over cells:

| grid (x 12 phi) | cells | producer / model, A | bin-centring, A | producer / cell average: A, Ac/A, Acc/A |
|---|---|---|---|---|
| 4^4 | 223 | 1.116 | 1.115 | 1.000, 1.002, 0.997 |
| 6^4 | 792 | 1.060 | 1.061 | 1.000, 0.994, 1.000 |
| 8^4 | 1601 | 1.036 | 1.036 | 1.000, 1.003, 1.007 |

* **MakePhiTable is right, absolutely.** Against the true cell average it gives
  A to 0.01%, and the modulations to within their fit errors, flat in every
  variable -- so the `2 M E nu` Jacobian, the cell volume, the luminosity and
  the phi fit are all correct. A missing nu factor would have been a factor ~3
  across the nu range.
* **Everything left is bin-centring:** the table holds cell *averages* and
  HAPRAD reads them as values at the cell *centre*. For A it is +11.5%, +6.1%,
  +3.6% on the three grids, shrinking roughly as (cell width)^2; the config
  grid is finer still. For the modulations it is ~1%, except in the lowest pt2
  bin: there Ac/A is low by ~14% on *every* grid, because Ac ~ p_t is not
  smooth at p_t = 0 and halving the bin does not help. A bin-centring
  correction (or storing the cell's mean kinematics) is a candidate for Phase 4,
  where the iteration can absorb it.

Three things this check found, now fixed in `MakePhiTable`:

1. **The fit was biased low by ~1/n** for n events per phi bin (0.1%, 0.5%,
   1.3% on the three grids). It took each bin's variance from its own content,
   so bins that fluctuated low got small errors and pulled the fit down -- on the
   full grid's sparse cells, several percent, varying cell to cell. The fit now
   takes the variances from its own prediction and iterates (quasi-Poisson).
   The producer closure above is unchanged: six seeds, pulls of Ac/A +0.004 and
   Acc/A +0.020, widths 0.93-1.05.
2. **The table's absolute scale is NOT arbitrary.** sigma_Born and the
   inelastic tail scale with the table; the exclusive tail (`tai[1]`, MAID) is
   added in absolute units and does not. So the scale sets the exclusive tail's
   share of the RC factor. The producer's comment claiming the constant
   "cancels in every RC factor" was wrong. **A table built without a luminosity
   effectively switches the exclusive tail off.** `MakePhiTable` now takes an
   optional `luminosity` (events per nb; default 1 = arbitrary units, recorded
   as such in the table). For real data this needs the acceptance correction,
   the integrated luminosity, and a decision on per-nucleon normalisation for
   the nuclear targets (MAID is a free-proton model).
3. **phi in the cell volume is now in radians.** HAPRAD's sigma is per radian
   of phi_h: at x = 0.2, Q2 = 2, z = 0.5, pt2 = 0.2 its B0 is 36 nb/GeV^2,
   against a rough leading-order estimate of ~11 per radian and ~0.2 per degree.
   Degrees would have put a factor 57 into the absolute scale.

---

## Phase 2 — Closure tests (local, no data needed)

The point of this phase is to prove the chain is right *before* acceptance
correction and the full dataset exist.

1. **Consumer closure.** Write a synthetic table with smooth, known amplitudes
   that carry the physical p_t scaling (`Ac ~ p_t`, `Acc ~ p_t^2`). Run HAPRAD,
   project sigma_Born(phi) back into harmonics, recover the inputs. Extends the
   round-trip check already done, but through the real file-loading path.
   Tier 0's four skipped checks should now run and pass, and the phi = 180 pole
   (H4) should be gone.
2. **Producer closure.** Generate toy events in the `DT` format with a known
   azimuthal modulation, run `MakePhiTable`, recover the amplitudes within
   their fitted errors.
   **DONE** (`validation/closure_producer.py`). 187 cells, six seeds: pulls of
   Ac/A -0.031 +- 0.030 and Acc/A +0.025 +- 0.030, widths 0.93-1.05. Before the
   bin-averaged fit the two had opposite-sign offsets of about 0.1.
3. **A Tier 3 that can actually agree.** Write a Born-only FORTRAN driver
   (`conkin` + `bornin`; no tail integrals, so it is fast), use it to tabulate
   A/Ac/Acc *from HAPRAD 2.0's own PDF x FF model* on our grid, and feed that
   table to the C++. Same model both sides, so Tier 3 compares like with like.
   Expect agreement at the level of the table's discretisation, with a floor
   from H11 (measured below at 0.2-0.4%).
4. Add all three to `validation/compare.py`.

**Step 1: done in 0.2.** With a reduced toy table all five Tier 0 checks pass.

**Step 3: done, with three findings -- one of which needs a decision.**

Tools (`validation/README.md`): `born_harmonics`, HAPRAD 2.0's Born alone,
matching the full calculation to 6e-13; `make_model_table`, which turns its
harmonics into a table using the C++'s own kinematic factor and checks it at
every cell centre (agreement to <=1e-5, set by the inversion dropping the
`-m_e^2 lambda_q` term of lambda); `table_convergence.py`, Tier 3 against grid
resolution; and `TRadCor::CalculateBorn`, sigma_Born without the tails,
bit-identical to the full calculation's.

*Finding 1 -- interpolate log(A), not A.* A carries a steep, roughly 1/Q^4,
cross-section-like dependence. Interpolated linearly, sigma_Born at off-centre
points was off by up to 53% on practical grids; refining Q2 helped most,
refining p_t^2 alone not at all. The reader now accepts
`interpolation = "log"`: log(A) and the ratios Ac/A, Acc/A, which are smooth.
That cut the base-grid error from 53% to 14%, and to 0.7% on finer grids.
**The producer should write log tables.**

*Finding 2 -- the tail carries a ~2% bias that refinement cannot remove.* With
the same physics on both sides, the RC factor converged to 2.2-2.5% and stayed
there. Splitting f1 = sigma_obs/sigma_Born + tail/sigma_Born:

```
x=0.234 Q2=2.5 z=0.34 phi=60    delta part -0.37%   tail part -1.25%
x=0.234 Q2=2.5 z=0.34 phi=180   delta part -0.34%   tail part -1.94%
x=0.3   Q2=4   z=0.5  phi=30    delta part -0.42%   tail part -0.72%
x=0.234 Q2=2.5 z=0.7  phi=60    delta part -0.30%   tail part -0.58%
```

The delta part is H11 and does not depend on the table. The tail part has the
same sign at every point. The cause is `TStructFunctionArray.cxx:79-82`: at
shifted kinematics the yield-to-structure-function inversion does not use the
beam energy E but reconstructs an "effective" one,

```cpp
Double_t a = S/(2M) * (S/(2M) - nu_Born) * tldQ2 / Q2;    // depends on the Born point
Double_t tldE = 0.5 * (tldNu + Sqrt(SQ(tldNu) + 4 * a));
```

so one table cell yields different structure functions depending on which Born
point is asking. At Born kinematics it reduces to E exactly, which is why every
Born-level check passed. Line 93 of the same expression (`tldXg`) uses the real
E. In a scratch copy with `tldE = E`, the tail part falls to +0.04%, -0.19%,
+0.03%, -0.17% and the RC factor converges with refinement to 0.6% -- close to
the H11 floor -- instead of levelling off at 2.2-2.5%:

```
grid (Q2 x nu x z x pt2)   f1 error, current   f1 error, tldE = E
 9x 8x10x24                    2.59%               0.58%
18x16x20x24                    2.29%               0.95%
27x24x30x36 (wider)            2.49%               0.60%
```

**This is a change to the core calculation, so it is not made.** It needs a
decision -- ideally from the port's author, together with the next point.

*Finding 3 -- one test point exposes the reference model, not the C++.* At
x = 0.15, Q2 = 2, z = 0.3 the FORTRAN's cos(2 phi) term sits exactly on the
authors' cap: `2 Bcc/B0 = 0.9000`, as `h4.f` runs away at low x. Under the cap
H4 goes like 1/p_t^2 and depends on y, so it is not a function of the table
variables and no table can follow it -- its error is large on a coarse p_t^2
grid (20% in sigma_Born) and only shrinks with refinement. A 45% cos(2 phi)
modulation is not physical; real data will not have it.

**Questions for the port's author:**
1. Should the inversion use the beam energy E (`tldE = E`) rather than the
   reconstructed effective energy? The evidence above says yes.
2. `TStructFunctionArray.cxx:99`: `N = Q^4 * Sqrt(tldQ2 + SQ(tldY)) / SQ(tldY)`
   adds Q^2 (GeV^2) to y^2 (dimensionless). Line 61 already computes
   `tld_sq = Sqrt(tldQ2 + SQ(tldNu)) = |q|`; was `tldNu` meant? *No longer
   blocking:* Phase 1 defines A so that this formula cancels, and swapping it
   for `|q|` was tested to make no difference. Still worth an answer, since the
   code should say what it means.

**Resolution guidance, from the converging runs:** with log interpolation,
~24 bins in p_t^2 and ~18 x 16 x 20 in (Q2, nu, z) the table contributes below
1% to the RC factor once Finding 2 is resolved.

**Done when:** all five Tier 0 checks pass with no SKIPs, and Tier 3 agrees with
HAPRAD 2.0 on the model-derived table.

The Born-only driver from step 3 is also the tool that later generates training
data for the equation fits.

---

## Phase 3 — Applying the correction as weights

1. **`MakeRCGrid`** — run HAPRAD over a grid in `(x, Q2, z, p_t, phi)`, one grid
   per beam energy (RG-E has three: 10.3894, 10.4057, 10.5473 GeV). Written as
   independent chunks so it parallelises with `xargs -P` locally or as a job
   array on the cluster.
2. **`ApplyRC`** — read events, interpolate the RC factor at each event's own
   kinematics, write the weight as a friend tree. Any analysis binning can then
   be applied afterwards.
3. **Interpolation error budget.** Call HAPRAD directly at random off-grid
   points and compare with the interpolated value. The grid is fine enough when
   this sits well below the statistical errors.

**Done when:** weights exist for run 020026 and the interpolation error is
quantified.

---

## Phase 4 — The iteration loop

```
table_k  ->  RC grid_k  ->  weights_k  ->  table_{k+1} from RC-corrected data
```

* Scripted end to end, **per target** — fitting the amplitudes on each target's
  own data absorbs the nuclear structure-function differences automatically,
  which is why EG2 produced a table per target.
* Iteration 0's table comes from RC-uncorrected (but acceptance-corrected) data.
* Stop when `max |delta_k - delta_{k-1}|` over the grid falls below a tolerance.
  Record the number of iterations; failure to converge is itself a sign that the
  grid or the cuts are wrong.

**Done when:** the loop converges on run 020026.

---

## Phase 5 — Production

Full dataset on the cluster, all targets, real acceptance. `GetRC` and the
centroid path are superseded by the weights and can be retired.

---

## Phase 6 — After it works

* **Equation fit for the amplitudes (A)** — the converged tables are the
  training data. The hard part is a functional form valid over the *integration
  domain*, not just the measured one: the existing `h4.f` form reaches 3.7e53 at
  x = 0.015, so EG2-era shapes do not transfer.
* **Equation fit for the RC factor (B)** — the RC grids are the training data.
  Fit `delta_0`, `delta_c`, `delta_cc` in five variables. Refit every iteration.
* **pi-** — needs a MAID gamma* n -> pi- p exclusive grid. None exists on this
  machine (`haprad3` has pi+ and pi0 only).
* **MAID2003 -> 2007** — a verified drop-in, but swap it on *both* sides at once,
  or Tier 2 loses the reference it currently agrees with.

---

## Blocked on things outside this repo

| | needed for | notes |
|---|---|---|
| Acceptance maps | Phase 1 onward, for real numbers | `acc_corr` in `clas12-rge-analysis` produces them; needs simulation |
| Full dataset | Phase 1 binning choice, Phase 5 | `runs_all.txt` lists 384 runs; 1 is processed locally |
| Target separation | Phase 1 | RG-E constants hold only a global `vz` window (-40, 26.12 cm); per-target windows for the LD2 cell vs the solid foil are still to be defined |
| Luminosity | absolute table scale, hence the exclusive tail (Phase 1) | integrated charge and target areal density per run; for nuclear targets, also a per-nucleon convention to match MAID's free proton |
| Nuclear corrections beyond the SFs | Phase 5 | the exclusive tail (MAID is free-proton), Coulomb distortion (`RGE_RC_CC` already exists), and external radiation in the target. Whether they matter depends on whether the observable is a ratio in which they cancel — a call for the professor |

## Order of work

Phase 0 and Phase 2 need no data at all and should come first — they take the
chain from "validated pieces" to "validated whole". Phase 1 can be developed on
run 020026 in parallel and moved to the cluster unchanged. Phases 3-4 follow
once a table loads.
