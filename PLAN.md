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
3. **A Tier 3 that can actually agree.** Write a Born-only FORTRAN driver
   (`conkin` + `bornin`; no tail integrals, so it is fast), use it to tabulate
   A/Ac/Acc *from HAPRAD 2.0's own PDF x FF model* on our grid, and feed that
   table to the C++. Same model both sides, so Tier 3 compares like with like.
   Expect agreement at the level of the table's discretisation, with a floor
   near 1e-3 from H11.
4. Add all three to `validation/compare.py`.

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
| Nuclear corrections beyond the SFs | Phase 5 | the exclusive tail (MAID is free-proton), Coulomb distortion (`RGE_RC_CC` already exists), and external radiation in the target. Whether they matter depends on whether the observable is a ratio in which they cancel — a call for the professor |

## Order of work

Phase 0 and Phase 2 need no data at all and should come first — they take the
chain from "validated pieces" to "validated whole". Phase 1 can be developed on
run 020026 in parallel and moved to the cluster unchanged. Phases 3-4 follow
once a table loads.
