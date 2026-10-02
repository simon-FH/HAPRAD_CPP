# H1 — the missing semi-inclusive structure-function model

Written after closing H2 (CERNLIB), H12 (validation harness), the p_t Jacobian
bug, and H4. H1 is now the only thing blocking every remaining check.

---

## 1. What is actually missing

`TStructFunctionArray::Evaluate()` is supposed to supply the four semi-inclusive
structure functions `H1..H4` at the photon-shifted kinematics. It does not.

```cpp
// TStructFunctionArray.cxx:83
Double_t A;
Double_t Ac;
Double_t Acc;
// HapradUtils::SemiInclusiveModel(tldQ2, tldX, tldY, tldZ, tldPt2, tldPx2, tldPl, A, Ac, Acc);
```

Three uninitialised doubles, and the call that was meant to fill them commented
out. Everything downstream reads stack garbage. The only way the code has ever
produced numbers is the `sed` hack in `Utilities/exec_rad-corr_chain.sh`, which
substitutes three fitted **constants** into the source, recompiles, runs, and
reverts — once per target.

`TSemiInclusiveModel.cxx` still contains the intended implementation: it loads
`newphihist.root`, an `NTuple` of `(Q2, Xb, Zh, Pt, A, Ac, Acc)`, into three 4-D
`THnD` histograms. That file is absent — but it is **not** a missing asset to go
looking for. Those are fitted azimuthal-modulation amplitudes, so it is an
*output of the analysis*, produced from your own data. The real statement is
that **nothing in the chain currently produces it**; see §6.

## 2. What this blocks

| | status |
|---|---|
| Tier 3 (cross sections vs HAPRAD 2.0) | cannot run — `sigma_Born` is garbage |
| Tier 0 phi symmetry | SKIP |
| Tier 0 `sigma_Born` depends on z | SKIP |
| Tier 0 `sigma_Born` depends on p_t | SKIP |
| Tier 0 no spikes in phi | SKIP |
| H4 (phi=180 pole) | symptom of H1 — see §4 |

Tier 1 (kinematics) and Tier 2 (exclusive tail) pass and are unaffected: neither
touches this model.

## 3. Why the constants are worse than they look

With `A`, `Ac`, `Acc` held constant the Born cross section loses **all** hadronic
dependence. Measured directly:

```
z=0.34 pt=0.424 -> sigb=85670.4     z=0.70 pt=0.30 -> sigb=85670.4
z=0.50 pt=0.30  -> sigb=85670.4     z=0.90 pt=0.30 -> sigb=85670.4
```

Identical to six figures across the whole hadronic phase space. This is
structural, not a tuning problem: every `h_i` in the conversion carries a leading
`tldZ` that cancels the `1/tldZ` in `coef`, so constant amplitudes give a
z-independent `sigma_Born` by construction.

Meanwhile the *tail* integrals do see z and p_t, so the RC factor is a ratio of
a shape-carrying numerator to a shape-free denominator. Akushevich–Shumeiko–Soroko
(1999) state plainly that RC depends strongly on the p_t^2 slope; a flat model is
precisely the wrong input.

## 4. The p_t singularity, and its link to H4

This model does not evaluate structure functions — it **inverts** measured
azimuthal amplitudes. In the Born cross section the modulations enter as

```
sigma_0 ~ A + cos(phi) * A_c + cos(2 phi) * A_cc,   A_c ~ b ~ p_t,  A_cc ~ b^2 ~ p_t^2
```

so solving for `H3`/`H4` necessarily divides it back out:

```cpp
H3z = ... / Sqrt(tldPt2);
H4z = ... / tldPt2;
```

That is harmless **only if the amplitudes vanish with the right power as
p_t -> 0**. Constants do not, so `H4z` diverges. This is exactly what H4 turned
out to be: at `phi_h = 180 deg` the shifted p_t passes through zero inside the
integration domain, `tldPt2 * H4z` becomes `0 * inf`, and NaN poisoned the whole
inelastic tail. Measured:

```
flat amplitudes            A_c ~ p_t, A_cc ~ p_t^2
phi=179     1.978841e+04   1.654784e+04
phi=179.9   2.679992e+04   1.654917e+04
phi=179.999 8.555769e+05   1.654919e+04
phi=180     5.327075e+07   1.654919e+04
GSL warnings: 43586        GSL warnings: 0
```

**Consequence for the data-driven route.** `newphihist.root` was binned in p_t
(5 bins over [0,1]), so the design did carry p_t dependence — the `sed` reduction
to constants is what destroyed it. But a *piecewise-constant* p_t dependence is
not enough either: the lowest bin is flat over `p_t in [0, 0.2]`, so `A_cc(0) != 0`
and `H4z ~ 1/p_t^2` still diverges, only less violently. Any data-driven route
must enforce the scaling `A_c ~ p_t`, `A_cc ~ p_t^2` in the inversion itself.

Note the `<cos phi>` / `<cos 2phi>` caps in `semi_inclusive_model.f` do **not**
rescue this: capping bounds the modulation but leaves `H4z ~ 1/p_t^2` at the cap.

## 5. The machinery is correct — verified

Confirmed with the author of the C++ port: this is a deliberately data-driven
design, with **no external model**. `TSemiInclusiveModel` exists to get the
parameters out of our own data. Three measurables constrain three quantities,
with `H1` and `H2` tied together because only three are independent:

```cpp
// TStructFunctionArray.cxx:90-100        FORTRAN: semi_inclusive_model.f:151
Double_t rlt = 0.14;                   // DATA rlt/0.14d0/
RelH1H2 = (1 + 4*SQ(M*tldX)/tldQ2) / 2 / tldX / (1 + rlt);
H1z = H2z * RelH1H2;                   // H1 = H2/(2x(1+rlt)) * (1+4mp^2x^2/q2)
```

That is the sigma_L/sigma_T closure, identical to the FORTRAN. `A` then fixes
`H2` (and `H1` through it), `A_c` fixes `H3`, `A_cc` fixes `H4`.

**The inversion round-trips exactly.** Injecting known amplitudes, computing
`sigma_Born(phi)` at twelve angles and projecting out the harmonics returns what
went in:

```
              injected      recovered
   Ac/A       -0.137000     -0.137000
   Acc/A      +0.041000     +0.041000
```

So the whole chain — the `RelH1H2` closure, the `H3z`/`H4z` formulas, the
Mulders-to-Akushevich conversion, and `TBorn`'s theta_B contraction — is
self-consistent. **Nothing here needs porting or replacing. The only thing
missing is the input.**

A useful consequence: a global rescaling of `A` cancels in the RC factor, since
`sigma_Born` and both tails scale together. **The fit does not need absolute
normalisation or luminosity** — only the correct *shape* in
`(Q2, x_B, z_h, p_t^2)` and the correct harmonic ratios.

An earlier version of these notes recommended linking the FORTRAN PDF x FF model
instead. That was based on a misreading of the design and is withdrawn. PDFLIB
remains built and available, and is still what makes `haprad2` runnable as the
Tier 1/2/3 reference, which is its own justification.

## 6. The actual work: producing the amplitudes

The pipeline, end to end:

1. Bin the data in `(Q2, x_B, z_h, p_t^2)`.
2. In each bin, histogram `phi_PQ` — the azimuthal angle about the virtual
   photon.
3. Acceptance-correct that histogram.
4. Fit `A + A_c cos(phi) + A_cc cos(2 phi)`; the three parameters are this bin's
   amplitudes.
5. Write `(Q2, x_B, z_h, p_t, A, A_c, A_cc)` per bin — that is
   `newphihist.root`.
6. `TSemiInclusiveModel` loads it; `TStructFunctionArray` looks up and inverts
   to `H1..H4` (§5).
7. HAPRAD computes the RC.
8. Apply the RC to the data, refit, repeat — the iteration both papers require.

### What the chain actually produces today

`newphihist.root` is an analysis product, so the question is why the chain does
not emit one. Because `FitPhiPQ` fits the **integrated** sample:

```cpp
// FitPhiPQ.cxx:66 -- CutPID && CutDIS && CutVertex, and no bin cut
dataChain->Draw(Form("PhiPQ>>data(%i, -180., 180.)", NbinsPhiPQ), ...);
```

One fit over everything, giving three numbers per target. Those three numbers are
what `exec_rad-corr_chain.sh` seds into the source, and they are the origin of
the flat model in §3. Nothing is binned in `(Q2, x, z, p_t)` at any point.

The missing piece is small, because `GetCentroids` **already builds exactly the
histograms required** — acceptance-corrected `phi_PQ` per 5-D bin — and then
throws the shape away:

```cpp
// GetCentroids.cxx:264
histPhiPQ_Corr[i]->Divide(histPhiPQ_Data[i], histPhiPQ_Acceptance[i], 1, 1);
meanPhiPQ[i] = histPhiPQ_Corr[i]->GetMean();   // only the mean is kept
```

So the work is to fit `A + Ac cos(phi) + Acc cos(2 phi)` to `histPhiPQ_Corr[i]`
in each bin and write the amplitudes out, rather than only their mean.

### Why a data table cannot stand alone

This is the part that decides the architecture. The tail integrals evaluate the
structure functions at **shifted** kinematics — `Q2 + R*tau`, `W2 - R(1+tau)`,
and the corresponding shifted `z` and `p_t` — which range far outside the region
where you have acceptance, and by construction cover the collinear peaks. A table
built from data has nothing to say there, and the present code handles that by
**silently clamping to the edge bin**:

```cpp
// TSemiInclusiveModel.cxx
if (bin[i] > xmaxs[i]) bin[i] = xmaxs[i] - 0.001;
if (bin[i] < xmins[i]) bin[i] = xmins[i] + 0.001;
```

which is precisely the silent-wrong-answer mode documented in
`~/externals/CLAS12_MIGRATION_PLAN.md`.

So model and data are not alternatives; they are layers. The model has to be
defined over the whole integration domain, with the data constraining it where
statistics exist. **This is why Route A is a prerequisite for Route B rather than
a competitor to it** — and it is also the natural seed for iteration 0, since the
first RC has to be computed before there is any corrected data to fit.

Requirements:

* Add a per-bin fit step (above), writing `(Q2, x, z, p_t, A, Ac, Acc)`.
* Widen the binning. The `THnD` axes — `Q2 in [1,4]`, `x in [0.1,0.55]`,
  6x5x10x5 bins — are far too small for RG-E, where `Q2` reaches 8.5 and `x`
  falls to 0.015. `Binning.hxx` is currently a 2x2x2x2 placeholder with 5 phi
  points per cell, which is 2 degrees of freedom against 3 fit parameters.
* Use the data only inside its support; fall back to the model outside it,
  rather than clamping.
* Enforce `A_c ~ p_t`, `A_cc ~ p_t^2` in the inversion (§4), otherwise H4 returns.
* Note `TSemiInclusiveModel.cxx` uses `THnD::Fill(values, weight)`, which
  *accumulates*. That is only correct if the tuple holds exactly one entry per
  bin; otherwise it needs averaging.
* Delete the `sed`-into-source mechanism in `exec_rad-corr_chain.sh` entirely.

## 7. Definition of done

Measurable, no judgement calls:

```
./validation/compare.py --grid rge
```

* Tier 0: all five checks PASS, none SKIP
* Tier 3: `sib`, `sig_obs`, `tai_in`, `f1..f3` agree with HAPRAD 2.0. Expect
  agreement at the few-1e-3 level, not machine precision — the C++ forms
  `sigma_B * exp(delta_inf) * (1 + delta_VR + delta_vac)` where the FORTRAN
  forces `delta_inf = 0` and uses `sigma_B * (1 + alpha/pi * delta)` (hurdle H11).
  Equivalent to O(alpha^2), but not digit-for-digit: measured at 0.2-0.4% of
  the RC factor in PLAN.md 2.3 (an earlier estimate here said ~1e-3).
* H4 closes on its own: the `phi=180` spike disappears when the amplitudes carry
  their proper p_t dependence.

## 8. Open questions

1. **Coverage outside the measured region — the open one.** A single Born point
   at `(Q2=2.5, x=0.234, z=0.34, p_t^2=0.18)` has its tails sample
   `Q2 1.06..8.46`, `x 0.234..0.830`, `z 0.34..0.83`, `p_t^2 0.18..0.48`.
   RG-E pi+ statistics reach `x ~ 0.695`, so the amplitudes are needed where
   there is no data, and high-x/high-z is exactly where the statistics thin out.
   Today that is handled by silently clamping to the edge bin. Does "fit
   `func(Q2,xB,Zh,Pt2)`" mean a per-bin table, or a smooth functional form in
   those four variables that can be evaluated outside the measured region? The
   second would resolve this cleanly. Worth settling before building the fitter.
2. **Which hadron?** The model is pi+/proton throughout. RG-E wants pi-, and the
   exclusive grid would need `pim`. `haprad3` ships one.
3. **Nuclear targets.** This is a free-proton calculation; `GetFactor3`'s
   `Z/A` weighting of the exclusive tail is ad hoc. Does the RC cancel in your
   multiplicity ratios? That decision belongs to you, not to the code.
