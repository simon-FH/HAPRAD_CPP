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
`THnD` histograms. **That file does not exist anywhere on this machine.**

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

## 5. Route A — use the FORTRAN model (recommended first step)

`haprad2/semi_inclusive_model.f` computes `H1..H4` from GRV94-LO PDFs x PKH
fragmentation functions x a Gaussian k_T, plus the empirical `h3.f`/`h4.f` fits
for the cos(phi) and cos(2 phi) terms. It is regular as p_t -> 0 (no inversion,
therefore no reciprocals), which removes H4's pole as a side effect.

**This is now buildable.** Everything it needs is already present:

* `PDFSET`/`STRUCTM` — `~/cernlib/lib/libpdflib804.so`, built in commit `4e67f12`
* `PKHFF` — `pkhff.f`, already compiled into `libTRadCor.so`
* `init_pdf.f` — already provides `exec_structm` / `exec_pkhff` wrappers
* `h3.f`, `h4.f`, `partons.inc`, `constants8.inc` — in `haprad2/`

The change is surgical. The conversion from the Mulders basis `H1z..H4z` to the
Akushevich basis `h1..h4` in `TStructFunctionArray` is **already a faithful port**
— I checked it line by line against `strf()` in `ihaprad.f`, including the `aa`
term. Only the *source* of `H1z..H4z` differs. So the work is: compile
`semi_inclusive_model.f`, `h3.f`, `h4.f` into the library, declare

```cpp
extern "C" void semi_inclusive_model_(double* q2, double* x, double* y, double* z,
                                      double* pt2, double* mx2, double* pl,
                                      double* H1z, double* H2z, double* H3z, double* H4z);
```

and replace the A/Ac/Acc block with that call.

Linking the reference implementation verbatim rather than re-porting it has a
large advantage: Tier 3 becomes a true apples-to-apples comparison — same model
both sides — so any residual disagreement is a genuine defect rather than a
modelling difference.

### Landmines

1. **`semi_inclusive_model.f:178` is an unconditional `stop`:**

   ```fortran
   if (abs(4.d0 * m_cos_phi) .gt. 0.9d0) stop      ! <-- aborts
   if (abs(4.d0 * m_cos_phi) .gt. 0.9d0) then      ! <-- cap, unreachable
      H3m = 0.9d0 * sign(1.d0,m_cos_phi) * ...
   ```

   A debugging leftover that makes the intended cap dead code and kills the
   process instead. It must be removed before linking, or any RG-E point where
   the modulation exceeds the bound takes the whole job down. (This is presumably
   part of why Klimenko's group reported making the code "return zero XSEC
   instead of crashing".)

2. **`nc`-counted one-time init.** `init_pdf` and the FF tables are initialised
   on the first call via a saved counter, and `IFINI` is passed through
   `COMMON /FRAGINI/`. Fine single-threaded; do not call it from threads.

3. **Data files by relative path.** `pkhff.f` opens `plo.grid`, `pnlo.grid`,
   `klo.grid`, `knlo.grid`, `hlo.grid`, `hnlo.grid` from the cwd, as
   `exclusive_model` does with `pi_n_maid.dat`. `validation/data/` already
   symlinks all of these.

4. **Hard-wired to pi+ on a proton.** `ISET=1, ICHARGE=1` in the FF call.

### Coverage at RG-E — the thing to check before trusting it

| ingredient | validity | RG-E (pi+, DIS cuts) |
|---|---|---|
| GRV94 LO | `SCALE = max(sqrt(Q2), 1)` | Q2 up to 8.5 |
| Gaussian width `sgmpt(x,z)` | clamped to [0.02, 0.15] | fitted to EG2-era data |
| `mx2 < (mp+mpi)^2` | returns zero | fires near threshold |
| `h3`/`h4` fits | `Q0=1`, `lambda=0.25` | `(log(Q2/L^2)/log(Q0/L^2))^(bb/x)` grows fast at small x |

`h4.f` in particular raises a logarithm to the power `bb/x` with `bb = 6.88`; at
`x = 0.015` that exponent is ~460. Worth plotting over the RG-E range before
trusting it.

**The lesson from `~/externals/CLAS12_MIGRATION_PLAN.md` applies directly here:**
*a replacement model must be valid over the integration domain, not just over the
kinematic points you ask for.* The tail integrals sample shifted kinematics well
outside the Born point.

## 6. Route B — data-driven amplitudes (the iteration step)

Both papers insist the RC procedure must be **iterative**: fit the corrected
data, feed that fit back as the model, repeat. Route B is that step, and it is
what the `Utilities` chain was built for. It should be layered on top of a
working Route A, not instead of it.

Requirements:

* Rebuild the `newphihist.root` equivalent from RG-E, binned in `(Q2, x, z, p_t)`.
  The existing `THnD` binning — `Q2 in [1,4]`, `x in [0.1,0.55]`, 6x5x10x5 bins —
  is far too small for RG-E, where `Q2` reaches 8.5 and `x` falls to 0.015.
  Out-of-range values are currently **clamped to the edge bin**, silently.
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
  agreement at the 1e-3 level, not machine precision — the C++ forms
  `sigma_B * exp(delta_inf) * (1 + delta_VR + delta_vac)` where the FORTRAN
  forces `delta_inf = 0` and uses `sigma_B * (1 + alpha/pi * delta)` (hurdle H11).
  Equivalent to O(alpha^2), but not digit-for-digit.
* H4 closes on its own: the `phi=180` spike disappears when the amplitudes carry
  their proper p_t dependence.

## 8. Open questions

1. **Route A or straight to B?** A is a day's work, unblocks every check, and
   gives a validated baseline. B is what the physics ultimately needs. My
   recommendation is A first, precisely so B has something to be checked against.
2. **Which hadron?** The model is pi+/proton throughout. RG-E wants pi-, and the
   exclusive grid would need `pim`. `haprad3` ships one.
3. **Nuclear targets.** This is a free-proton calculation; `GetFactor3`'s
   `Z/A` weighting of the exclusive tail is ad hoc. Does the RC cancel in your
   multiplicity ratios? That decision belongs to you, not to the code.
