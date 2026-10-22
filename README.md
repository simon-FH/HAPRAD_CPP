HAPRAD_CPP
==========

C++ port of HAPRAD 2.0, the FORTRAN code of Akushevich et al. for the
lowest-order QED radiative corrections (RC) to semi-inclusive hadron
leptoproduction, e N → e′ h X. Here it is used for CLAS12 RG-E π⁺.

## What it computes

For each kinematic point (x, Q², z, pₜ, φₕ), the ratio of the observed to the
Born cross section:

    η = σ_obs / σ_B = r_vr + r_in + (Z/A) r_ex,      event weight w = 1/η

* `r_vr` = e^δinf (1 + δ_VR + δ_vac): soft-photon emission, the vertex
  correction and vacuum polarisation, which multiply the Born cross section.
* `r_in`: the inelastic radiative tail, hard-photon emission integrated over
  the photon's phase space.
* `r_ex`: the exclusive radiative tail from e p → e′ π⁺ n γ, built from MAID2007.

The hadronic input are the semi-inclusive structure functions. The library
reads them from a table of the azimuthal amplitudes A, Ac, Acc of the Born
cross section in (Q², ν, z, pₜ²) cells, interpolated between cell centres.
The table is measured from the data, so the correction is iterated: the table
is refilled with RC-weighted events until the weights stop changing.

## The chain

1. `MakePhiTable fill` / `fit`: count selected events in (Q², ν, z, pₜ², φ) and
   fit A + Ac cos φ + Acc cos 2φ in each cell, giving the table.
2. `MakeRCGrid`: run HAPRAD at the nodes of a fixed 5-D grid with that table.
3. `ApplyRC`: interpolate η from the grid to every event; the weights go to a
   friend tree next to each ntuple.
4. Refill the table with the weights and repeat: `Utilities/rc_iterate.py` on
   one machine, `Utilities/farm/` as a swif2 workflow on the JLab farm.

## Layout

| Directory | Contents |
|---|---|
| `.` | the library `libTRadCor` (`TRadCor` is the entry point) |
| `Utilities/` | the RG-E tools above, their config (`config/rge_pip.cfg`) and the farm workflow; see [Utilities/README.md](Utilities/README.md) |
| `haprad2/` | the original HAPRAD 2.0 FORTRAN (the reference), the fragmentation-function grids and the MAID grids |
| `validation/` | the C++ against the FORTRAN, closure tests on toy data, and the figures; see [validation/README.md](validation/README.md) |

## Build

Needs ROOT built with MathMore (GSL) and `gfortran`.

```bash
make                # slib/libTRadCor.so
make -C Utilities   # the tools
```

On the JLab farm, `module load root/6.28.10-gcc11.4.0` (the 6.30 builds there
lack MathMore). HAPRAD reads `pi_n_maid.dat` from the working directory; the
scripts run in `haprad2/` or link the file.

CERNLIB is not needed for the library. The FORTRAN reference in `validation/`
needs PDFLIB 8.04 (`export CERN_LIB=/path/to/cernlib/lib`), which still builds
standalone from the community CERNLIB sources.

## Status

* The data table is not acceptance-corrected yet, so its φ modulations carry
  the detector acceptance. Most of the kinematics converges in two iterations;
  a tail of a few percent of the grid does not.
* The exclusive tail stays out of the weights until the table is absolutely
  normalised (it is the one part that does not cancel in the ratio).
* Not included, as in HAPRAD 2.0: radiation from the hadrons and two-photon
  exchange.

## Authors

C++ port: H. Hakobyan, R. Oyarzun and S. Mancilla. Adaptation for this
repository: A. Bórquez (2021). RG-E tools, validation and farm workflow:
S. Fonseca (2026).

## References

* I. Akushevich, N. Shumeiko, A. Soroko, Eur. Phys. J. C **10**, 681 (1999)
* I. Akushevich, A. Ilyichev, M. Osipenko, Phys. Lett. B **672**, 35 (2009)
* D. Drechsel, S. S. Kamalov, L. Tiator, Eur. Phys. J. A **34**, 69 (2007) (MAID2007)
