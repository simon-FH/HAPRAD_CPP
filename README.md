HAPRAD_CPP
==========

C++ program for calculation of radiative corrections to semi-inclusive hadron leptoproduction, based in the original Fortran code `HAPRAD2`.

## Requirements

* [**ROOT**](https://root.cern.ch/) (built with `MathMore`, i.e. GSL)

* A Fortran compiler (`gfortran`) -- used for `pkhff.f`.

**CERNLIB is not required.** It used to be, because `init_pdf.f` wraps the PDFLIB
entry points `PDFSET` and `STRUCTM`; but nothing in the C++ code calls it, and
building it in only dragged two unresolvable symbols into `libTRadCor.so`, which
then broke the link of every program using the library. It is now left out of
the build by default.

CERNLIB is only needed if the PDF x FF semi-inclusive structure-function model
is restored in `TSemiInclusiveModel.cxx` (it currently reads pre-fitted
`A`/`Ac`/`Acc` parameters instead). See "Optional: CERNLIB/PDFLIB" below.

## Compilation

Compile by running `make`.

### Optional: CERNLIB / PDFLIB

Modern Linux distributions no longer package CERNLIB (Debian dropped it, and
there is no Fedora package), but PDFLIB 8.04 still ships inside the
community-maintained CERNLIB sources and builds standalone. Only the `mclibs/pdf`
subtree is needed -- this project has its own `dfint` and integrates via GSL, so
packlib/kernlib/mathlib are not required beyond a handful of support routines.

Once you have a `libpdflib804`:

```bash
export CERN_LIB=/path/to/cernlib/lib
make USE_CERNLIB=1
```

`USE_CERNLIB=1` compiles `init_pdf.f` and links `$(CERNLIB_LIBS)` (default
`-lpdflib804`); override `CERNLIB_LIBS` if your install needs more. The same
flag works for `Utilities/`.

## Usage

If the compilation was succesful, a file called `slib/libTRadCor.so` should now exist.

You can link this dynamic library to your own programs in order to calculate radiative corrections. Please, refer to [Utilities/README.md](Utilities/README.md) to know how to use the programs present in this repository.

## Authors

H. Hakobyan, R. Oyarzun and S. Mancilla.

## References

* I. Akushevich, N. Shumeiko, A. Soroko. Eur. Phys. J. C **10**, 681 (1999)

* I. Akushevich, A. Ilyichev, M. Osipenko. Physics Letters B **672**, 35-44 (2009)
