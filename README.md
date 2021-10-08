HAPRAD_CPP
==========

C++ program for calculation of radiative corrections to semi-inclusive hadron leptoproduction, based in the original Fortran code `HAPRAD2`.

## Requirements

* [**ROOT**](https://root.cern.ch/)

* **cernlib**, preferably the official package that comes with Debian/Ubuntu, and not the one located in http://github.com/utfsm-eg2-data-analysis/cernlib.

## Compilation

Compile by running `make`.

## Usage

If the compilation was succesful, a file called `slib/libTRadCor.so` should now exist.

You can link this dynamic library to your own programs in order to calculate radiative corrections. Please, refer to [Utilities/README.md](Utilities/README.md) to know how to use the programs present in this repository.

## Authors

H. Hakobyan, R. Oyarzun and S. Mancilla.

## References

* I. Akushevich, N. Shumeiko, A. Soroko. Eur. Phys. J. C **10**, 681 (1999)

* I. Akushevich, A. Ilyichev, M. Osipenko. Physics Letters B **672**, 35-44 (2009)
