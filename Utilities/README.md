HAPRAD_CPP/Utilities
====================

ROOT/C++ programs that use the dynamic library `TRadCor` to do radiative corrections calculations.

## MakePhiTable -- RG-E data to the structure-function table (CLAS12)

HAPRAD needs a table of the measured azimuthal amplitudes A, Ac, Acc in
(Q2, nu, z, pt2) cells. `MakePhiTable` builds it from RG-E ntuples (the `DT`
trees made by `clas12-rge-analysis`):

```bash
MakePhiTable fill config/rge_pip.cfg counts.root "/path/to/run/*_ntuples_dc.root"
MakePhiTable fit  config/rge_pip.cfg counts.root table.root
```

`fill` applies the cuts and fills 5-D counts; it can run on any subset of the
files, and the outputs of several jobs merge with `hadd`, so it is meant to run
on the cluster where the data lives. `fit` fits each cell and writes the table.
Pass the table to `GetRC -m`, or to `TRadCor::LoadSemiInclusiveTable()`.

Everything -- beam energy, cuts including the missing-mass cut, binning, fit
thresholds -- is in the config file; `config/rge_pip.cfg` documents each key.
One table is built per beam energy.

The table records `acceptance_corrected = no`: acceptance correction is not
applied yet. See `PLAN.md` (Phase 1) for what A is and how it is normalised,
and `validation/closure_producer.py` for the closure test.

## The EG2 chain (CLAS6)

The programs below are the original EG2 analysis chain.

## Environment

Feel free to modify the variables in `set_env.sh` according to your own environment.

## Usage

**Important:** to use these programs your data files must have been filtered by [**GetSimpleTuple**](http://github.com/utfsm-eg2-data-analysis/GetSimpleTuple).

1. Modify `include/Binning.hxx` to the binning of your preference.

2. **Compile** all programs (see below).

3. Execute `./exec_rad-corr_chain.sh --pid <pid>` where `<pid> = (211, -211, 2212)`.

If everything went well, you should had six files with the Radiative Correction Factors called `RCFactor_<target>.txt` and six images corresponding to the fit of the `PhiPQ` distributions for each target.

## Compilation

Do `source set_env.sh` and compile the programs by running `make`.
