HAPRAD_CPP/Utilities
====================

Collection of four ROOT/C++ programs that use the dynamic library `TRadCor` to do radiative corrections calculations.

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
