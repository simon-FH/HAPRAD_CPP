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
applied yet. An optional `weight_branch` (e.g. 1/acceptance) weights each event.

**The table's absolute scale matters.** HAPRAD adds the exclusive radiative
tail in absolute units (MAID), while the Born cross section and the inelastic
tail scale with the table. Without the optional `luminosity` key (events per
nb) the table is in arbitrary units and the exclusive tail is effectively
switched off; the table records which. See `PLAN.md` (Phase 1) for what A is
and how it is normalised, and `validation/closure_producer.py` and
`validation/closure_normalisation.py` for the closure tests.

## MakeRCGrid, ApplyRC -- radiative-correction weights (CLAS12)

HAPRAD takes ~0.5 s per call, so it is not run per event. It is evaluated once
on a grid of nodes in (Q2, nu, z, p_t, phi), and every event gets the RC factor
interpolated at its own kinematics, written as a friend tree:

```bash
./run_rc_grid.sh config/rge_pip.cfg table.root grid.root [jobs]            # parallel chunks + merge
./check_rc_grid.py config/rge_pip.cfg grid.root table.root 2000 <ntuples>  # interpolation error budget
bin/ApplyRC config/rge_pip.cfg grid.root run_file.root run_file_rc.root   # one job per input file
```

```cpp
dt->AddFriend("RC", "run_file_rc.root");
dt->Draw("Q2", "RC.w * (RC.selected && RC.rc_status <= 1)");
```

* The nodes (`rc_nodes_*`), beam energy and particle come from the same config
  as the table. A grid belongs to one table: each iteration needs a new one.
* On the cluster, `MakeRCGrid run <config> <table> <chunk> <nchunks> <out>` is
  one job of an array; merge the chunks with `hadd`. The grid refuses to load
  if a chunk is missing or merged twice.
* HAPRAD reads its MAID grid `pi_n_maid.dat` from the **working directory**;
  the scripts run in `HAPRAD_DATA_DIR` (default `../haprad2`). Without the file
  the library now stops with an error -- it used to return an exclusive tail of
  exactly zero.
* The grid stores the parts of the RC factor separately (Born x virtual/soft,
  inelastic tail, exclusive tail). `ApplyRC` leaves the exclusive tail out
  unless `rc_exclusive_tail = yes`, which it refuses for a table in arbitrary
  units (see the luminosity note above).
* Friend branches: `rc`, `rc_noex`, `w = 1/rc`, `rc_status` (0 full grid
  support, 1 partial, 2 outside the grid, 3 none, 4 other particle or beam),
  `rc_cover`, `selected` (passes the config's cuts).

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
