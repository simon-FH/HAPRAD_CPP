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

## rc_iterate.py -- the RC iteration (CLAS12)

HAPRAD needs Born-level amplitudes; the data give observed ones. The loop
removes the difference:

```bash
./rc_iterate.py config/rge_pip.cfg work_LD2/ "/path/to/run/*_ntuples_dc.root" [--max-iter 4] [--tol 0.002]
```

Iteration 0 builds the table from the data as they are; iteration k fills it
with the RC weights of iteration k-1 (`MakePhiTable fill --rc-friend DIR`,
which multiplies in `RC.w` from `DIR/<name>_rc.root` for each input
`<name>.root`, on top of any `weight_branch`). It stops when the RC factor
changes by less than `--tol` at 95% of the grid nodes; the last iteration's
`friends/` hold the weights for the analysis. Run it per target. Each iteration
costs one RC grid; `--resume` keeps what is already there.
`validation/closure_iteration.py` is its closure test.

## On the JLab farm (swif2)

`farm/` runs the same iteration as a swif2 workflow, one per target: a fill
job per run, one fit, the RC grid in `NCHUNKS` chunk jobs, one merge, and an
`ApplyRC` job per run, with dependencies so that each iteration follows the
last. Every job records its timing; `timing_report.py` turns those records into
a benchmark.

1. Build on ifarm with the ROOT that `farm/farm.cfg` loads in the jobs
   (`module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles`,
   `module load root/6.28.10-gcc11.4.0`; then `make` in `..` and here, with
   `HAPRAD_CPP` set), and check the paths in `farm/farm.cfg`.
2. Check the job list without submitting:
   `farm/submit_rc_workflow.sh C_D2 farm/runs/runs_bench_C_D2.txt -i 2 --dry-run`.
   The first argument fills `%TARGET%` in `DATA_PATTERN`; `FIELD` (inb/out)
   and `TRACKING` (dc/dc_fmt, default dc_fmt) fill the rest and can be changed
   per workflow with `--field` and `--tracking`. A run without input files
   stops the submission; `--skip-missing` leaves such runs out instead.
3. Benchmark on 2-3 runs with the production config and grid:
   `farm/submit_rc_workflow.sh C_D2 farm/runs/runs_bench_C_D2.txt -i 2`, then, when
   `swif2 status <workflow>` shows it done,
   `farm/timing_report.py <OUT_BASE>/<workflow> --runs <number of runs in the full list>`.
   It prints per-step wall time, CPU, memory and queue wait, the expected cost
   of the full list, and resource settings for `farm.cfg`.
4. Then the full run list.

Run lists are in `farm/runs/`: `rge_runs.tsv` is the RG-E run table (run, events,
type, cooking workflow, cooked status) and `make_run_lists.py` writes one list
per target group from it, cooked runs only (`--exclude lumi,test_cook,maybe` to
drop run types). `runs_bench_C_D2.txt` holds the three benchmark runs.

`--local [jobs]` runs the same jobs on the current machine instead, for testing
the workflow without the farm. Outputs land in `<OUT_BASE>/<workflow>/iter_k/`;
the weights for the analysis are the last iteration's `friends/<run>/`.

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
