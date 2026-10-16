#!/bin/bash
# submit_rc_workflow.sh -- the RC iteration as a swif2 workflow on the JLab farm.
#
#   ./submit_rc_workflow.sh <target> <run list> [options]
#
#   <target> names the workflow and fills %TARGET% in DATA_PATTERN, e.g. C_D2.
#
#     -s <farm.cfg>    settings (default: farm.cfg next to this script)
#     -C <config>      analysis config, overriding CONFIG in the settings
#     -i <n>           RC iterations, overriding ITERATIONS
#     -c <n>           grid chunks per iteration, overriding NCHUNKS
#     -n <name>        workflow name (default: rc_<field>_<target>_<tracking>_<date>_<time>)
#     --tracking <t>   dc or dc_fmt, overriding TRACKING
#     --field <f>      inb or out, overriding FIELD
#     -d <pattern>     input pattern, overriding DATA_PATTERN
#     --skip-missing   leave out runs without input files, with a warning,
#                      instead of stopping
#     --dry-run        print the swif2 commands, submit nothing
#     --local [jobs]   run the same jobs on this machine instead, stage by
#                      stage, up to <jobs> at a time (default: all cores)
#
# The run list has one run number per line; blank lines and # comments are
# ignored. Per iteration k the jobs are, with their dependencies:
#
#   fill_ik_<run>   one per run      after apply_i(k-1)_<run> (k > 0)
#   fit_ik          one              after every fill_ik_*
#   grid_ik_<c>     NCHUNKS          after fit_ik
#   merge_ik        one              after every grid_ik_*
#   apply_ik_<run>  one per run      after merge_ik
#
# so a run's next fill starts as soon as its own weights exist. Everything is
# written under OUT_BASE/<workflow name>: the settings, config and run list as
# used, jobs.tsv (every job and its dependencies), iter_k/ (counts, table,
# grid, friends), logs/ and timing/. When the workflow is done,
#
#   ./timing_report.py OUT_BASE/<workflow name> [--runs N]
#
# gives the benchmark. Watch progress with `swif2 status <workflow name>`.

set -euo pipefail
HERE=$(dirname "$(readlink -f "$0")")

# ==== Arguments ====
if [ $# -lt 2 ]; then
  sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//' >&2
  exit 2
fi
TARGET=$1
RUNLIST=$(readlink -f "$2")
shift 2
SETTINGS=$HERE/farm.cfg
OPT_CONFIG="" OPT_ITER="" OPT_CHUNKS="" OPT_PATTERN="" OPT_TRACKING="" OPT_FIELD="" NAME="" MODE=swif2 LOCAL_JOBS=$(nproc) SKIP_MISSING=0
while [ $# -gt 0 ]; do
  case $1 in
    -s) SETTINGS=$(readlink -f "$2"); shift 2 ;;
    -C) OPT_CONFIG=$(readlink -f "$2"); shift 2 ;;
    -i) OPT_ITER=$2; shift 2 ;;
    -c) OPT_CHUNKS=$2; shift 2 ;;
    -n) NAME=$2; shift 2 ;;
    -d) OPT_PATTERN=$2; shift 2 ;;
    --tracking) OPT_TRACKING=$2; shift 2 ;;
    --field) OPT_FIELD=$2; shift 2 ;;
    --skip-missing) SKIP_MISSING=1; shift ;;
    --dry-run) MODE=dry; shift ;;
    --local)
      MODE=local
      shift
      if [ $# -gt 0 ] && [[ $1 =~ ^[0-9]+$ ]]; then LOCAL_JOBS=$1; shift; fi
      ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done

source "$SETTINGS"
[ -n "$OPT_CONFIG" ] && CONFIG=$OPT_CONFIG
[ -n "$OPT_ITER" ] && ITERATIONS=$OPT_ITER
[ -n "$OPT_CHUNKS" ] && NCHUNKS=$OPT_CHUNKS
[ -n "$OPT_PATTERN" ] && DATA_PATTERN=$OPT_PATTERN
[ -n "$OPT_TRACKING" ] && TRACKING=$OPT_TRACKING
[ -n "$OPT_FIELD" ] && FIELD=$OPT_FIELD
NAME=${NAME:-rc_${FIELD}_${TARGET}_${TRACKING}_$(date +%Y%m%d_%H%M)}
OUT=$OUT_BASE/$NAME
RC_JOB=$HERE/rc_job.sh

# ==== Checks before anything is created ====
fail() { echo "ERROR: $*" >&2; exit 1; }
[ -r "$RUNLIST" ] || fail "cannot read run list $RUNLIST"
[ -r "$CONFIG" ] || fail "cannot read config $CONFIG"
[ -r "$HAPRAD_DATA_DIR/pi_n_maid.dat" ] || fail "no pi_n_maid.dat in $HAPRAD_DATA_DIR"
for b in MakePhiTable MakeRCGrid ApplyRC; do
  [ -x "$HAPRAD_CPP/Utilities/bin/$b" ] || fail "$HAPRAD_CPP/Utilities/bin/$b missing -- build first (make, then make in Utilities)"
done
[ "$MODE" = swif2 ] && ! command -v swif2 > /dev/null && fail "swif2 not found (run on ifarm, or use --dry-run / --local)"
[ -e "$OUT" ] && fail "$OUT already exists -- choose another name with -n"

# ---- Runs: six-digit numbers, each with input files ----
LISTED=()
while read -r line; do
  line=${line%%#*}
  line=${line//[[:space:]]/}
  [ -z "$line" ] && continue
  [[ $line =~ ^[0-9]+$ ]] || fail "not a run number in $RUNLIST: '$line'"
  LISTED+=("$(printf "%06d" $((10#$line)))")
done < "$RUNLIST"
[ ${#LISTED[@]} -gt 0 ] || fail "no runs in $RUNLIST"
RUNS=() MISSING=()
shopt -s nullglob
for r in "${LISTED[@]}"; do
  pat=${DATA_PATTERN//%RUN%/$r}
  pat=${pat//%RUNNO%/$((10#$r))}
  pat=${pat//%TARGET%/$TARGET}
  pat=${pat//%FIELD%/$FIELD}
  pat=${pat//%TRACKING%/$TRACKING}
  # A pattern without wildcards expands to itself: keep existing files only.
  files=()
  for f in $pat; do [ -e "$f" ] && files+=("$f"); done
  if [ ${#files[@]} -eq 0 ]; then
    MISSING+=("$r")
    echo "run $r: NO FILES match $pat"
  else
    RUNS+=("$r")
    echo "run $r: ${#files[@]} file(s)"
  fi
done
shopt -u nullglob
if [ ${#MISSING[@]} -gt 0 ]; then
  [ "$SKIP_MISSING" = 1 ] || fail "${#MISSING[@]} run(s) without input files: ${MISSING[*]} (--skip-missing leaves them out)"
  echo "WARNING: leaving out ${#MISSING[@]} run(s) without input files: ${MISSING[*]}"
fi
[ ${#RUNS[@]} -gt 0 ] || fail "no run has input files"

# ==== Workflow directory: everything the jobs read, frozen ====
# A dry run writes these to a temporary directory instead, removed on exit.
META=$OUT
if [ "$MODE" = dry ]; then
  META=$(mktemp -d)
  trap 'rm -rf "$META"' EXIT
fi
mkdir -p "$META/logs" "$META/timing"
cp "$SETTINGS" "$META/farm.cfg"
printf "\n# Set by submit_rc_workflow.sh\nTARGET=%s\nFIELD=%s\nTRACKING=%s\nDATA_PATTERN='%s'\nCONFIG=%s\nITERATIONS=%s\nNCHUNKS=%s\n" \
  "$TARGET" "$FIELD" "$TRACKING" "$DATA_PATTERN" "$OUT/config.cfg" "$ITERATIONS" "$NCHUNKS" >> "$META/farm.cfg"
cp "$CONFIG" "$META/config.cfg"
printf "%s\n" "${RUNS[@]}" > "$META/runs.txt"
printf "name\tstep\titeration\targ\tafter\n" > "$META/jobs.tsv"
date +%s.%N > "$META/timing/submitted"

# ==== Adding jobs ====
# add_job <name> <step> <iter> <arg> <"time ram disk"> [antecedents ...]
STAGE=()
add_job() {
  local name=$1 step=$2 iter=$3 arg=$4 res=$5
  shift 5
  local t ram disk
  read -r t ram disk <<< "$res"
  printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$step" "$iter" "${arg:--}" "$(IFS=,; echo "${*:-}")" >> "$META/jobs.tsv"
  local cmd="$RC_JOB $OUT $step $iter${arg:+ $arg}"
  case $MODE in
    local) STAGE+=("$cmd") ;;
    *)
      local a=(swif2 add-job "$NAME" -name "$name" -partition "$PARTITION" -cores 1 -time "$t" -ram "$ram" -disk "$disk" -shell /bin/bash)
      [ -n "$ACCOUNT" ] && a+=(-account "$ACCOUNT")
      local x
      for x in "$@"; do a+=(-antecedent "$x"); done
      if [ "$MODE" = dry ]; then
        echo "${a[*]} \"$cmd\""
      else
        "${a[@]}" "$cmd" > /dev/null
      fi
      ;;
  esac
}

# Local mode: run the jobs added since the last call, in parallel, and stop at
# the first failure (swif2 would hold back everything that depends on it).
end_stage() {
  [ "$MODE" = local ] || return 0
  [ ${#STAGE[@]} -gt 0 ] || return 0
  echo "local: ${#STAGE[@]} job(s): ${STAGE[0]##*rc_job.sh $OUT }${STAGE[1]:+ ...}"
  printf "%s\n" "${STAGE[@]}" | xargs -P "$LOCAL_JOBS" -I{} sh -c '{}' ||
    fail "a job failed; see $OUT/logs"
  STAGE=()
}

case $MODE in
  swif2) swif2 create "$NAME" ;;
  dry) echo "swif2 create $NAME" ;;
esac

# ==== The iterations ====
for ((k = 0; k < ITERATIONS; k++)); do
  fills=()
  for r in "${RUNS[@]}"; do
    if [ "$k" -gt 0 ]; then
      add_job "fill_i${k}_$r" fill "$k" "$r" "$RES_FILL" "apply_i$((k - 1))_$r"
    else
      add_job "fill_i${k}_$r" fill "$k" "$r" "$RES_FILL"
    fi
    fills+=("fill_i${k}_$r")
  done
  end_stage
  add_job "fit_i$k" fit "$k" "" "$RES_FIT" "${fills[@]}"
  end_stage
  chunks=()
  for ((c = 0; c < NCHUNKS; c++)); do
    add_job "grid_i${k}_$c" grid "$k" "$c" "$RES_GRID" "fit_i$k"
    chunks+=("grid_i${k}_$c")
  done
  end_stage
  add_job "merge_i$k" merge "$k" "" "$RES_MERGE" "${chunks[@]}"
  end_stage
  for r in "${RUNS[@]}"; do
    add_job "apply_i${k}_$r" apply "$k" "$r" "$RES_APPLY" "merge_i$k"
  done
  end_stage
done

# ==== Start ====
njobs=$(($(wc -l < "$META/jobs.tsv") - 1))
case $MODE in
  swif2)
    swif2 run "$NAME"
    echo "workflow $NAME: $njobs jobs submitted; output in $OUT"
    echo "  swif2 status $NAME        progress"
    echo "  $HERE/timing_report.py $OUT   benchmark, when done"
    ;;
  dry)
    echo "swif2 run $NAME"
    echo "dry run: $njobs jobs; nothing submitted, nothing written"
    ;;
  local)
    echo "local workflow done: $njobs jobs; output in $OUT"
    ;;
esac
