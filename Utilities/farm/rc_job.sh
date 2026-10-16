#!/bin/bash
# rc_job.sh -- one job of the RC workflow built by submit_rc_workflow.sh, timed.
#
#   rc_job.sh <workflow dir> <step> <iteration> [run | chunk]
#
# Steps, for iteration k (output under <workflow dir>/iter_k):
#   fill  <run>    MakePhiTable fill over the run's ntuples -> counts/<run>.root;
#                  from k = 1 on, weighted by the run's RC friends of k-1
#   fit            hadd the per-run counts, MakePhiTable fit -> table.root
#   grid  <chunk>  MakeRCGrid run on one chunk of the nodes -> grid_chunks/
#   merge          hadd the chunks -> grid.root; from k = 1 on, MakeRCGrid
#                  compare against the grid of k-1 -> compare.txt
#   apply <run>    ApplyRC on each ntuple of the run -> friends/<run>/
#
# Each job appends one line to timing/<job>.tsv: name, step, iteration,
# argument, host, start and end (epoch s), wall, user and system CPU (s),
# maximum RSS (kB), exit status, and two work counts (fill: events read and
# selected; fit: cells written; grid: nodes and valid nodes; apply: entries and
# selected). timing_report.py summarises them. The job's output goes to
# logs/<job>.log.

set -uo pipefail
if [ $# -lt 3 ]; then
  echo "usage: $0 <workflow dir> <step> <iteration> [run | chunk]" >&2
  exit 2
fi
OUT=$1 STEP=$2 ITER=$3 ARG=${4:-}

# ==== Inner call: the step itself, run under /usr/bin/time by the outer call ====
if [ "${RC_JOB_INNER:-}" = 1 ]; then
  source "$OUT/farm.cfg"
  BIN=$HAPRAD_CPP/Utilities/bin
  CFG=$OUT/config.cfg
  I=$OUT/iter_$ITER
  P=$OUT/iter_$((ITER - 1))

  # The ntuples of run $1, from DATA_PATTERN.
  run_files() {
    shopt -s nullglob
    local pat=${DATA_PATTERN//%RUN%/$1}
    pat=${pat//%RUNNO%/$((10#$1))}
    pat=${pat//%TARGET%/$TARGET}
    # A pattern without wildcards expands to itself: keep existing files only.
    FILES=()
    local f
    for f in $pat; do [ -e "$f" ] && FILES+=("$f"); done
    shopt -u nullglob
    if [ ${#FILES[@]} -eq 0 ]; then
      echo "no input files match $pat" >&2
      exit 1
    fi
  }

  set -e
  case $STEP in
    fill)
      run_files "$ARG"
      mkdir -p "$I/counts"
      friend=()
      [ "$ITER" -gt 0 ] && friend=(--rc-friend "$P/friends/$ARG")
      "$BIN/MakePhiTable" fill "$CFG" "$I/counts/$ARG.root" "${friend[@]}" "${FILES[@]}"
      ;;
    fit)
      hadd -f "$I/counts.root" "$I"/counts/*.root > /dev/null
      "$BIN/MakePhiTable" fit "$CFG" "$I/counts.root" "$I/table.root"
      ;;
    grid)
      mkdir -p "$I/grid_chunks"
      "$BIN/MakeRCGrid" run "$CFG" "$I/table.root" "$ARG" "$NCHUNKS" "$I/grid_chunks/chunk_$ARG.root"
      ;;
    merge)
      hadd -f "$I/grid.root" "$I"/grid_chunks/chunk_*.root > /dev/null
      # The change since the previous grid is a diagnostic: report, never fail.
      if [ "$ITER" -gt 0 ]; then
        "$BIN/MakeRCGrid" compare "$P/grid.root" "$I/grid.root" > "$I/compare.txt" 2>&1 ||
          echo "WARNING: grid comparison with iteration $((ITER - 1)) failed" >> "$I/compare.txt"
        cat "$I/compare.txt"
      fi
      ;;
    apply)
      run_files "$ARG"
      mkdir -p "$I/friends/$ARG"
      for f in "${FILES[@]}"; do
        b=$(basename "$f" .root)
        "$BIN/ApplyRC" "$CFG" "$I/grid.root" "$f" "$I/friends/$ARG/${b}_rc.root"
      done
      ;;
    *)
      echo "unknown step $STEP" >&2
      exit 2
      ;;
  esac
  exit 0
fi

# ==== Outer call: environment, scratch directory, timing ====
source "$OUT/farm.cfg"
eval "$ENV_SETUP"
NAME=${STEP}_i${ITER}${ARG:+_$ARG}
mkdir -p "$OUT/logs" "$OUT/timing"
LOG=$OUT/logs/$NAME.log
TIMES=$OUT/timing/$NAME.time

# HAPRAD reads pi_n_maid.dat from the working directory: run in a scratch
# directory that links it.
WORK=$(mktemp -d "${TMPDIR:-/tmp}/rcjob.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
ln -s "$HAPRAD_DATA_DIR/pi_n_maid.dat" "$WORK/pi_n_maid.dat"
cd "$WORK"

{
  echo "# $NAME on $(uname -n) at $(date '+%F %T')"
  command -v root-config > /dev/null || echo "WARNING: root-config is not on PATH (set ENV_SETUP in farm.cfg)"
} > "$LOG"

start=$(date +%s.%N)
if [ -x /usr/bin/time ]; then
  RC_JOB_INNER=1 /usr/bin/time -f "%e %U %S %M" -o "$TIMES" "$0" "$@" >> "$LOG" 2>&1
  status=$?
  read -r wall user sys rss < <(tail -n 1 "$TIMES") || true
else
  RC_JOB_INNER=1 "$0" "$@" >> "$LOG" 2>&1
  status=$?
  wall="" user=NA sys=NA rss=NA
fi
end=$(date +%s.%N)
wall=$(awk -v a="$start" -v b="$end" 'BEGIN { printf "%.1f", b - a }')
rm -f "$TIMES"

# ---- Work counts from the step's own output ----
w1=NA w2=NA
case $STEP in
  fill)
    w1=$(awk '$1 == "read" { print $2 }' "$LOG" | tail -n 1)
    w2=$(awk '$1 == "inside" && $2 == "grid" { print $3 }' "$LOG" | tail -n 1)
    ;;
  fit) w1=$(awk '$1 == "written" { print $2 }' "$LOG" | tail -n 1) ;;
  grid)
    w1=$(sed -n 's/.* \([0-9]*\) nodes -- .*/\1/p' "$LOG" | tail -n 1)
    w2=$(sed -n 's/.* rejected by HAPRAD, \([0-9]*\) valid.*/\1/p' "$LOG" | tail -n 1)
    ;;
  apply)
    w1=$(sed -n 's/.*: \([0-9]*\) entries, \([0-9]*\) selected$/\1/p' "$LOG" | awk '{ s += $1 } END { print s + 0 }')
    w2=$(sed -n 's/.*: \([0-9]*\) entries, \([0-9]*\) selected$/\2/p' "$LOG" | awk '{ s += $1 } END { print s + 0 }')
    ;;
esac

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$NAME" "$STEP" "$ITER" "${ARG:--}" "$(uname -n)" \
  "$start" "$end" "$wall" "${user:-NA}" "${sys:-NA}" "${rss:-NA}" "$status" "${w1:-NA}" "${w2:-NA}" > "$OUT/timing/$NAME.tsv"
echo "# exit $status after $wall s" >> "$LOG"
exit $status
