#!/bin/bash
# run_rc_grid.sh -- build an RC grid in parallel chunks and merge them
# (PLAN.md Phase 3).
#
#   ./run_rc_grid.sh <config> <table.root> <grid.root> [jobs, default: all cores]
#
# Each chunk is an independent `MakeRCGrid run`, so on the cluster the same
# chunks can go to a job array instead, followed by the same `hadd`.
#
# HAPRAD reads its MAID grid (pi_n_maid.dat) from the working directory, so the
# chunks run in HAPRAD_DATA_DIR (default: haprad2/ in this repository).

set -euo pipefail
if [ $# -lt 3 ]; then
  echo "usage: $0 <config> <table.root> <grid.root> [jobs]" >&2
  exit 2
fi
cfg=$(readlink -f "$1")
table=$(readlink -f "$2")
out=$(readlink -f "$3")
jobs=${4:-$(nproc)}
here=$(dirname "$(readlink -f "$0")")
bin=$here/bin
data=$(readlink -f "${HAPRAD_DATA_DIR:-$here/../haprad2}")
[ -r "$data/pi_n_maid.dat" ] || { echo "no pi_n_maid.dat in $data (set HAPRAD_DATA_DIR)" >&2; exit 1; }
tmp=$(mktemp -d "${out%.root}_chunks.XXXX")
cd "$data"

echo "RC grid: $jobs chunks in parallel, work in $tmp"
seq 0 $((jobs - 1)) | xargs -P "$jobs" -I{} sh -c \
  "'$bin/MakeRCGrid' run '$cfg' '$table' {} $jobs '$tmp/chunk_{}.root' > '$tmp/chunk_{}.log' 2>&1 \
   || { echo 'chunk {} failed, see $tmp/chunk_{}.log' >&2; exit 255; }"
hadd -f "$out" "$tmp"/chunk_*.root > "$tmp/hadd.log" 2>&1 || { echo "hadd failed, see $tmp/hadd.log" >&2; exit 1; }
grep -h " chunk .* of " "$tmp"/chunk_*.log | sed "s|^$tmp/||"
rm -r "$tmp"
echo "wrote $out"
