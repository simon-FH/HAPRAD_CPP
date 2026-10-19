#!/bin/bash
# export_results.sh -- pack a finished workflow's results for plotting on
# another machine: event histograms of every iteration (event_histograms.C),
# and each iteration's RC grid and table, into one tarball.
#
#   ./export_results.sh <workflow dir> [out.tgz]
#
# It takes a while, so run it detached from the terminal (tcsh suspends
# background jobs that touch it). In tcsh, >&! overwrites an existing log even
# with noclobber set; plain >& would fail, and the error message alone
# suspends the job:
#
#   setsid nohup ./export_results.sh <workflow dir> < /dev/null >&! export.log &
#
# Writes <workflow dir>/export/ (event_histograms.root, grid_<k>.root,
# table_<k>.root, the workflow's settings and run list) and the tarball,
# by default <workflow dir>/export.tgz. Needs ROOT on PATH. The histograms
# read every selected event of every run, so this takes a while (tens of
# minutes for a full target).

set -euo pipefail
if [ $# -lt 1 ]; then
  sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//' >&2
  exit 2
fi
W=$(readlink -f "$1")
OUT=$(readlink -f "${2:-$W/export.tgz}")
HERE=$(dirname "$(readlink -f "$0")")
E=$W/export
mkdir -p "$E"

# ---- Grids and tables of every iteration, under distinct names ----
for d in "$W"/iter_*; do
  k=${d##*_}
  cp "$d/grid.root" "$E/grid_$k.root"
  cp "$d/table.root" "$E/table_$k.root"
done
cp "$W/farm.cfg" "$W/runs.txt" "$W/jobs.tsv" "$E/"
[ -d "$W/timing" ] && tar -C "$W" -cf "$E/timing.tar" timing

# ---- Event histograms (the long part) ----
echo "event histograms: started $(date '+%F %T')"
# stdin from /dev/null: ROOT otherwise sets up the terminal for input, which
# suspends the script when it runs in the background ("tty output").
root -l -b -q "$HERE/event_histograms.C(\"$W\", \"$E/event_histograms.root\")" < /dev/null | grep -v "^  \|^Processing"
echo "event histograms: done $(date '+%F %T')"

# ---- One file to copy ----
tar -C "$W" -czf "$OUT" export
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
