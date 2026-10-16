#!/usr/bin/env python3
"""Per-target run lists for submit_rc_workflow.sh, from rge_runs.tsv.

    ./make_run_lists.py [--exclude lumi,test_cook,maybe]

Keeps the runs with cooked = YES and writes one list per target group, the
cooking workflow without its trailing block number (inb_C_D2_1 .. _4 ->
inb_C_D2): runs_<group>.txt, one six-digit run per line, with the run type
and the number of events as a comment. --exclude drops runs of those types.
"""

import argparse
import collections
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exclude", default="", help="comma-separated run types to leave out")
    args = ap.parse_args()
    exclude = {t for t in args.exclude.split(",") if t}

    # ---- Read the table: cooked runs, grouped by target ----
    groups = collections.defaultdict(list)
    with open(os.path.join(HERE, "rge_runs.tsv")) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if l.strip() and not l.startswith("#")]
    for run, events, rtype, workflow, cooked in lines[1:]:
        if cooked != "YES" or rtype in exclude:
            continue
        groups[re.sub(r"_\d+$", "", workflow)].append((int(run), int(events), rtype))

    # ---- One list per group ----
    for g, runs in sorted(groups.items()):
        path = os.path.join(HERE, "runs_%s.txt" % g)
        types = collections.Counter(t for _, _, t in runs)
        with open(path, "w") as out:
            out.write("# %s: %d cooked runs, %.2f G events (%s)\n" % (
                g, len(runs), sum(e for _, e, _ in runs) / 1e9, ", ".join("%d %s" % (n, t) for t, n in sorted(types.items()))))
            for run, events, rtype in sorted(runs):
                out.write("%06d  # %-9s %6.1f M\n" % (run, rtype, events / 1e6))
        print("%-22s %4d runs  %6.2f G events" % (os.path.basename(path), len(runs), sum(e for _, e, _ in runs) / 1e9))


if __name__ == "__main__":
    main()
