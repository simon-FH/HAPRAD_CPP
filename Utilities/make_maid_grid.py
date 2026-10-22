#!/usr/bin/env python3
"""Build the exclusive-tail input pi_n_maid.dat from the MAID2007 web calculator.

    ./make_maid_grid.py [--out ../haprad2/pi_n_maid.dat] [--cache maid2007_blocks]

The exclusive radiative tail (TExclusiveModel, and exclusive_model.f in the
FORTRAN) reads the gamma* p -> pi+ n 2-fold cross sections from
pi_n_maid.dat: one line per (Q2, W, theta) with dsT, dsL, dsTT, dsTL, dsTL' in
microbarn/sr (pi N centre-of-mass frame, MAID convention: the longitudinal
parts without epsilon_L), on the grid

    Q2     0, 0.3, ..., 4.8 GeV^2 and 5.0      18 values
    W      1.08, 1.10, ..., 2.00 GeV           47 values
    theta  0, 3, ..., 180 deg                  61 values

looping theta fastest, then W, then Q2 (51,606 lines, no header: both codes read
the numbers in that order). This script asks the MAID2007 "observables" page
(maid.kph.uni-mainz.de, "Differential Cross Sections", channel pi+ n, default
model parameters) for one theta scan per (Q2, W), checks each answer, keeps it
in --cache so an interrupted run resumes, and writes the file. The server
prints four decimals, as in the MAID2003 file it replaces.
"""

import argparse
import os
import re
import sys
import time
import urllib.parse
import urllib.request

URL = "http://maid.kph.uni-mainz.de/cgi-bin/maid1"
Q2_GRID = [round(0.3 * i, 1) for i in range(17)] + [5.0]
W_GRID = [round(1.08 + 0.02 * i, 2) for i in range(47)]
N_THETA = 61

# The MAID2007 observables form (maid2007/observ.html) with its defaults: every
# model ingredient on (param11..26), standard couplings (value11..13, 51..84),
# and the 2-fold cross sections (switch 211).
DEFAULTS = ([("switch", "211"), ("param99", "0")] + [("param%d" % i, "1") for i in range(11, 27)] +
            [("value%d" % i, "1.0") for i in list(range(11, 14)) + list(range(51, 85))])


def theta_scan(q2, w):
    """The 61 rows (theta, dsT, dsL, dsTT, dsTL, dsTL') at Q2 [GeV^2], W [GeV]."""
    query = DEFAULTS + [("param2", "3"), ("param50", "3"), ("value35", "%g" % q2), ("value36", "%g" % round(w * 1000, 3)),
                        ("value37", "0"), ("value41", "3"), ("value42", "180")]
    url = URL + "?" + urllib.parse.urlencode(query)
    page = urllib.request.urlopen(url, timeout=120).read().decode("latin-1")
    text = re.sub(r"<[^>]+>", "", page)
    if "pi+ + neutron" not in text:
        raise RuntimeError("not a pi+ n answer")
    rows = []
    for line in text.splitlines():
        f = line.split()
        if len(f) < 10:
            continue
        try:
            rows.append([float(x) for x in f[:6]])
        except ValueError:
            continue
    if len(rows) != N_THETA or any(abs(r[0] - 3 * i) > 1e-6 for i, r in enumerate(rows)):
        raise RuntimeError("expected theta = 0, 3, ..., 180, got %d rows" % len(rows))
    return rows, url


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(here, "..", "haprad2", "pi_n_maid.dat"))
    ap.add_argument("--cache", default="maid2007_blocks", help="directory for the per-(Q2, W) answers")
    ap.add_argument("--pause", type=float, default=0.5, help="seconds between requests")
    args = ap.parse_args()
    os.makedirs(args.cache, exist_ok=True)

    # ---- One theta scan per (Q2, W), cached ----
    blocks = {}
    for iq, q2 in enumerate(Q2_GRID):
        for iw, w in enumerate(W_GRID):
            path = os.path.join(args.cache, "q%02d_w%02d.txt" % (iq, iw))
            if not os.path.exists(path):
                for attempt in range(5):
                    try:
                        rows, url = theta_scan(q2, w)
                        break
                    except Exception as e:  # network or server hiccup: wait and retry
                        print("Q2 %g W %g: %s, retrying" % (q2, w, e), file=sys.stderr)
                        time.sleep(10)
                else:
                    sys.exit("giving up at Q2 = %g, W = %g" % (q2, w))
                with open(path + ".tmp", "w") as f:
                    f.write("# MAID2007 pi+ n, Q2=%g GeV2, W=%g GeV; %s\n" % (q2, w, url))
                    for r in rows:
                        f.write(" ".join("%.4f" % x for x in r) + "\n")
                os.rename(path + ".tmp", path)
                time.sleep(args.pause)
            blocks[iq, iw] = [[float(x) for x in l.split()] for l in open(path) if not l.startswith("#")]
        print("Q2 = %g: %d W values" % (q2, len(W_GRID)), flush=True)

    # ---- The file, theta fastest, then W, then Q2 ----
    with open(args.out, "w") as out:
        for iq in range(len(Q2_GRID)):
            for iw in range(len(W_GRID)):
                for r in blocks[iq, iw]:
                    out.write("".join("%9.4f" % x for x in r[1:6]) + "\n")
    print("wrote %s: %d lines" % (args.out, len(Q2_GRID) * len(W_GRID) * N_THETA))


if __name__ == "__main__":
    main()
