// convergence_report.C -- how far the RC iteration has converged, checked on
// the grid, the table and the events, with flags for divergence.
//
//   root -l -b -q 'convergence_report.C("<dir>")'
//   root -l -b -q 'convergence_report.C("<dir>", 0.002, 0.5, 0.002)'
//
// <dir> is a workflow directory (iter_<k>/grid.root, iter_<k>/table.root) or
// the export directory export_results.sh writes (grid_<k>.root, table_<k>.root,
// event_histograms.root). Every pair of consecutive iterations k-1 -> k is
// checked:
//
//   grid    change of eta = r_vr + r_in at the nodes valid in both;
//           converged when the 95th percentile of |change| < tolGrid
//   table   change of A, Ac, Acc in the cells fitted in both, divided by the
//           fit error of iteration k; converged when the 95th percentile of
//           |change| / error < tolPull for all three
//   events  relative change of each selected event's weight, |w_k/w_(k-1) - 1|
//           (event_histograms.C, export directory only); converged when its
//           95th percentile < tolWeight
//
// Divergence flags, from k >= 2: the largest grid change grows; the 99th
// percentile of the grid change shrinks by less than half; more than 2% of
// the nodes change more than in the step before; the share of events with
// rc > 2 grows (export directory only).

// ---- Files of iteration k in either layout ----
TString IterFile(const TString& dir, const char* what, int k) {
  TString flat = Form("%s/%s_%d.root", dir.Data(), what, k);
  if (!gSystem->AccessPathName(flat)) return flat;
  TString nested = Form("%s/iter_%d/%s.root", dir.Data(), k, what);
  if (!gSystem->AccessPathName(nested)) return nested;
  return "";
}

// The q-quantile of v, by rank (v is sorted in place).
double Quantile(std::vector<double>& v, double q) {
  if (v.empty()) return NAN;
  std::sort(v.begin(), v.end());
  return v[std::min(v.size() - 1, size_t(q * (v.size() - 1) + 0.5))];
}

// The q-quantile of a histogram of log10(x), back in x; the underflow counts as
// x = 0.
double HistQuantile(const TH1D* h, double q) {
  const double total = h->Integral(0, h->GetNbinsX() + 1);
  if (total <= 0.) return NAN;
  double sum = h->GetBinContent(0);
  if (sum >= q * total) return 0.;
  for (int b = 1; b <= h->GetNbinsX() + 1; ++b) {
    sum += h->GetBinContent(b);
    if (sum >= q * total) return std::pow(10., h->GetXaxis()->GetBinUpEdge(std::min(b, h->GetNbinsX())));
  }
  return std::pow(10., h->GetXaxis()->GetXmax());
}

const char* Verdict(bool ok) { return ok ? "CONVERGED" : "NOT converged"; }

void convergence_report(const char* dirIn, double tolGrid = 0.002, double tolPull = 0.5, double tolWeight = 0.002) {
  const TString dir(dirIn);
  int nIter = 0;
  while (IterFile(dir, "grid", nIter) != "") ++nIter;
  printf("Convergence of the RC iteration in %s: %d iteration(s) with a grid\n", dirIn, nIter);
  printf("Tolerances: grid 95%% |d eta| < %g; table 95%% |change|/error < %g; events 95%% |dw/w| < %g\n\n", tolGrid, tolPull,
         tolWeight);
  if (nIter < 2) {
    printf("Nothing to compare: at least two iterations are needed.\n");
    return;
  }

  // ==== Grids: eta = r_vr + r_in at every node, valid flag ====
  std::vector<std::vector<double>> eta(nIter);
  std::vector<std::vector<char>> valid(nIter);
  for (int k = 0; k < nIter; ++k) {
    TFile f(IterFile(dir, "grid", k));
    THnD *v = f.Get<THnD>("valid"), *vr = f.Get<THnD>("r_vr"), *in = f.Get<THnD>("r_in");
    const Long64_t n = v->GetNbins();
    eta[k].resize(n);
    valid[k].resize(n);
    for (Long64_t b = 0; b < n; ++b) {
      valid[k][b] = v->GetBinContent(b) > 0.5;
      eta[k][b] = vr->GetBinContent(b) + in->GetBinContent(b);
    }
  }

  // ==== Event histograms, if this is an export directory ====
  std::unique_ptr<TFile> ev;
  if (!gSystem->AccessPathName(dir + "/event_histograms.root")) ev.reset(TFile::Open(dir + "/event_histograms.root"));
  auto shareAbove2 = [&](int k) -> double {
    TH1D* h = ev ? ev->Get<TH1D>(Form("h_logrc_%d", k)) : nullptr;
    TH1D* le0 = ev ? ev->Get<TH1D>(Form("n_le0_%d", k)) : nullptr;
    if (!h) return NAN;
    const double total = h->Integral(0, h->GetNbinsX() + 1) + (le0 ? le0->GetBinContent(1) : 0.);
    return h->Integral(h->FindFixBin(std::log10(2.)), h->GetNbinsX() + 1) / total;
  };

  std::vector<double> prevChange;  // |d eta| of step k-1 -> k, by node, for the growth flag
  bool allConverged = true, anyDivergence = false;
  for (int k = 1; k < nIter; ++k) {
    printf("==== Iteration %d -> %d ====\n", k - 1, k);

    // ---- Grid ----
    std::vector<double> d, byNode(eta[k].size(), -1.);
    Long64_t oneOnly = 0;
    for (size_t b = 0; b < eta[k].size(); ++b) {
      if (valid[k][b] != valid[k - 1][b]) ++oneOnly;
      if (!valid[k][b] || !valid[k - 1][b]) continue;
      byNode[b] = std::fabs(eta[k][b] - eta[k - 1][b]);
      d.push_back(byNode[b]);
    }
    const size_t nNodes = d.size();
    Long64_t aboveTol = 0;
    for (double x : d) aboveTol += x > tolGrid;
    const double q50 = Quantile(d, 0.5), q95 = Quantile(d, 0.95), q99 = Quantile(d, 0.99), qMax = d.empty() ? NAN : d.back();
    const bool gridOk = q95 < tolGrid;
    printf("grid    %zu nodes (+%lld valid in only one): |d eta| median %.3g, 95%% %.3g, 99%% %.3g, max %.3g; "
           "%.1f%% above %g  -> %s\n",
           nNodes, oneOnly, q50, q95, q99, qMax, 100. * aboveTol / std::max<size_t>(nNodes, 1), tolGrid, Verdict(gridOk));

    // ---- Table ----
    TFile ta(IterFile(dir, "table", k - 1)), tb(IterFile(dir, "table", k));
    const char* amp[3] = {"A", "Ac", "Acc"};
    bool tableOk = true;
    for (int i = 0; i < 3; ++i) {
      THnD *a = ta.Get<THnD>(amp[i]), *b = tb.Get<THnD>(amp[i]), *e = tb.Get<THnD>(Form("%s_err", amp[i]));
      THnD *fa = ta.Get<THnD>("fitted"), *fb = tb.Get<THnD>("fitted");
      std::vector<double> pull;
      for (Long64_t c = 0; c < a->GetNbins(); ++c) {
        if (fa->GetBinContent(c) == 0. || fb->GetBinContent(c) == 0.) continue;
        const double err = e->GetBinContent(c);
        if (!(err > 0.)) continue;
        pull.push_back(std::fabs(b->GetBinContent(c) - a->GetBinContent(c)) / err);
      }
      Long64_t above1 = 0;
      for (double x : pull) above1 += x > 1.;
      const size_t nCells = pull.size();
      const double p50 = Quantile(pull, 0.5), p95 = Quantile(pull, 0.95);
      tableOk = tableOk && p95 < tolPull;
      printf("table   %-3s %zu cells fitted in both: |change|/error median %.3g, 95%% %.3g; %.1f%% above 1\n", amp[i], nCells,
             p50, p95, 100. * above1 / std::max<size_t>(nCells, 1));
    }
    printf("table   -> %s\n", Verdict(tableOk));

    // ---- Events ----
    bool eventsOk = true, haveEvents = false;
    TH1D* dw = ev ? ev->Get<TH1D>(Form("h_dw_%d", k)) : nullptr;
    if (dw) {
      haveEvents = true;
      TH1D* lost = ev->Get<TH1D>(Form("n_dw_%d", k));
      const double w50 = HistQuantile(dw, 0.5), w95 = HistQuantile(dw, 0.95), w99 = HistQuantile(dw, 0.99);
      const double total = dw->Integral(0, dw->GetNbinsX() + 1);
      const double above = dw->Integral(dw->FindFixBin(std::log10(tolWeight)), dw->GetNbinsX() + 1);
      eventsOk = w95 < tolWeight;
      printf("events  %.0f selected events with a factor in both (+%.0f in only one): |dw/w| median %.2g, 95%% %.2g, "
             "99%% %.2g; %.1f%% above %g  -> %s\n",
             total, lost ? lost->GetBinContent(1) : 0., w50, w95, w99, 100. * above / std::max(total, 1.), tolWeight,
             Verdict(eventsOk));
    } else {
      printf("events  no weight-change histogram (needs event_histograms.root from export_results.sh, made with the "
             "current event_histograms.C)\n");
    }

    // ---- Divergence flags ----
    std::vector<std::string> flags;
    if (k >= 2) {
      // Both steps on the same nodes: those valid in all three iterations.
      std::vector<double> prev, cur;
      Long64_t grew = 0, both = 0;
      for (size_t b = 0; b < byNode.size(); ++b) {
        if (byNode[b] < 0. || prevChange[b] < 0.) continue;
        ++both;
        grew += byNode[b] > prevChange[b] && byNode[b] > tolGrid;
        prev.push_back(prevChange[b]);
        cur.push_back(byNode[b]);
      }
      const double p99 = Quantile(prev, 0.99), c99 = Quantile(cur, 0.99);
      const double pMax = prev.empty() ? NAN : prev.back(), cMax = cur.empty() ? NAN : cur.back();
      if (cMax > pMax) flags.push_back(Form("largest change grows: %.3g -> %.3g", pMax, cMax));
      if (c99 > 0.5 * p99) flags.push_back(Form("99th percentile shrinks by less than half: %.3g -> %.3g", p99, c99));
      if (both > 0 && grew > 0.02 * both)
        flags.push_back(Form("%.1f%% of nodes change more than in the step before (and by more than %g)", 100. * grew / both,
                             tolGrid));
    }
    const double s0 = shareAbove2(k - 1), s1 = shareAbove2(k);
    if (std::isfinite(s0) && std::isfinite(s1) && s1 > s0)
      flags.push_back(Form("share of events with rc > 2 grows: %.2f%% -> %.2f%%", 100. * s0, 100. * s1));
    if (flags.empty()) {
      printf("flags   none\n");
    } else {
      for (const auto& f : flags) printf("flag    DIVERGENCE: %s\n", f.c_str());
      anyDivergence = true;
    }
    if (k == nIter - 1) allConverged = gridOk && tableOk && (!haveEvents || eventsOk);
    prevChange = byNode;
    printf("\n");
  }

  // ==== Summary for the last step ====
  printf("==== Last step (%d -> %d): %s%s ====\n", nIter - 2, nIter - 1, allConverged ? "CONVERGED" : "NOT converged",
         anyDivergence ? "; divergence flags raised (more iterations will not help)" : "");
  printf("Not checked here: the change of your observables in the analysis bins, which needs the analysis binning.\n");
}
