// weight_summary.C -- the RC weights of one workflow iteration, over all its
// friend files: how the selected events divide by rc_status, the quantiles of
// their RC factor, and the fractions far from 1 (< 0.5, > 1.5, > 2, > 10) --
// the share of the data that extreme grid nodes affect.
//
//   root -l -b -q 'weight_summary.C("<workflow dir>/iter_<k>")'

void weight_summary(const char* iterDir) {
  // ==== All friend files of the iteration, one directory per run ====
  TChain c("RC");
  const TString list = gSystem->GetFromPipe(Form("ls %s/friends/*/*_rc.root 2>/dev/null", iterDir));
  std::unique_ptr<TObjArray> files(list.Tokenize("\n"));
  for (TObject* f : *files) c.Add(((TObjString*)f)->GetString());
  if (c.GetNtrees() == 0) {
    printf("no friend files under %s/friends\n", iterDir);
    return;
  }

  // ==== Selected events by rc_status ====
  const double nSel = c.GetEntries("selected");
  const char* label[5] = {"full grid support", "partial support", "outside the grid", "no grid support", "other particle/beam"};
  printf("%s: %d files, %.0f selected events\n", iterDir, c.GetNtrees(), nSel);
  for (int s = 0; s < 5; ++s)
    printf("  rc_status %d %-20s %6.2f%%\n", s, label[s], 100. * c.GetEntries(Form("selected && rc_status == %d", s)) / nSel);

  // ==== The RC factor of the selected events that have one ====
  c.SetEstimate(c.GetEntries() + 1);
  const Long64_t n = c.Draw("rc", "selected && rc_status <= 1", "goff");
  if (n <= 0) return;
  std::vector<double> v(c.GetV1(), c.GetV1() + n);
  std::sort(v.begin(), v.end());
  auto q = [&](double x) { return v[std::min<size_t>(v.size() - 1, size_t(x * (v.size() - 1)))]; };
  auto frac = [&](auto pred) { return 100. * std::count_if(v.begin(), v.end(), pred) / v.size(); };
  printf("RC factor of the %lld selected events with one:\n", n);
  printf("  1%% %.3f  5%% %.3f  median %.3f  95%% %.3f  99%% %.3f  max %.4g\n", q(0.01), q(0.05), q(0.5), q(0.95), q(0.99),
         v.back());
  printf("  < 0.5: %.2f%%   > 1.5: %.2f%%   > 2: %.2f%%   > 10: %.2f%%\n", frac([](double x) { return x < 0.5; }),
         frac([](double x) { return x > 1.5; }), frac([](double x) { return x > 2.; }), frac([](double x) { return x > 10.; }));
}
