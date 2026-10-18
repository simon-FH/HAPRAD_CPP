// weight_summary.C -- the RC weights of one workflow iteration, over all its
// friend files: how the selected events divide by rc_status, the quantiles of
// their RC factor, and the fractions far from 1 (<= 0, < 0.5, > 1.5, > 2, > 10)
// -- the share of the data that extreme grid nodes affect.
//
//   root -l -b -q 'weight_summary.C("<workflow dir>/iter_<k>")'
//
// One pass over the files. Every entry's `selected` flag is read, the RC
// factor and status only for selected entries, and the factors are counted in
// a fine histogram of log10(rc) rather than stored, so memory does not grow
// with the data (a full target is billions of entries).

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

  // ==== One pass: selected flag for every entry, the rest for selected ones ====
  Float_t rc = 0.f;
  Int_t status = 0;
  Bool_t sel = false;
  TBranch *bSel = nullptr, *bRc = nullptr, *bStatus = nullptr;
  Int_t tree = -1;
  Long64_t nSel = 0, byStatus[5] = {0, 0, 0, 0, 0}, nRc = 0, nonPos = 0, lt05 = 0, gt15 = 0, gt2 = 0, gt10 = 0;
  Double_t maxRc = 0.;
  TH1D logRc("logRc", "", 80000, -4., 4.);  // log10(rc) for rc > 0; quantiles to 1e-4 in log10
  for (Long64_t i = 0;; ++i) {
    const Long64_t local = c.LoadTree(i);
    if (local < 0) break;
    if (c.GetTreeNumber() != tree) {  // a new file: re-bind the three branches
      tree = c.GetTreeNumber();
      TTree* t = c.GetTree();
      bSel = t->GetBranch("selected");
      bRc = t->GetBranch("rc");
      bStatus = t->GetBranch("rc_status");
      bSel->SetAddress(&sel);
      bRc->SetAddress(&rc);
      bStatus->SetAddress(&status);
    }
    bSel->GetEntry(local);
    if (!sel) continue;
    ++nSel;
    bStatus->GetEntry(local);
    if (status >= 0 && status < 5) ++byStatus[status];
    if (status > 1) continue;  // no RC factor
    bRc->GetEntry(local);
    ++nRc;
    if (rc <= 0.f) {
      ++nonPos;
      continue;
    }
    logRc.Fill(std::log10(rc));
    lt05 += rc < 0.5f;
    gt15 += rc > 1.5f;
    gt2 += rc > 2.f;
    gt10 += rc > 10.f;
    maxRc = std::max(maxRc, Double_t(rc));
  }

  // ==== Report ====
  const char* label[5] = {"full grid support", "partial support", "outside the grid", "no grid support", "other particle/beam"};
  printf("%s: %d files, %lld selected events\n", iterDir, c.GetNtrees(), nSel);
  for (int s = 0; s < 5; ++s) printf("  rc_status %d %-20s %6.2f%%\n", s, label[s], 100. * byStatus[s] / std::max(nSel, 1LL));
  if (nRc == 0) return;
  const int nq = 6;
  Double_t p[nq] = {0.01, 0.05, 0.5, 0.95, 0.99, 0.999}, x[nq];
  logRc.GetQuantiles(nq, x, p);
  printf("RC factor of the %lld selected events with one:\n", nRc);
  printf("  1%% %.3f  5%% %.3f  median %.3f  95%% %.3f  99%% %.3f  99.9%% %.3f  max %.4g\n", std::pow(10., x[0]),
         std::pow(10., x[1]), std::pow(10., x[2]), std::pow(10., x[3]), std::pow(10., x[4]), std::pow(10., x[5]), maxRc);
  auto pc = [&](Long64_t k) { return 100. * k / nRc; };
  printf("  <= 0: %.3f%%   < 0.5: %.3f%%   > 1.5: %.3f%%   > 2: %.3f%%   > 10: %.3f%%\n", pc(nonPos), pc(lt05), pc(gt15),
         pc(gt2), pc(gt10));
}
