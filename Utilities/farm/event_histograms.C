// event_histograms.C -- histograms of the selected events' RC factors over a
// whole workflow, every iteration, for plotting elsewhere: the distribution
// of the factor, and where in the kinematics the factor exceeds 2.
//
//   root -l -b -q 'event_histograms.C("<workflow dir>", "<out.root>")'
//
// It reads the workflow's own farm.cfg and runs.txt to find each run's
// ntuples, pairs every ntuple with its RC friend in each iter_<k>/friends/
// <run>/, and makes one pass per file: the `selected` flag of every entry
// (iteration 0's friend), and for selected entries the kinematics from DT and
// rc / rc_status from each iteration's friend. Written, per iteration k:
//
//   h_logrc_<k>          log10(rc) of selected events with a factor (rc > 0)
//   n_le0_<k>            how many had rc <= 0 (a 1-bin TH1D)
//   all_<v>_<k>, gt2_<v>_<k>
//                        selected events with a factor, and those with rc > 2,
//                        in the 2-D variables v = Q2_y, z_pt2, phi_y
//                        (phi in degrees, folded onto 0..180)

// The last `KEY=value` (quotes stripped) of a bash settings file.
TString SettingOf(const char* path, const char* key) {
  std::ifstream in(path);
  std::string line;
  TString value;
  const TString k = TString(key) + "=";
  while (std::getline(in, line)) {
    TString l(line.c_str());
    if (!l.BeginsWith(k)) continue;
    value = l(k.Length(), l.Length());
    const Ssiz_t hash = value.Index(" #");
    if (hash != kNPOS) value = value(0, hash);
    value = value.Strip(TString::kBoth);
    value.ReplaceAll("'", "");
    value.ReplaceAll("\"", "");
  }
  return value;
}

// The lines of `ls` on a pattern (shell globbing), or nothing.
std::vector<TString> Ls(const TString& pattern) {
  std::vector<TString> out;
  const TString list = gSystem->GetFromPipe("ls -d " + pattern + " 2>/dev/null");
  std::unique_ptr<TObjArray> t(list.Tokenize("\n"));
  for (TObject* o : *t) out.push_back(((TObjString*)o)->GetString());
  return out;
}

void event_histograms(const char* workflowDir, const char* outPath) {
  // ==== The workflow: settings, runs, iterations ====
  const TString W(workflowDir), cfg = W + "/farm.cfg";
  TString pattern = SettingOf(cfg, "DATA_PATTERN");
  pattern.ReplaceAll("%TARGET%", SettingOf(cfg, "TARGET"));
  pattern.ReplaceAll("%FIELD%", SettingOf(cfg, "FIELD"));
  pattern.ReplaceAll("%TRACKING%", SettingOf(cfg, "TRACKING"));
  std::vector<TString> runs;
  {
    std::ifstream in((W + "/runs.txt").Data());
    std::string r;
    while (in >> r) runs.push_back(r.c_str());
  }
  int nIter = 0;
  while (!gSystem->AccessPathName(Form("%s/iter_%d/friends", workflowDir, nIter))) ++nIter;
  printf("%s: %zu runs, %d iterations, inputs %s\n", workflowDir, runs.size(), nIter, pattern.Data());
  if (runs.empty() || nIter == 0) return;

  // ==== Histograms ====
  TFile out(outPath, "RECREATE");
  std::vector<TH1D*> hLog, hLe0;
  std::vector<std::array<TH2D*, 3>> hAll, hGt2;
  for (int k = 0; k < nIter; ++k) {
    hLog.push_back(new TH1D(Form("h_logrc_%d", k), "log10(rc);log10(rc);events", 1600, -2., 4.));
    hLe0.push_back(new TH1D(Form("n_le0_%d", k), "events with rc <= 0", 1, 0., 1.));
    std::array<TH2D*, 3> a, g;
    for (int s = 0; s < 2; ++s) {
      auto& h = s == 0 ? a : g;
      const char* p = s == 0 ? "all" : "gt2";
      h[0] = new TH2D(Form("%s_Q2_y_%d", p, k), ";y;Q^{2} [GeV^{2}]", 26, 0.2, 0.85, 36, 1., 10.);
      h[1] = new TH2D(Form("%s_z_pt2_%d", p, k), ";p_{T}^{2} [GeV^{2}];z", 30, 0., 1.5, 34, 0.15, 1.);
      h[2] = new TH2D(Form("%s_phi_y_%d", p, k), ";y;#phi_{h} [deg]", 26, 0.2, 0.85, 18, 0., 180.);
    }
    hAll.push_back(a);
    hGt2.push_back(g);
  }

  // ==== One pass per input file ====
  Long64_t nFiles = 0, nSel = 0;
  for (const TString& run : runs) {
    TString pat = pattern;
    pat.ReplaceAll("%RUN%", run);
    pat.ReplaceAll("%RUNNO%", Form("%d", run.Atoi()));
    for (const TString& input : Ls(pat)) {
      TString base = gSystem->BaseName(input);
      base.ReplaceAll(".root", "");
      std::unique_ptr<TFile> fd(TFile::Open(input));
      TTree* dt = fd ? fd->Get<TTree>("DT") : nullptr;
      if (!dt) {
        printf("  %s: no DT tree, skipped\n", input.Data());
        continue;
      }
      // The friends of this file in every iteration; all must exist and line up.
      std::vector<std::unique_ptr<TFile>> ff;
      std::vector<TTree*> rcT;
      for (int k = 0; k < nIter; ++k) {
        ff.emplace_back(TFile::Open(Form("%s/iter_%d/friends/%s/%s_rc.root", workflowDir, k, run.Data(), base.Data())));
        TTree* t = ff.back() ? ff.back()->Get<TTree>("RC") : nullptr;
        if (!t || t->GetEntries() != dt->GetEntries()) {
          rcT.clear();
          break;
        }
        rcT.push_back(t);
      }
      if ((int)rcT.size() != nIter) {
        printf("  %s: friends missing or not aligned, skipped\n", input.Data());
        continue;
      }
      Float_t Q2, yb, zh, pt2, phi;
      dt->SetBranchStatus("*", 0);
      const char* names[5] = {"Q2", "y_bjorken", "z_h", "p_T2", "phi_PQ"};
      Float_t* addr[5] = {&Q2, &yb, &zh, &pt2, &phi};
      std::vector<TBranch*> bk;
      for (int i = 0; i < 5; ++i) {
        dt->SetBranchStatus(names[i], 1);
        dt->SetBranchAddress(names[i], addr[i]);
        bk.push_back(dt->GetBranch(names[i]));
      }
      std::vector<Float_t> rc(nIter);
      std::vector<Int_t> st(nIter);
      Bool_t sel = false;
      TBranch* bSel = rcT[0]->GetBranch("selected");
      bSel->SetAddress(&sel);
      std::vector<TBranch*> bRc, bSt;
      for (int k = 0; k < nIter; ++k) {
        bRc.push_back(rcT[k]->GetBranch("rc"));
        bSt.push_back(rcT[k]->GetBranch("rc_status"));
        bRc[k]->SetAddress(&rc[k]);
        bSt[k]->SetAddress(&st[k]);
      }
      const Long64_t n = dt->GetEntries();
      for (Long64_t i = 0; i < n; ++i) {
        bSel->GetEntry(i);
        if (!sel) continue;
        ++nSel;
        for (TBranch* b : bk) b->GetEntry(i);
        const double phiDeg = std::fabs(std::remainder(double(phi), 2. * TMath::Pi())) * 180. / TMath::Pi();
        for (int k = 0; k < nIter; ++k) {
          bSt[k]->GetEntry(i);
          if (st[k] > 1) continue;  // no factor
          bRc[k]->GetEntry(i);
          if (rc[k] <= 0.f) {
            hLe0[k]->Fill(0.5);
            continue;
          }
          hLog[k]->Fill(std::log10(rc[k]));
          const double v[3][2] = {{yb, Q2}, {pt2, zh}, {yb, phiDeg}};
          for (int s = 0; s < 3; ++s) {
            hAll[k][s]->Fill(v[s][0], v[s][1]);
            if (rc[k] > 2.f) hGt2[k][s]->Fill(v[s][0], v[s][1]);
          }
        }
      }
      ++nFiles;
      printf("  %s: %lld entries\n", base.Data(), n);
    }
  }

  // ==== Write ====
  out.cd();
  TNamed("source", Form("%s: %lld files, %lld selected events, %d iterations", workflowDir, nFiles, nSel, nIter)).Write();
  out.Write();
  printf("wrote %s: %lld files, %lld selected events\n", outPath, nFiles, nSel);
}
