// table_shapes.C -- the phi shape of every fitted cell of a structure-function
// table, as a TSV for plotting: the cell centre, A, the modulations
// rc = Ac p_t / A and rcc = Acc p_t^2 / A at the centre (a reduced table
// stores Ac / p_t and Acc / p_t^2), the minimum over phi of
// 1 + rc cos(phi) + rcc cos(2 phi) -- the Born cross section's phi shape,
// which must stay positive -- and how many standard deviations the
// modulations are from zero.
//
//   root -l -b -q 'table_shapes.C("table.root", "shapes.tsv")'

void table_shapes(const char* tablePath, const char* outPath) {
  TFile f(tablePath);
  THnD *A = f.Get<THnD>("A"), *Ac = f.Get<THnD>("Ac"), *Acc = f.Get<THnD>("Acc"), *F = f.Get<THnD>("fitted");
  THnD *eAc = f.Get<THnD>("Ac_err"), *eAcc = f.Get<THnD>("Acc_err");
  const TNamed* scaling = f.Get<TNamed>("pt_scaling");
  const bool reduced = scaling && TString(scaling->GetTitle()) == "reduced";
  FILE* out = fopen(outPath, "w");
  fprintf(out, "# %s\nQ2\tnu\tz\tpt2\tA\trc\trcc\tmin_shape\tsignificance\n", tablePath);
  Int_t idx[4];
  for (Long64_t b = 0; b < A->GetNbins(); ++b) {
    if (F ? F->GetBinContent(b) == 0 : A->GetBinContent(b) == 0) continue;
    A->GetBinContent(b, idx);
    Double_t c[4];
    for (int d = 0; d < 4; ++d) c[d] = A->GetAxis(d)->GetBinCenter(idx[d]);
    const double pt2 = c[3], pt = std::sqrt(pt2), a = A->GetBinContent(b);
    const double rc = Ac->GetBinContent(b) * (reduced ? pt : 1.) / a, rcc = Acc->GetBinContent(b) * (reduced ? pt2 : 1.) / a;
    double m = 1e9;
    for (int k = 0; k <= 360; ++k) {
      const double p = k * TMath::Pi() / 180.;
      m = std::min(m, 1. + rc * std::cos(p) + rcc * std::cos(2. * p));
    }
    // Significance needs the fit errors, which only data tables carry.
    double sig = -1.;
    if (eAc && eAcc && eAc->GetBinContent(b) > 0 && eAcc->GetBinContent(b) > 0)
      sig = std::hypot(Ac->GetBinContent(b) / eAc->GetBinContent(b), Acc->GetBinContent(b) / eAcc->GetBinContent(b));
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", c[0], c[1], c[2], c[3], a, rc, rcc, m, sig);
  }
  fclose(out);
}
