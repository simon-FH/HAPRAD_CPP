// export_farm.C -- turn a workflow export (export_results.sh) into plain TSV
// tables for make_figures.py, which has no ROOT. For every iteration k:
//
//   nodes_<k>.tsv   every RC-grid node: Q2 nu z pt phi, valid, sigma_Born,
//                   r_vr r_in r_ex, and the fractions of table lookups that
//                   were outside the table, empty, or second-neighbour
//   cells_<k>.tsv   every table cell: centre, A Ac Acc (stored, reduced), the
//                   events in it, and its state: 0 fitted, 1 unfitted but
//                   filled from fitted neighbours (the second-neighbour
//                   copy), 2 unfitted and empty
//   logrc_<k>.tsv   the selected events' log10(rc) histogram, and rc <= 0
//   map_<v>_<k>.tsv the 2-D histograms of selected events with a factor, and
//                   of those with rc > 2, for v = Q2_y, z_pt2, phi_y
//
//   root -l -b -q 'export_farm.C("figures/farm/export", "figures/farm/tsv")'
//
// export_grid writes one grid's nodes alone, in the nodes_<k>.tsv format, e.g.
// the grid run_rc_grid.sh builds from the HAPRAD 2.0 model table:
//
//   root -l -b -q -e '.L export_farm.C' -e 'export_grid("tables/model_grid.root", "figures/farm/tsv/nodes_model.tsv")'

// Every RC-grid node: centre, valid flag, sigma_Born, r_vr r_in r_ex and the
// table-lookup fractions.
void export_grid(const char* gridPath, const char* outPath) {
  TFile g(gridPath);
  const char* names[8] = {"valid", "sigma_born", "r_vr", "r_in", "r_ex", "table_outside", "table_empty", "table_second"};
  THnD* h[8];
  for (int i = 0; i < 8; ++i) h[i] = g.Get<THnD>(names[i]);
  FILE* out = fopen(outPath, "w");
  fprintf(out, "Q2\tnu\tz\tpt\tphi\tvalid\tsigma_born\tr_vr\tr_in\tr_ex\toutside\tempty\tsecond\n");
  Int_t idx[5];
  for (Long64_t b = 0; b < h[0]->GetNbins(); ++b) {
    h[0]->GetBinContent(b, idx);
    bool inside = true;
    for (int d = 0; d < 5; ++d) inside = inside && idx[d] >= 1 && idx[d] <= h[0]->GetAxis(d)->GetNbins();
    if (!inside) continue;
    double c[5];
    for (int d = 0; d < 5; ++d) c[d] = h[0]->GetAxis(d)->GetBinCenter(idx[d]);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%d", c[0], c[1], c[2], c[3], c[4], h[0]->GetBinContent(b) > 0.5);
    for (int i = 1; i < 8; ++i) fprintf(out, "\t%.6g", h[i]->GetBinContent(b));
    fprintf(out, "\n");
  }
  fclose(out);
}

void export_farm(const char* inDir, const char* outDir) {
  gSystem->mkdir(outDir, true);
  for (int k = 0;; ++k) {
    TString gp = Form("%s/grid_%d.root", inDir, k);
    if (gSystem->AccessPathName(gp)) break;

    // ==== Grid nodes ====
    export_grid(gp, Form("%s/nodes_%d.tsv", outDir, k));

    // ==== Table cells, with the second-neighbour fill state ====
    {
      TFile t(Form("%s/table_%d.root", inDir, k));
      THnD *A = t.Get<THnD>("A"), *Ac = t.Get<THnD>("Ac"), *Acc = t.Get<THnD>("Acc"), *F = t.Get<THnD>("fitted"), *N = t.Get<THnD>("events");
      Int_t n[4];
      for (int d = 0; d < 4; ++d) n[d] = A->GetAxis(d)->GetNbins();
      // usable, as TSemiInclusiveModel uses it for a log table: fitted and A > 0
      auto usable = [&](Int_t* i) {
        const Long64_t b = A->GetBin(i);
        return F->GetBinContent(b) > 0 && A->GetBinContent(b) > 0;
      };
      FILE* out = fopen(Form("%s/cells_%d.tsv", outDir, k), "w");
      fprintf(out, "Q2\tnu\tz\tpt2\tA\tAc_red\tAcc_red\tevents\tstate\n");
      Int_t i[4], nb[4];
      for (i[0] = 1; i[0] <= n[0]; ++i[0])
        for (i[1] = 1; i[1] <= n[1]; ++i[1])
          for (i[2] = 1; i[2] <= n[2]; ++i[2])
            for (i[3] = 1; i[3] <= n[3]; ++i[3]) {
              int state = 2;
              if (usable(i))
                state = 0;
              else
                for (int m = 0; m < 81 && state == 2; ++m) {
                  int r = m, inside = 1, dist = 0;
                  for (int d = 0; d < 4; ++d) {
                    const int step = r % 3 - 1;
                    r /= 3;
                    nb[d] = i[d] + step;
                    dist += step * step;
                    inside = inside && nb[d] >= 1 && nb[d] <= n[d];
                  }
                  if (inside && dist > 0 && usable(nb)) state = 1;
                }
              const Long64_t b = A->GetBin(i);
              fprintf(out, "%g\t%g\t%g\t%g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\n", A->GetAxis(0)->GetBinCenter(i[0]),
                      A->GetAxis(1)->GetBinCenter(i[1]), A->GetAxis(2)->GetBinCenter(i[2]), A->GetAxis(3)->GetBinCenter(i[3]),
                      A->GetBinContent(b), Ac->GetBinContent(b), Acc->GetBinContent(b), N ? N->GetBinContent(b) : 0., state);
            }
      fclose(out);
    }
  }

  // ==== Event histograms ====
  TFile e(Form("%s/event_histograms.root", inDir));
  for (int k = 0;; ++k) {
    TH1D* h = e.Get<TH1D>(Form("h_logrc_%d", k));
    if (!h) break;
    FILE* out = fopen(Form("%s/logrc_%d.tsv", outDir, k), "w");
    fprintf(out, "# rc <= 0: %.0f events\nlog10rc_low\tlog10rc_high\tevents\n", e.Get<TH1D>(Form("n_le0_%d", k))->GetBinContent(1));
    for (int b = 0; b <= h->GetNbinsX() + 1; ++b)
      fprintf(out, "%g\t%g\t%.0f\n", b == 0 ? -1e9 : h->GetXaxis()->GetBinLowEdge(b),
              b > h->GetNbinsX() ? 1e9 : h->GetXaxis()->GetBinUpEdge(b), h->GetBinContent(b));
    fclose(out);
    for (const char* v : {"Q2_y", "z_pt2", "phi_y"}) {
      TH2D *a = e.Get<TH2D>(Form("all_%s_%d", v, k)), *g = e.Get<TH2D>(Form("gt2_%s_%d", v, k));
      FILE* o = fopen(Form("%s/map_%s_%d.tsv", outDir, v, k), "w");
      fprintf(o, "x_low\tx_high\ty_low\ty_high\tall\tgt2\n");
      for (int bx = 1; bx <= a->GetNbinsX(); ++bx)
        for (int by = 1; by <= a->GetNbinsY(); ++by)
          fprintf(o, "%g\t%g\t%g\t%g\t%.0f\t%.0f\n", a->GetXaxis()->GetBinLowEdge(bx), a->GetXaxis()->GetBinUpEdge(bx),
                  a->GetYaxis()->GetBinLowEdge(by), a->GetYaxis()->GetBinUpEdge(by), a->GetBinContent(bx, by),
                  g->GetBinContent(bx, by));
      fclose(o);
    }
  }
  printf("wrote %s\n", outDir);
}
