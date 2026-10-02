// ApplyRC -- radiative-correction weights for RG-E events (PLAN.md Phase 3).
//
//   ApplyRC <config> <grid.root> <ntuple.root> <friend.root>
//
// Writes a friend tree "RC" with one entry per entry of the input's `DT` tree,
// so any analysis binning can be applied afterwards:
//
//   TTree* dt = ...;  dt->AddFriend("RC", "friend.root");
//   dt->Draw("Q2", "RC.w * (RC.rc_status <= 1 && RC.selected)");
//
// Branches:
//   rc         sigma_obs / sigma_Born at the event, interpolated in the grid
//   rc_noex    the same without the exclusive radiative tail
//   w          1 / rc: the weight that turns an observed yield into a Born one
//              (0 when there is no rc)
//   rc_status  0 all grid corners had a HAPRAD result; 1 some did not, the
//              rest were used (see rc_cover); 2 outside the grid; 3 no corner
//              had a result; 4 not the configured particle or beam energy
//   rc_cover   interpolation weight that had a result (1 = full)
//   selected   the event passes the config's selection (the table's cuts)
//
// rc includes the exclusive tail (times Z/A, `target_naz`) only with
// `rc_exclusive_tail = yes`, which is refused for a table in arbitrary units:
// the tail would be on the wrong scale (PLAN.md Phase 1). The grid must come
// from the same config -- beam energy and particle are checked.
//
// One input file per call keeps the friend aligned entry by entry; run one job
// per file on the cluster.

#include "RCGrid.hxx"
#include "RGESelection.hxx"

#include "TChain.h"
#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <string>

int main(int argc, char** argv) {
  if (argc != 5) {
    fprintf(stderr, "usage: ApplyRC <config> <grid.root> <ntuple.root> <friend.root>\n");
    return 2;
  }
  rge::Config c;
  if (!rge::LoadConfig(argv[1], c)) return 1;
  rge::RCGrid grid;
  std::string why;
  if (!grid.Load(argv[2], why)) {
    fprintf(stderr, "%s\n", why.c_str());
    return 1;
  }
  if (std::fabs(std::atof(grid.Info("beam_energy").c_str()) - c.E) > 1e-6) {
    fprintf(stderr, "the grid is for E = %s GeV, the config for %g GeV\n", grid.Info("beam_energy").c_str(), c.E);
    return 1;
  }
  const bool absolute = rge::TableIsAbsolute(grid.Info("table_luminosity"));
  if (c.rcExclusiveTail && !absolute) {
    fprintf(stderr,
            "rc_exclusive_tail = yes, but the grid's table is in arbitrary units (luminosity: %s).\n"
            "The exclusive tail would be on the wrong scale. Build the table with a luminosity, or set "
            "rc_exclusive_tail = no.\n",
            grid.Info("table_luminosity").c_str());
    return 1;
  }
  const Double_t exScale = c.rcExclusiveTail ? c.targetNAZ : 0.;

  TChain chain("DT");
  if (chain.Add(argv[3]) != 1) {
    fprintf(stderr, "ApplyRC takes exactly one input file, so that the friend lines up with it\n");
    return 1;
  }
  rge::Event e;
  if (!rge::BindEvent(chain, e)) return 1;

  TFile out(argv[4], "RECREATE");
  if (out.IsZombie()) return 1;
  TTree tree("RC", "radiative-correction weights, friend of DT");
  Float_t rc, rcNoEx, w, cover;
  Int_t status;
  Bool_t selected;
  tree.Branch("rc", &rc, "rc/F");
  tree.Branch("rc_noex", &rcNoEx, "rc_noex/F");
  tree.Branch("w", &w, "w/F");
  tree.Branch("rc_status", &status, "rc_status/I");
  tree.Branch("rc_cover", &cover, "rc_cover/F");
  tree.Branch("selected", &selected, "selected/O");

  Long64_t count[5] = {0, 0, 0, 0, 0}, nSelected = 0, selectedByStatus[5] = {0, 0, 0, 0, 0};
  const Long64_t n = chain.GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    chain.GetEntry(i);
    rc = rcNoEx = w = cover = 0.f;
    selected = rge::Stages(c, e) == rge::kNStages;
    if (Int_t(e.pid) != c.pid || std::fabs(e.Eb - c.E) > c.Etol) {
      status = 4;
    } else {
      Double_t v[5], r[3], cv;
      rge::GridPoint(e, v);
      status = grid.Interpolate(v, r, cv);
      cover = cv;
      if (status <= rge::RCGrid::kPartial) {
        rcNoEx = r[0] + r[1];
        rc = r[0] + r[1] + exScale * r[2];
        w = rc > 0.f ? 1.f / rc : 0.f;
      }
    }
    ++count[status];
    if (selected) {
      ++nSelected;
      ++selectedByStatus[status];
    }
    tree.Fill();
  }
  tree.Write();
  TNamed("grid", argv[2]).Write();
  TNamed("table", grid.Info("table").c_str()).Write();
  TNamed("rc_definition", c.rcExclusiveTail ? Form("r_vr + r_in + %g r_ex (exclusive tail x Z/A)", exScale)
                                            : "r_vr + r_in (no exclusive tail)")
      .Write();
  TNamed("config", c.text.c_str()).Write();
  out.Close();

  const char* label[5] = {"full grid support", "partial support", "outside the grid", "no grid support",
                          "other particle / beam"};
  printf("%s: %lld entries, %lld selected\n", argv[4], n, nSelected);
  printf("  %-22s %12s %12s\n", "rc_status", "all", "selected");
  for (Int_t s = 0; s < 5; ++s) printf("  %d %-20s %12lld %12lld\n", s, label[s], count[s], selectedByStatus[s]);
  if (!absolute)
    printf("  NOTE: the table is in arbitrary units, so rc leaves out the exclusive tail (PLAN.md Phase 1).\n");
  return 0;
}
