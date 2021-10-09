/*****************************************/
/*  FitPhiPQ.cxx                         */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program fits the PhiPQ distributions of a certain particle

// October 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main(int argc, char **argv) {

  gProgram = "FitPhiPQ";

  /*** INPUT ***/

  parseCommandLine(argc, argv);
  printOptions();

  TChain *dataChain = new TChain();
  if (gTargetOption == "D_C" || gTargetOption == "C") {
    dataChain->Add(gDataDir + "/C/*.root/ntuple_data");
  } else if (gTargetOption == "D_Fe" || gTargetOption == "Fe") {
    dataChain->Add(gDataDir + "/Fe/*.root/ntuple_data");
  } else if (gTargetOption == "D_Pb" || gTargetOption == "Pb") {
    dataChain->Add(gDataDir + "/Pb/*.root/ntuple_data");
  }

  TChain *simChain = new TChain();
  if (gTargetOption == "D_C" || gTargetOption == "D_Fe" || gTargetOption == "D_Pb") {
    simChain->Add(gSimDir + "/D/*.root/ntuple_sim");
  } else if (gTargetOption == "C") {
    simChain->Add(gSimDir + "/C/*.root/ntuple_sim");
  } else if (gTargetOption == "Fe") {
    simChain->Add(gSimDir + "/Fe/*.root/ntuple_sim");
  } else if (gTargetOption == "Pb") {
    simChain->Add(gSimDir + "/Pb/*.root/ntuple_sim");
  }

  // define cuts for data and sim. rec.
  TCut CutPID = Form("pid == %i", gPID);
  TCut CutDIS = "Q2 > 1. && W > 2. && Yb < 0.85";
  TCut CutVertex;
  if (gTargetOption == "D_C" || gTargetOption == "D_Fe" || gTargetOption == "D_Pb") {
    CutVertex = "TargType == 1 && vyec > -1.4 && vyec < 1.4";
  } else {  // in case of solid targets: C, Fe, Pb
    CutVertex = "TargType == 2 && vyec > -1.4 && vyec < 1.4";
  }

  // define cuts for MC
  TCut CutPID_MC = Form("mc_pid == %i", gPID);
  TCut CutDIS_MC = "mc_Q2 > 1. && mc_W > 2. && mc_Yb < 0.85";

  /*** GET NUMBER OF BINS ***/

  Int_t NbinsPhiPQ = (Int_t)(sizeof(kEdgesPhiPQ) / sizeof(kEdgesPhiPQ[0])) - 1;

  /*** HISTOGRAMS ***/

  // data histogram
  dataChain->Draw(Form("PhiPQ>>data(%i, -180., 180.)", NbinsPhiPQ), CutPID && CutDIS && CutVertex, "goff");
  TH1D *dataHist = (TH1D *)gROOT->FindObject("data");

  // MC gen. / thrown histogram
  simChain->Draw(Form("mc_PhiPQ>>mc(%i, -180., 180.)", NbinsPhiPQ), CutPID_MC && CutDIS_MC, "goff");
  TH1D *mcHist = (TH1D *)gROOT->FindObject("mc");

  // sim. rec. / accepted histogram
  simChain->Draw(Form("PhiPQ>>sim(%i, -180., 180.)", NbinsPhiPQ), CutPID && CutDIS, "goff");
  TH1D *simHist = (TH1D *)gROOT->FindObject("sim");

  // acceptance histogram
  TH1D *accHist = new TH1D("acc", "", NbinsPhiPQ, -180., 180.);
  accHist->Divide(simHist, mcHist, 1, 1, "B");

  // corrected data histogram
  TH1D *corrHist = new TH1D("corr", "", NbinsPhiPQ, -180., 180.);
  corrHist->Divide(dataHist, accHist, 1, 1);

  corrHist->SetTitle("");
  corrHist->GetYaxis()->SetTitle("Counts");
  corrHist->GetYaxis()->SetTitleSize(0.04);
  corrHist->GetYaxis()->SetTitleOffset(1.2);
  corrHist->GetYaxis()->SetMaxDigits(3);

  corrHist->GetXaxis()->SetTitle("#phi_{PQ} [deg]");
  corrHist->GetXaxis()->SetTitleSize(0.04);
  corrHist->GetXaxis()->SetTitleOffset(1.2);

  // style
  corrHist->SetMarkerColor(kBlack);
  corrHist->SetLineColor(kBlack);
  corrHist->SetLineWidth(2);
  corrHist->SetFillStyle(0);

  /*** CANVAS ***/

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(2);

  // set margin sizes
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);

  gStyle->SetOptStat(0);

  TString CanvasName = "PhiPQ_" + gTargetOption;
  TCanvas *can1 = new TCanvas(CanvasName, CanvasName, 1080, 1080);

  /*** FIT & DRAW ***/

  TF1 *func = new TF1("fit", "[0]+[1]*TMath::Cos(x*TMath::DegToRad())+[2]*TMath::Cos(2*x*TMath::DegToRad())");
  corrHist->Fit(func, "Q");
  corrHist->Draw("E");

  // pave
  TString Chi2String = Form("#chi^{2}/ndf = %.4f", func->GetChisquare() / (Double_t)func->GetNDF());
  TString A_String = Form("A = %.3f #pm %.3f", func->GetParameter(0), func->GetParError(0));
  TString Ac_String = Form("A_{c} = %.3f #pm %.6f", func->GetParameter(1), func->GetParError(1));
  TString Acc_String = Form("A_{cc} = %.3f #pm %.6f", func->GetParameter(2), func->GetParError(2));

  TPaveText *pav2 = new TPaveText(0.45, 0.7, 0.65, 0.9, "NDC NB");  // no border
  pav2->AddText("Fit on integrated " + gTargetOption + " Acc. Corr. Data");
  pav2->AddText(Chi2String);
  pav2->AddText(A_String);
  pav2->AddText(Ac_String);
  pav2->AddText(Acc_String);
  pav2->SetBorderSize(0);
  pav2->SetTextAlign(12);
  pav2->SetFillStyle(0);
  pav2->Draw();

  /*** OUTPUT ***/

  can1->Print(CanvasName + ".png");

  A_String = Form("export PARAM_A=%.3f", func->GetParameter(0));
  Ac_String = Form("export PARAM_AC=%.3f", func->GetParameter(1));
  Acc_String = Form("export PARAM_ACC=%.3f", func->GetParameter(2));

  TString OutFileName = "fit-results_" + gTargetOption + ".sh";
  std::ofstream OutFile(OutFileName, std::ios::out);

  OutFile << "#!/bin/bash" << std::endl;
  OutFile << std::endl;
  OutFile << A_String << std::endl;
  OutFile << Ac_String << std::endl;
  OutFile << Acc_String << std::endl;
  OutFile << std::endl;

  OutFile.close();

  std::cout << "The following file has been created: " << OutFileName << std::endl;
}
