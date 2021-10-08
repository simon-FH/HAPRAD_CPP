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

  TChain *DataChain = new TChain();
  if (gTargetOption == "D_C" || gTargetOption == "C") {
    DataChain->Add(gDataDir + "/C/*.root/ntuple_data");
  } else if (gTargetOption == "D_Fe" || gTargetOption == "Fe") {
    DataChain->Add(gDataDir + "/Fe/*.root/ntuple_data");
  } else if (gTargetOption == "D_Pb" || gTargetOption == "Pb") {
    DataChain->Add(gDataDir + "/Pb/*.root/ntuple_data");
  }

  // define cuts
  TCut CutPID = Form("pid == %i", gPID);
  TCut CutDIS = "Q2 > 1. && W > 2. && Yb < 0.85";
  TCut CutVertex;
  if (gTargetOption == "D_C" || gTargetOption == "D_Fe" || gTargetOption == "D_Pb") {
    CutVertex = "TargType == 1 && vyec > -1.4 && vyec < 1.4";
  } else {  // in case of solid targets: C, Fe, Pb
    CutVertex = "TargType == 2 && vyec > -1.4 && vyec < 1.4";
  }

  Int_t NbinsPhiPQ = (Int_t)(sizeof(kEdgesPhiPQ) / sizeof(kEdgesPhiPQ[0]));

  TH1D *theHist1;
  DataChain->Draw(Form("PhiPQ>>data(%i, -180, 180)", NbinsPhiPQ), CutPID && CutDIS && CutVertex, "goff");
  theHist1 = (TH1D *)gROOT->FindObject("data");

  theHist1->SetTitle("");
  theHist1->GetYaxis()->SetTitle("Counts");
  theHist1->GetYaxis()->SetTitleSize(0.04);
  theHist1->GetYaxis()->SetTitleOffset(1.2);
  theHist1->GetYaxis()->SetMaxDigits(3);

  theHist1->GetXaxis()->SetTitle("#phi_{PQ} [deg]");
  theHist1->GetXaxis()->SetTitleSize(0.04);
  theHist1->GetXaxis()->SetTitleOffset(1.2);

  // style
  theHist1->SetMarkerColor(kBlack);
  theHist1->SetLineColor(kBlack);
  theHist1->SetLineWidth(2);
  theHist1->SetFillStyle(0);

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
  theHist1->Fit(func, "Q");
  theHist1->Draw("E");

  // pave
  TString Chi2String = Form("#chi^{2}/ndf = %.4f", func->GetChisquare() / (Double_t)func->GetNDF());
  TString A_String = Form("A = %.3f #pm %.3f", func->GetParameter(0), func->GetParError(0));
  TString Ac_String = Form("A_{c} = %.3f #pm %.6f", func->GetParameter(1), func->GetParError(1));
  TString Acc_String = Form("A_{cc} = %.3f #pm %.6f", func->GetParameter(2), func->GetParError(2));

  TPaveText *pav2 = new TPaveText(0.45, 0.7, 0.65, 0.9, "NDC NB");  // no border
  pav2->AddText("Fit on integrated " + gTargetOption + " data");
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
