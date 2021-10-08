/*****************************************/
/*  GetCentroids.cxx                     */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program gets the centroids of (PhiPQ) bins based on a (Q2, Nu, Zh, Pt2, PhiPQ) binning given by a CSV file
// June 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main(int argc, char **argv) {

  gProgram = "GetCentroids";

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

  /*** READ CSV FILE ***/

  std::ifstream BinningFile;
  BinningFile.open("binning.csv", std::ios::in);
  std::vector<Double_t> Q2_low;
  std::vector<Double_t> Q2_up;
  std::vector<Double_t> Nu_low;
  std::vector<Double_t> Nu_up;
  std::vector<Double_t> Z_low;
  std::vector<Double_t> Z_up;
  std::vector<Double_t> Pt2_low;
  std::vector<Double_t> Pt2_up;
  std::vector<Double_t> PhiPQ_low;
  std::vector<Double_t> PhiPQ_up;

  // exit program if ifstream could not open file
  if (!BinningFile) {
    std::cerr << "ERROR: File could not be opened" << std::endl;
    exit(EXIT_FAILURE);
  }

  TString auxLine;
  TString auxString[10];
  while (BinningFile >> auxLine) {
    for (Int_t i = 0; i < 10; i++) {
      auxString[i] = ((TObjString *)(auxLine.Tokenize(",")->At(i)))->String();
    }
    Q2_low.push_back(auxString[0].Atof());
    Q2_up.push_back(auxString[1].Atof());
    Nu_low.push_back(auxString[2].Atof());
    Nu_up.push_back(auxString[3].Atof());
    Z_low.push_back(auxString[4].Atof());
    Z_up.push_back(auxString[5].Atof());
    Pt2_low.push_back(auxString[6].Atof());
    Pt2_up.push_back(auxString[7].Atof());
    PhiPQ_low.push_back(auxString[8].Atof());
    PhiPQ_up.push_back(auxString[9].Atof());
  }

  /*** HISTOGRAMS ***/

  Int_t NTotalBins = (Int_t)PhiPQ_low.size();

  Double_t meanQ2[NTotalBins];
  TH1D *histQ2[NTotalBins];
  Double_t meanNu[NTotalBins];
  TH1D *histNu[NTotalBins];
  Double_t meanZ[NTotalBins];
  TH1D *histZ[NTotalBins];
  Double_t meanPt2[NTotalBins];
  TH1D *histPt2[NTotalBins];
  Double_t meanPhiPQ[NTotalBins];
  TH1D *histPhiPQ[NTotalBins];

  TCut CutBin;

  for (Int_t i = 0; i < NTotalBins; i++) {
    CutBin = Form(
        "Q2 > %.2f && Q2 < %.2f && Nu > %.2f && Nu < %.2f && Zh > %.2f && Zh < %.2f && Pt2 > %.2f && Pt2 < %.2f && PhiPQ > %.2f && "
        "PhiPQ < %.2f",
        Q2_low[i], Q2_up[i], Nu_low[i], Nu_up[i], Z_low[i], Z_up[i], Pt2_low[i], Pt2_up[i], PhiPQ_low[i], PhiPQ_up[i]);
    // Q2
    DataChain->Draw(Form("Q2>>h_Q2_%i(100, 1.0, 4.0)", i), CutPID && CutDIS && CutVertex && CutBin, "goff");
    histQ2[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("h_Q2_%i", i));
    meanQ2[i] = histQ2[i]->GetMean();
    // Nu
    DataChain->Draw(Form("Nu>>h_Nu_%i(100, 2.2, 4.2)", i), CutPID && CutDIS && CutVertex && CutBin, "goff");
    histNu[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("h_Nu_%i", i));
    meanNu[i] = histNu[i]->GetMean();
    // Z
    DataChain->Draw(Form("Zh>>h_Z_%i(100, 0.5, 0.9)", i), CutPID && CutDIS && CutVertex && CutBin, "goff");
    histZ[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("h_Z_%i", i));
    meanZ[i] = histZ[i]->GetMean();
    // Pt2
    DataChain->Draw(Form("Pt2>>h_Pt2_%i(100, 0.0, 1.5)", i), CutPID && CutDIS && CutVertex && CutBin, "goff");
    histPt2[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("h_Pt2_%i", i));
    meanPt2[i] = histPt2[i]->GetMean();
    // PhiPQ
    DataChain->Draw(Form("PhiPQ>>h_PhiPQ_%i(36, -180, 180)", i), CutPID && CutDIS && CutVertex && CutBin, "goff");
    histPhiPQ[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("h_PhiPQ_%i", i));
    meanPhiPQ[i] = histPhiPQ[i]->GetMean();
  }

  /*** OUTPUT FILE ***/

  TString OutFileName = "centroids_" + gTargetOption + ".csv";
  std::ofstream OutFile(OutFileName, std::ios::out);
  for (Int_t i = 0; i < NTotalBins; i++) {
    OutFile << meanQ2[i] << "," << meanNu[i] << "," << meanZ[i] << "," << meanPt2[i] << "," << meanPhiPQ[i] << std::endl;
  }

  OutFile.close();
  std::cout << "The following file has been created: " << OutFileName << std::endl;

  return 0;
}
