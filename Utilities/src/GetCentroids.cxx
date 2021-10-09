/*****************************************/
/*  GetCentroids.cxx                     */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program gets the centroids of (PhiPQ) bins based on a (Q2, Nu, Zh, Pt2, PhiPQ) binning given by a CSV file

// October 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main(int argc, char **argv) {

  gProgram = "GetCentroids";

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

  /*** GET NUMBER OF BINS ***/

  Int_t NbinsQ2 = (Int_t)(sizeof(kEdgesQ2) / sizeof(kEdgesQ2[0])) - 1;
  Int_t NbinsNu = (Int_t)(sizeof(kEdgesNu) / sizeof(kEdgesNu[0])) - 1;
  Int_t NbinsZ = (Int_t)(sizeof(kEdgesZ) / sizeof(kEdgesZ[0])) - 1;
  Int_t NbinsPt2 = (Int_t)(sizeof(kEdgesPt2) / sizeof(kEdgesPt2[0])) - 1;
  Int_t NbinsPhiPQ = (Int_t)(sizeof(kEdgesPhiPQ) / sizeof(kEdgesPhiPQ[0])) - 1;

  Int_t NTotalBins = (Int_t)PhiPQ_low.size(); // = number of lines of the file

  /*** HISTOGRAMS ***/

  // define data histograms
  TH1D *histQ2_Data[NTotalBins];
  TH1D *histNu_Data[NTotalBins];
  TH1D *histZ_Data[NTotalBins];
  TH1D *histPt2_Data[NTotalBins];
  TH1D *histPhiPQ_Data[NTotalBins];

  // define MC generated/thrown histograms
  TH1D *histQ2_MC[NTotalBins];
  TH1D *histNu_MC[NTotalBins];
  TH1D *histZ_MC[NTotalBins];
  TH1D *histPt2_MC[NTotalBins];
  TH1D *histPhiPQ_MC[NTotalBins];

  // define simulations reconstructed/accepted histograms
  TH1D *histQ2_Sim[NTotalBins];
  TH1D *histNu_Sim[NTotalBins];
  TH1D *histZ_Sim[NTotalBins];
  TH1D *histPt2_Sim[NTotalBins];
  TH1D *histPhiPQ_Sim[NTotalBins];

  // define acceptance histograms
  TH1D *histQ2_Acceptance[NTotalBins];
  TH1D *histNu_Acceptance[NTotalBins];
  TH1D *histZ_Acceptance[NTotalBins];
  TH1D *histPt2_Acceptance[NTotalBins];
  TH1D *histPhiPQ_Acceptance[NTotalBins];

  // define corrected data histograms
  TH1D *histQ2_Corr[NTotalBins];
  TH1D *histNu_Corr[NTotalBins];
  TH1D *histZ_Corr[NTotalBins];
  TH1D *histPt2_Corr[NTotalBins];
  TH1D *histPhiPQ_Corr[NTotalBins];

  Double_t meanQ2[NTotalBins];
  Double_t meanNu[NTotalBins];
  Double_t meanZ[NTotalBins];
  Double_t meanPt2[NTotalBins];
  Double_t meanPhiPQ[NTotalBins];

  TCut CutBin;
  TCut CutBin_MC;

  for (Int_t i = 0; i < NTotalBins; i++) {

    // set bin for data and sim. rec.
    CutBin = Form(
        "Q2 > %.2f && Q2 < %.2f && Nu > %.2f && Nu < %.2f && Zh > %.2f && Zh < %.2f && Pt2 > %.2f && Pt2 < %.2f && PhiPQ > %.2f && "
        "PhiPQ < %.2f",
        Q2_low[i], Q2_up[i], Nu_low[i], Nu_up[i], Z_low[i], Z_up[i], Pt2_low[i], Pt2_up[i], PhiPQ_low[i], PhiPQ_up[i]);

    // set bin edges for mc
    CutBin_MC = Form(
        "mc_Q2 > %.2f && mc_Q2 < %.2f && mc_Nu > %.2f && mc_Nu < %.2f && mc_Zh > %.2f && mc_Zh < %.2f && mc_Pt2 > %.2f && mc_Pt2 < %.2f && "
        "mc_PhiPQ > %.2f && "
        "mc_PhiPQ < %.2f",
        Q2_low[i], Q2_up[i], Nu_low[i], Nu_up[i], Z_low[i], Z_up[i], Pt2_low[i], Pt2_up[i], PhiPQ_low[i], PhiPQ_up[i]);

    /*** DATA HISTOGRAMS ***/

    // Q2
    dataChain->Draw(Form("Q2>>data_Q2_%i(100, %.2f, %.2f)", i, kEdgesQ2[0], kEdgesQ2[NbinsQ2]), CutPID && CutDIS && CutVertex && CutBin,
                    "goff");
    histQ2_Data[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("data_Q2_%i", i));
    // Nu
    dataChain->Draw(Form("Nu>>data_Nu_%i(100, %.2f, %.2f)", i, kEdgesNu[0], kEdgesNu[NbinsNu]), CutPID && CutDIS && CutVertex && CutBin,
                    "goff");
    histNu_Data[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("data_Nu_%i", i));
    // Z
    dataChain->Draw(Form("Zh>>data_Z_%i(100, %.2f, %.2f)", i, kEdgesZ[0], kEdgesZ[NbinsZ]), CutPID && CutDIS && CutVertex && CutBin,
                    "goff");
    histZ_Data[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("data_Z_%i", i));
    // Pt2
    dataChain->Draw(Form("Pt2>>data_Pt2_%i(100, %.2f, %.2f)", i, kEdgesPt2[0], kEdgesPt2[NbinsPt2]),
                    CutPID && CutDIS && CutVertex && CutBin, "goff");
    histPt2_Data[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("data_Pt2_%i", i));
    // PhiPQ
    dataChain->Draw(Form("PhiPQ>>data_PhiPQ_%i(36, -180., 180.)", i), CutPID && CutDIS && CutVertex && CutBin, "goff");
    histPhiPQ_Data[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("data_PhiPQ_%i", i));

    /*** MC GEN. HISTOGRAMS ***/

    // Q2
    simChain->Draw(Form("mc_Q2>>mc_Q2_%i(100, %.2f, %.2f)", i, kEdgesQ2[0], kEdgesQ2[NbinsQ2]), CutPID_MC && CutDIS_MC && CutBin_MC,
                   "goff");
    histQ2_MC[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("mc_Q2_%i", i));
    // Nu
    simChain->Draw(Form("mc_Nu>>mc_Nu_%i(100, %.2f, %.2f)", i, kEdgesNu[0], kEdgesNu[NbinsNu]), CutPID_MC && CutDIS_MC && CutBin_MC,
                   "goff");
    histNu_MC[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("mc_Nu_%i", i));
    // Z
    simChain->Draw(Form("mc_Zh>>mc_Z_%i(100, %.2f, %.2f)", i, kEdgesZ[0], kEdgesZ[NbinsZ]), CutPID_MC && CutDIS_MC && CutBin_MC,
                   "goff");
    histZ_MC[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("mc_Z_%i", i));
    // Pt2
    simChain->Draw(Form("mc_Pt2>>mc_Pt2_%i(100, %.2f, %.2f)", i, kEdgesPt2[0], kEdgesPt2[NbinsPt2]),
                   CutPID_MC && CutDIS_MC && CutBin_MC, "goff");
    histPt2_MC[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("mc_Pt2_%i", i));
    // PhiPQ
    simChain->Draw(Form("mc_PhiPQ>>mc_PhiPQ_%i(36, -180., 180.)", i), CutPID_MC && CutDIS_MC && CutBin_MC, "goff");
    histPhiPQ_MC[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("mc_PhiPQ_%i", i));

    /*** SIM. REC. HISTOGRAMS ***/

    // Q2
    simChain->Draw(Form("Q2>>sim_Q2_%i(100, %.2f, %.2f)", i, kEdgesQ2[0], kEdgesQ2[NbinsQ2]), CutPID && CutDIS && CutBin, "goff");
    histQ2_Sim[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("sim_Q2_%i", i));
    // Nu
    simChain->Draw(Form("Nu>>sim_Nu_%i(100, %.2f, %.2f)", i, kEdgesNu[0], kEdgesNu[NbinsNu]), CutPID && CutDIS && CutBin, "goff");
    histNu_Sim[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("sim_Nu_%i", i));
    // Z
    simChain->Draw(Form("Zh>>sim_Z_%i(100, %.2f, %.2f)", i, kEdgesZ[0], kEdgesZ[NbinsZ]), CutPID && CutDIS && CutBin, "goff");
    histZ_Sim[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("sim_Z_%i", i));
    // Pt2
    simChain->Draw(Form("Pt2>>sim_Pt2_%i(100, %.2f, %.2f)", i, kEdgesPt2[0], kEdgesPt2[NbinsPt2]), CutPID && CutDIS && CutBin, "goff");
    histPt2_Sim[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("sim_Pt2_%i", i));
    // PhiPQ
    simChain->Draw(Form("PhiPQ>>sim_PhiPQ_%i(36, -180., 180.)", i), CutPID && CutDIS && CutBin, "goff");
    histPhiPQ_Sim[i] = (TH1D *)gDirectory->GetList()->FindObject(Form("sim_PhiPQ_%i", i));

    /*** CALCULATE ACCEPTANCE ***/

    // Q2
    histQ2_Acceptance[i] = new TH1D(Form("acc_Q2_%i", i), "", 100, kEdgesQ2[0], kEdgesQ2[NbinsQ2]);
    histQ2_Acceptance[i]->Divide(histQ2_Sim[i], histQ2_MC[i], 1, 1, "B");
    // Nu
    histNu_Acceptance[i] = new TH1D(Form("acc_Nu_%i", i), "", 100, kEdgesNu[0], kEdgesNu[NbinsNu]);
    histNu_Acceptance[i]->Divide(histNu_Sim[i], histNu_MC[i], 1, 1, "B");
    // Z
    histZ_Acceptance[i] = new TH1D(Form("acc_Z_%i", i), "", 100, kEdgesZ[0], kEdgesZ[NbinsZ]);
    histZ_Acceptance[i]->Divide(histZ_Sim[i], histZ_MC[i], 1, 1, "B");
    // Pt2
    histPt2_Acceptance[i] = new TH1D(Form("acc_Pt2_%i", i), "", 100, kEdgesPt2[0], kEdgesPt2[NbinsPt2]);
    histPt2_Acceptance[i]->Divide(histPt2_Sim[i], histPt2_MC[i], 1, 1, "B");
    // PhiPQ
    histPhiPQ_Acceptance[i] = new TH1D(Form("acc_PhiPQ_%i", i), "", 36, kEdgesPhiPQ[0], kEdgesPhiPQ[NbinsPhiPQ]);
    histPhiPQ_Acceptance[i]->Divide(histPhiPQ_Sim[i], histPhiPQ_MC[i], 1, 1, "B");

    /*** CORRECT DATA ***/

    // Q2
    histQ2_Corr[i] = new TH1D(Form("corr_Q2_%i", i), "", 100, kEdgesQ2[0], kEdgesQ2[NbinsQ2]);
    histQ2_Corr[i]->Divide(histQ2_Data[i], histQ2_Acceptance[i], 1, 1);
    // Nu
    histNu_Corr[i] = new TH1D(Form("corr_Nu_%i", i), "", 100, kEdgesNu[0], kEdgesNu[NbinsNu]);
    histNu_Corr[i]->Divide(histNu_Data[i], histNu_Acceptance[i], 1, 1);
    // Z
    histZ_Corr[i] = new TH1D(Form("corr_Z_%i", i), "", 100, kEdgesZ[0], kEdgesZ[NbinsZ]);
    histZ_Corr[i]->Divide(histZ_Data[i], histZ_Acceptance[i], 1, 1);
    // Pt2
    histPt2_Corr[i] = new TH1D(Form("corr_Pt2_%i", i), "", 100, kEdgesPt2[0], kEdgesPt2[NbinsPt2]);
    histPt2_Corr[i]->Divide(histPt2_Data[i], histPt2_Acceptance[i], 1, 1);
    // PhiPQ
    histPhiPQ_Corr[i] = new TH1D(Form("corr_PhiPQ_%i", i), "", 36, kEdgesPhiPQ[0], kEdgesPhiPQ[NbinsPhiPQ]);
    histPhiPQ_Corr[i]->Divide(histPhiPQ_Data[i], histPhiPQ_Acceptance[i], 1, 1);

    /*** GET CENTROIDS ***/

    meanQ2[i] = histQ2_Corr[i]->GetMean();
    meanNu[i] = histNu_Corr[i]->GetMean();
    meanZ[i] = histZ_Corr[i]->GetMean();
    meanPt2[i] = histPt2_Corr[i]->GetMean();
    meanPhiPQ[i] = histPhiPQ_Corr[i]->GetMean();
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
