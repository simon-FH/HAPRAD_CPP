/*****************************************/
/*  GetRC.cxx                            */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// October 2021

#include "Headers.hxx"
#include "UX.hxx"

int main(int argc, char **argv) {

  gProgram = "GetRC";

  /*** INPUT ***/

  parseCommandLine(argc, argv);
  printOptions();

  Double_t NAZ = 0.5;
  if (gTargetOption == "Pb") {
    NAZ = 82. / 208.;
  }

  /*** READ CSV FILE ***/

  std::vector<Double_t> q2_centroid;
  std::vector<Double_t> nu_centroid;
  std::vector<Double_t> zh_centroid;
  std::vector<Double_t> pt2_centroid;
  std::vector<Double_t> phi_centroid;

  std::ifstream CentroidsFile;
  CentroidsFile.open("centroids_" + gTargetOption + ".csv", std::ios::in);

  // exit program if ifstream could not open file
  if (!CentroidsFile) {
    std::cerr << "ERROR: File could not be opened" << std::endl;
    exit(EXIT_FAILURE);
  }

  TString auxLine;
  TString auxString[5];
  while (CentroidsFile >> auxLine) {
    for (Int_t i = 0; i < 5; i++) {
      auxString[i] = ((TObjString *)(auxLine.Tokenize(",")->At(i)))->String();
    }
    q2_centroid.push_back(auxString[0].Atof());
    nu_centroid.push_back(auxString[1].Atof());
    zh_centroid.push_back(auxString[2].Atof());
    pt2_centroid.push_back(auxString[3].Atof());
    phi_centroid.push_back(auxString[4].Atof());
  }

  /*** RECALCULATE SOME CENTROIDS ***/

  Double_t MassProton = 0.938272;
  std::vector<Double_t> xb_centroid;
  std::vector<Double_t> pt_centroid;
  for (Int_t i = 0; i < (Int_t)phi_centroid.size(); i++) {
    xb_centroid.push_back(q2_centroid[i] / (2 * MassProton * nu_centroid[i]));
    pt_centroid.push_back(TMath::Sqrt(pt2_centroid[i]));
  }

  /*** CALCULATE RC FACTORS ***/

  // this "m" is necessary for the RC calculations
  // the model is prepared for pions, that explains the kMassPion term here
  Double_t m = TMath::Power((kMassNeutron + kMassPion), 2);
  Double_t f1, f3;
  Double_t a1, a2, a3;

  TRadCor rc;

  a1 = 1;
  a2 = 0;
  a3 = a1 / a2;

  // define and create output file
  std::ofstream out;
  out.open("RCFactor_" + gTargetOption + ".txt");

  // write first line
  out << std::setw(15) << std::fixed << std::setprecision(5) << "XbCentroid" << std::setw(15) << "Q2Centroid" << std::setw(15)
      << "ZhCentroid" << std::setw(15) << "PtCentroid" << std::setw(15) << "PhiPQCentroid" << std::setw(15) << "Fact_noex" << std::setw(15)
      << "Fact_ex" << std::endl;

  // loop over bins
  for (Int_t i = 0; i < (Int_t)phi_centroid.size(); i++) {
    // calculate rc factors
    rc.CalculateRCFactor(5.015, xb_centroid[i], q2_centroid[i], zh_centroid[i], pt_centroid[i], phi_centroid[i], m, NAZ);
    f1 = rc.GetFactor1();
    f3 = rc.GetFactor3();
    if (TMath::IsNaN(f1) || f1 == a3) f1 = 0;
    if (TMath::IsNaN(f3) || f3 == a3) f3 = 0;
    // if any of the centroids is zero, null the rc factor
    Bool_t emptyBin = xb_centroid[i] == 0. || q2_centroid[i] == 0. || zh_centroid[i] == 0. || pt_centroid[i] == 0. || phi_centroid[i] == 0.;
    if (emptyBin) {
      f1 = 0;
      f3 = 0;
    }
    // write centroids and rc factors into output file
    out << std::setw(15) << std::fixed << std::setprecision(5) << xb_centroid[i] << std::setw(15) << q2_centroid[i] << std::setw(15)
        << zh_centroid[i] << std::setw(15) << pt_centroid[i] << std::setw(15) << phi_centroid[i] << std::setw(15) << f1 << std::setw(15)
        << f3 << std::endl;
  }
  out.close();

  return 0;
}
