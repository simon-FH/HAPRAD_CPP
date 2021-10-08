/*****************************************/
/*  GetRC.cxx                            */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// June 2021

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

  Double_t f1, f3;
  Double_t m = TMath::Power((kMassNeutron + kMassPion), 2);  // ...ok?
  Double_t a1, a2, a3;
  Int_t sysReturn;
  Int_t dataSL;

  TRadCor rc;

  a1 = 1;
  a2 = 0;
  a3 = a1 / a2;

  std::ofstream out;
  out.open("RCFactor_" + gTargetOption + ".txt");
  out << "Phi\tSigmaB\tSigmaOb\tTail1\tTaile2\tFact_noex\tFact_ex" << std::endl;
  for (Int_t i = 0; i < (Int_t)phi_centroid.size(); i++) {
    rc.CalculateRCFactor(5.015, xb_centroid[i], q2_centroid[i], zh_centroid[i], pt_centroid[i], phi_centroid[i], m, NAZ);
    f1 = rc.GetFactor1();
    f3 = rc.GetFactor3();
    if (TMath::IsNaN(f1) || f1 == a3) f1 = 0;
    if (TMath::IsNaN(f3) || f3 == a3) f3 = 0;
    out << i << "\t0\t0\t0\t0\t" << f1 << "\t" << f3 << std::endl;
  }
  out.close();

  return 0;
}
