/*****************************************/
/*  GetBinning.cxx                       */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program prints the binning of (Q2, Nu, Zh, Pt2, PhiPQ) into a CSV file
// from arrays defined in include/Headers.hxx

// October 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main() {

  gProgram = "GetBinning";

  /*** WRITE OUTPUT FILE ***/

  TString OutFileName = "binning.csv";
  std::ofstream OutFile(OutFileName, std::ios::out);

  Int_t NbinsQ2 = (Int_t)(sizeof(kEdgesQ2) / sizeof(kEdgesQ2[0])) - 1;
  Int_t NbinsNu = (Int_t)(sizeof(kEdgesNu) / sizeof(kEdgesNu[0])) - 1;
  Int_t NbinsZ = (Int_t)(sizeof(kEdgesZ) / sizeof(kEdgesZ[0])) - 1;
  Int_t NbinsPt2 = (Int_t)(sizeof(kEdgesPt2) / sizeof(kEdgesPt2[0])) - 1;
  Int_t NbinsPhiPQ = (Int_t)(sizeof(kEdgesPhiPQ) / sizeof(kEdgesPhiPQ[0])) - 1;

  std::cout << "Binning : (Q2, Nu, Z, Pt2, PhiPQ) = " << NbinsQ2 << " x " << NbinsNu << " x " << NbinsZ << " x " << NbinsPt2 << " x "
            << NbinsPhiPQ << std::endl;

  for (Int_t q = 0; q < NbinsQ2; q++) {
    for (Int_t n = 0; n < NbinsNu; n++) {
      for (Int_t z = 0; z < NbinsZ; z++) {
        for (Int_t t = 0; t < NbinsPt2; t++) {
          for (Int_t p = 0; p < NbinsPhiPQ; p++) {
            OutFile << kEdgesQ2[q] << "," << kEdgesQ2[q + 1] << "," << kEdgesNu[n] << "," << kEdgesNu[n + 1] << "," << kEdgesZ[z] << ","
                    << kEdgesZ[z + 1] << "," << kEdgesPt2[t] << "," << kEdgesPt2[t + 1] << "," << kEdgesPhiPQ[p] << ","
                    << kEdgesPhiPQ[p + 1] << std::endl;
          }
        }
      }
    }
  }

  OutFile.close();
  std::cout << "The following file has been created: " << OutFileName << std::endl;

  return 0;
}
