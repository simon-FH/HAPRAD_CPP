/*****************************************/
/*  GetBinning.cxx                       */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program prints the binning of (Q2, Nu, Zh, Pt2, PhiPQ) into a CSV file
// from arrays defined in include/Headers.hxx

// June 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main() {

  gProgram = "GetBinning";

  /*** WRITE OUTPUT FILE ***/

  TString OutFileName = "binning.csv";
  std::ofstream OutFile(OutFileName, std::ios::out);

  Int_t NbinsQ2 = (Int_t)(sizeof(kEdgesQ2) / sizeof(kEdgesQ2[0]));
  Int_t NbinsNu = (Int_t)(sizeof(kEdgesNu) / sizeof(kEdgesNu[0]));
  Int_t NbinsZ = (Int_t)(sizeof(kEdgesZ) / sizeof(kEdgesZ[0]));
  Int_t NbinsPt2 = (Int_t)(sizeof(kEdgesPt2) / sizeof(kEdgesPt2[0]));
  Int_t NbinsPhiPQ = (Int_t)(sizeof(kEdgesPhiPQ) / sizeof(kEdgesPhiPQ[0]));

  std::cout << "Binning : (Q2, Nu, Z, Pt2, PhiPQ) = " << NbinsQ2 - 1 << " x " << NbinsNu - 1 << " x " << NbinsZ - 1 << " x " << NbinsPt2 - 1
            << " x " << NbinsPhiPQ - 1 << std::endl;

  for (Int_t q = 0; q < NbinsQ2 - 1; q++) {
    for (Int_t n = 0; n < NbinsNu - 1; n++) {
      for (Int_t z = 0; z < NbinsZ - 1; z++) {
        for (Int_t t = 0; t < NbinsPt2 - 1; t++) {
          for (Int_t p = 0; p < NbinsPhiPQ - 1; p++) {
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
