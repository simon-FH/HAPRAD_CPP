#ifndef UX_HXX
#define UX_HXX

// getopt(), optarg. Older libstdc++ headers pulled this in transitively;
// current ones do not, so include it explicitly.
#include <unistd.h>

#ifndef HEADERS_HXX
#include "Headers.hxx"
#endif

/*** Global variables ***/

TString gDataDir = getenv("DATA_DIR");
TString gSimDir = getenv("SIM_DIR");

TString gTargetOption;
Int_t gPID;

// Beam energy in GeV. Default is the CLAS/EG2 value; RG-E runs are 10.3894,
// 10.4057 or 10.5473 depending on the run, and the ntuples carry E_beam
// per event.
Double_t gBeamEnergy = 5.015;

// Semi-inclusive structure-function table (TSemiInclusiveModel format). GetRC
// needs one: without it the amplitudes, sigma_Born and every RC factor are 0.
TString gTablePath = "";

TString gProgram = "";

/*** Input-related functions ***/

void printUsage() {
  std::cout << gProgram << " program. Usage is:" << std::endl;
  std::cout << std::endl;
  std::cout << "./" << gProgram << " -h" << std::endl;
  std::cout << "    prints this message and exits program" << std::endl;
  std::cout << std::endl;
  if (gProgram != "GetBinning") {
    std::cout << "./" << gProgram << " -t[D_C, D_Fe, D_Pb, C, Fe, Pb]" << std::endl;
    std::cout << "    filters respective target" << std::endl;
    std::cout << std::endl;
    if (gProgram != "GetRC") {
      std::cout << "./" << gProgram << " -p[211, -211, 2212]" << std::endl;
      std::cout << "    filters respective particle" << std::endl;
      std::cout << std::endl;
    } else {
      std::cout << "./" << gProgram << " -e<beam energy in GeV>" << std::endl;
      std::cout << "    default is " << gBeamEnergy << " (CLAS/EG2); "
                << "RG-E is 10.3894 / 10.4057 / 10.5473 by run" << std::endl;
      std::cout << std::endl;
      std::cout << "./" << gProgram << " -m<table.root>" << std::endl;
      std::cout << "    structure-function table (see TSemiInclusiveModel.h); required" << std::endl;
      std::cout << std::endl;
    }
  }
}

void parseCommandLine(int argc, char* argv[]) {
  Int_t c;
  if (argc == 1) {
    std::cerr << "Empty command line. Execute ./" << gProgram << " -h to print help." << std::endl;
    exit(0);
  }
  while ((c = getopt(argc, argv, "ht:p:e:m:")) != -1) switch (c) {
      case 'h':
        printUsage();
        exit(0);
        break;
      case 't':
        gTargetOption = optarg;
        break;
      case 'p':
        gPID = atoi(optarg);
        break;
      case 'e':
        gBeamEnergy = atof(optarg);
        break;
      case 'm':
        gTablePath = optarg;
        break;
      default:
        std::cerr << "Unrecognized argument. Execute ./" << gProgram << " -h to print help." << std::endl;
        exit(0);
        break;
    }
}

void printOptions() {
  std::cout << "Executing " << gProgram << " program. ";
  if (gProgram != "GetBinning") {
    std::cout << "The chosen parameters are: " << std::endl;
    std::cout << "  gTargetOption   = " << gTargetOption << std::endl;
    if (gProgram != "GetRC") {
      std::cout << "  gPID            = " << gPID << std::endl;
    } else {
      std::cout << "  gBeamEnergy     = " << gBeamEnergy << " GeV" << std::endl;
      std::cout << "  gTablePath      = " << gTablePath << std::endl;
    }
  }
  std::cout << std::endl;
}

#endif
