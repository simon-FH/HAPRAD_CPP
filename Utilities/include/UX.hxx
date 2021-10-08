#ifndef UX_HXX
#define UX_HXX

#ifndef HEADERS_HXX
#include "Headers.hxx"
#endif

/*** Global variables ***/

TString gDataDir = getenv("DATA_DIR");
TString gSimDir = getenv("SIM_DIR");

TString gTargetOption;
Int_t gPID;

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
    }
  }
}

void parseCommandLine(int argc, char* argv[]) {
  Int_t c;
  if (argc == 1) {
    std::cerr << "Empty command line. Execute ./" << gProgram << " -h to print help." << std::endl;
    exit(0);
  }
  while ((c = getopt(argc, argv, "ht:p:")) != -1) switch (c) {
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
    }
  }
  std::cout << std::endl;
}

#endif
