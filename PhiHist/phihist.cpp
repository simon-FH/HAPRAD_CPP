#include <fstream>
#include <iostream>
#include <string>
#include "TString.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

Double_t PHI_MIN;
Double_t PHI_MAX;
Int_t N_PHI;
Int_t N_METAL;
const Int_t MIN_BIN = 4;

int main(int argc, char **argv)
{
	TFile *f;
	TFile *newf;
	TNtuple *ntuple;
	TNtuple *newntuple;
	TH1F *h1;
	TF1 *func;
	TString Metal;
	TString dataLoc;
	Int_t dataSL;
	
	if(argc < 5 || argc > 6){
		std::cout << "The number of arguments is incorrect" << std::endl;
		return 0;
	}
	
	dataLoc = argv[1];
	PHI_MIN = (Double_t) std::stod(argv[2]);
	PHI_MAX = (Double_t) std::stod(argv[3]);
	N_PHI = (Int_t) std::stoi(argv[4]);
	if(argc >= 6) Metal = argv[5];
		
	dataSL = dataLoc.Length();
	if(dataLoc(dataSL-1, 1) != "/"){
		dataLoc = dataLoc + "/";
	}	
	
	std::cout << "The following settings are being used" << std::endl;
	std::cout << "Data Directory = " << dataLoc << std::endl;
	std::cout << PHI_MIN << " < PHI < " << PHI_MAX << ", N_PHI = " << N_PHI << std::endl;
	if(argc == 6) std::cout << "Running only for Metal "<< Metal << std::endl;
	else std::cout << "Running all Metals" << std::endl;	
	
	Float_t Q2, Xb, Zh, Pt, Phi, Val, Err;
	Float_t A, Ac, Acc;
	Float_t AErr, AcErr, AccErr;
	Float_t ChiSQ;

	Int_t nentries, empty;

	if(Metal == "") N_METAL = 6;
	else N_METAL = 1;
	
	for(Int_t met = 0; met < N_METAL; met++){
		if(Metal == ""){
			if(met == 0) Metal = "C";
			else if(met == 1) Metal = "Fe";
			else if(met == 2) Metal = "Pb";
			else if(met == 3) Metal = "D_C";
			else if(met == 4) Metal = "D_Fe";
			else if(met == 5) Metal = "D_Pb";
		}
		if(dataLoc == "")
			f = new TFile("../Gen5DimData/" + Metal + "_5_dim_dist.root", "READ");
		else
			f = new TFile(dataLoc + Metal + "_5_dim_dist.root", "READ");
		ntuple = (TNtuple*) f->Get("fit_data");

		nentries = ntuple->GetEntries();
		ntuple->SetBranchAddress("Xb", &Xb);
		ntuple->SetBranchAddress("Q2", &Q2);
		ntuple->SetBranchAddress("Xb", &Xb);
		ntuple->SetBranchAddress("Zh", &Zh);
		ntuple->SetBranchAddress("Pt", &Pt);
		ntuple->SetBranchAddress("Phi", &Phi);
		ntuple->SetBranchAddress("Val", &Val);
		ntuple->SetBranchAddress("Err", &Err);

		newntuple = new TNtuple("AAcAcc_data", "AAcAcc_data", "Q2:Xb:Zh:Pt:A:AErr:Ac:AcErr:Acc:AccErr:ChiSQ");
		func = new TF1("fit", "[0]+TMath::Cos(x*TMath::Pi()/180)*[1]+TMath::Cos(2*x*TMath::Pi()/180)*[2]");

		newf = new TFile(Metal + "newphihist.root", "RECREATE");
		newf->cd();

		for(Int_t i = 0; i < nentries; i = i + N_PHI){
			ntuple->GetEntry(i);
			h1 = new TH1F((const char*) Form("PhiDist Q2=%.3f Xb=%.3f Zh=%.3f Pt=%.3f", Q2, Xb, Zh, Pt), 
							(const char*) Form("PhiDist Q2=%.3f Xb=%.3f Zh=%.3f Pt=%.3f", Q2, Xb, Zh, Pt), N_PHI, PHI_MIN, PHI_MAX);
			empty = 0;
			for(Int_t j = 1; j <= N_PHI; j++){
				ntuple->GetEntry(i+j-1);
				h1->SetBinContent(j, Val);
				if(h1->GetBinContent(j) == 0)
					empty++;
				h1->SetBinError(j, Err*1.04);
			}
			if(empty <= (N_PHI - MIN_BIN)){
				h1->Fit(func, "q");
				h1->Write();
				if(func->GetNDF() != 0){
					ChiSQ = func->GetChisquare();
					A = func->GetParameter(0);
					AErr = func->GetParError(0);
					Ac = func->GetParameter(1);
					AcErr = func->GetParError(1);
					Acc = func->GetParameter(2);
					AccErr = func->GetParError(2);
					newntuple->Fill(Q2, Xb, Zh, Pt, A, AErr, Ac, AcErr, Acc, AccErr, ChiSQ);
				}
			}
			h1->Delete();
		}
		
		Metal = "";
		
		delete ntuple;
		
		newf->cd();
		newntuple->Write();
		newntuple->Delete();
		newf->Close();
		delete newf;
		f->Close();
		delete f;
	}
	return 0;
}
