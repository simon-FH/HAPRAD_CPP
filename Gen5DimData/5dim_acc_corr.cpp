#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1.h"

TString dataLoc;
TString fDataExt;
TString fSimuExt;
TString elecExt;
TString pionExt;

Int_t nSimuFiles;

Double_t Q2_MIN;
Double_t Q2_MAX;
Double_t XB_MIN;
Double_t XB_MAX;
Double_t NU_MIN;
Double_t NU_MAX;
Double_t ZH_MIN;
Double_t ZH_MAX;
Double_t PT_MIN;
Double_t PT_MAX;
Double_t PHI_MIN;
Double_t PHI_MAX;

Int_t N_Q2;
Int_t N_XB;
Int_t N_NU;
Int_t N_ZH;
Int_t N_PT;
Int_t N_PHI;

Int_t XF_POS = 0;

Double_t delta_Q2;
Double_t delta_XB;
Double_t delta_NU;
Double_t delta_ZH;
Double_t delta_PT;
Double_t delta_PHI;

Int_t N_METAL = 6;

TFile *plots;
TNtuple *fit_data;

void run_file(TString Metal, Double_t q2_min, Double_t q2_max, Double_t xb_min, Double_t xb_max, Double_t zh_min, Double_t zh_max, Double_t pt_min, Double_t pt_max);

int main(int argc, char **argv)
{
	TString Metal;
	TString rootFName;
	Int_t dataSL;

	if(argc < 25 || argc > 28){
		std::cout << "The number of arguments is incorrect" << std::endl;
		return 0;
	}

	Q2_MIN = (Double_t) std::stod(argv[1]);
	Q2_MAX = (Double_t) std::stod(argv[2]);
	N_Q2 = (Int_t) std::stoi(argv[3]);
	XB_MIN = (Double_t) std::stod(argv[4]);
	XB_MAX = (Double_t) std::stod(argv[5]);
	N_XB = (Int_t) std::stoi(argv[6]);
	NU_MIN = (Double_t) std::stod(argv[7]);
	NU_MAX = (Double_t) std::stod(argv[8]);
	N_NU = (Int_t) std::stoi(argv[9]);
	ZH_MIN = (Double_t) std::stod(argv[10]);
	ZH_MAX = (Double_t) std::stod(argv[11]);
	N_ZH = (Int_t) std::stoi(argv[12]); 
	PT_MIN = (Double_t) std::stod(argv[13]);
	PT_MAX = (Double_t) std::stod(argv[14]);
	N_PT = (Int_t) std::stoi(argv[15]);
	PHI_MIN = (Double_t) std::stod(argv[16]);
	PHI_MAX = (Double_t) std::stod(argv[17]);
	N_PHI = (Int_t) std::stoi(argv[18]);

	dataLoc = (TString) argv[19];
	fDataExt = (TString) argv[20];
	fSimuExt = (TString) argv[21];
	nSimuFiles = (Int_t) std::stoi(argv[22]);
	elecExt = (TString) argv[23];
	pionExt = (TString) argv[24];

	if(argc >= 26) XF_POS = (Int_t) std::stoi(argv[25]);
	if(argc >= 27) Metal = argv[26];
	
	dataSL = dataLoc.Length();
	if(dataLoc(dataSL-1, 1) != "/") dataLoc = dataLoc + "/";
	
	std::cout << "The following settings are being used" << std::endl;
	std::cout << Q2_MIN << " < Q2 < " << Q2_MAX << ", N_Q2 = " << N_Q2 << std::endl;
	std::cout << XB_MIN << " < Xb < " << XB_MAX << ", N_XB = " << N_XB << std::endl;
	std::cout << ZH_MIN << " < Zh < " << ZH_MAX << ", N_ZH = " << N_ZH << std::endl;	
	std::cout << PT_MIN << " < Pt < " << PT_MAX << ", N_PT = " << N_PT << std::endl;	
	std::cout << PHI_MIN << " < Phi < " << PHI_MAX << ", N_PHI = " << N_PHI << std::endl;	
	std::cout << std::endl;
	std::cout << "Data Directory = " << dataLoc << std::endl;
	std::cout << "Data Files Extension = " << fDataExt << std::endl;
	std::cout << "Simul Files Extension = " << fSimuExt << std::endl;
	std::cout << "Number of Simu Files = " << nSimuFiles << std::endl;
	std::cout << "Electron Files Extension = " << elecExt << std::endl;
	std::cout << "Pion Files Extension = " << pionExt << std::endl;
	std::cout << std::endl;
	if(XF_POS == 1) std::cout << "Xf Cut = Xf > 0" << std::endl;
	else if(XF_POS == -1) std::cout << "Xf Cut = Xf < 0" << std::endl;
	else std::cout << "No Xf Cut" << std::endl; 
	if(argc == 27) std::cout << "Running only for Metal " << Metal << std::endl;

	delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
	delta_XB = (XB_MAX-XB_MIN)/N_XB;
	delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
	delta_PT = (PT_MAX-PT_MIN)/N_PT;
	delta_PHI = (PHI_MAX-PHI_MIN)/N_PHI;
	delta_NU = (NU_MAX-NU_MIN)/N_NU;
	
	if(Metal != "") N_METAL = 1;

	for(Int_t met = 0; met < N_METAL; met++){
		if(N_METAL == 6){
			if(met == 0) Metal = "C";
			else if(met == 1) Metal = "Fe";
			else if(met == 2) Metal = "Pb";
			else if(met == 3) Metal = "D_C";
			else if(met == 4) Metal = "D_Fe";
			else if(met == 5) Metal = "D_Pb";
		}
		rootFName = Metal + "_5_dim_dist.root";
		plots = new TFile(rootFName, "RECREATE");
		fit_data = new TNtuple("fit_data", "DATA FOR 5 DIM FIT", "Q2:Xb:Zh:Pt:Phi:Val:Err");

		for(Int_t i = 0; i < N_Q2; i++){
		    for(Int_t j = 0; j < N_XB; j++){
		        for(Int_t k = 0; k < N_ZH; k++){
		            for(Int_t l = 0; l < N_PT; l++){
		                run_file(Metal, Q2_MIN+i*delta_Q2, Q2_MIN+(i+1)*delta_Q2,
		                            XB_MIN+j*delta_XB, XB_MIN+(j+1)*delta_XB,
		                            ZH_MIN+k*delta_ZH, ZH_MIN+(k+1)*delta_ZH,
		                            PT_MIN+l*delta_PT, PT_MIN+(l+1)*delta_PT);
		                plots->cd();
		            }
		        }
		    }
		}
		plots->cd();
		fit_data->Write();
		fit_data->Delete();
		plots->Close();

		delete plots;
	}

    return 0;
}

void run_file(TString Metal, Double_t q2_min, Double_t q2_max, Double_t xb_min, Double_t xb_max, Double_t zh_min, Double_t zh_max, Double_t pt_min, Double_t pt_max)
{	
  TString MetalD;
  TCut Q2_cut = Form("Q2>%f && Q2<%f", q2_min, q2_max);
  TCut Xb_cut = Form("Xb>%f && Xb<%f", xb_min, xb_max);
  TCut Zh_cut = Form("Zh>%f && Zh<%f", zh_min, zh_max);
  TCut Pt_cut = Form("Pt>%f && Pt<%f", pt_min, pt_max);
  TCut Xf_cut;
  if(XF_POS == 0) Xf_cut = Form("Xf<0 || Xf>0");
  if(XF_POS == 1) Xf_cut = Form("Xf>0");
  if(XF_POS == -1) Xf_cut = Form("Xf<0");
  TCut Target_cut;
  if(Metal(0,1) == "D")
    Target_cut = Form("TargType==%d", 1);
  else
    Target_cut = Form("TargType==%d", 2);
  
  TCut cuts = Q2_cut && Xb_cut && Zh_cut && Pt_cut && Xf_cut && Target_cut;
  TCut cuts_simul = Q2_cut && Xb_cut && Zh_cut&& Pt_cut && Xf_cut;
  
  TChain *fntuple, *faccept, *fthrown;
  
  if(Metal(0,1) != "D"){
    //		fntuple = new TChain("data_pion");
    fntuple = new TChain("ntuple_data");
    
    fntuple->Add(dataLoc + Metal + fDataExt + pionExt);
    //	faccept = new TChain("accept_pion");
    faccept = new TChain("ntuple_accept");
    for(Int_t q = 0; q < nSimuFiles; q++)
      faccept->Add(dataLoc + Metal + std::to_string(static_cast<long long>(q+1)) + fSimuExt + pionExt);
    //		fthrown = new TChain("thrown_pion");
    fthrown = new TChain("ntuple_thrown");
    for(Int_t q = 0; q < nSimuFiles; q++)
      fthrown->Add(dataLoc + Metal + std::to_string(static_cast<long long>(q+1)) + fSimuExt + pionExt);		
  }
  else{
    MetalD = Metal(1,2);
    //		fntuple = new TChain("data_pion");
    fntuple = new TChain("ntuple_data");
    fntuple->Add(dataLoc + MetalD + fDataExt + pionExt);
    //		faccept = new TChain("accept_pion");
    faccept = new TChain("ntuple_accept");
    for(Int_t q = 0; q < nSimuFiles; q++)
      faccept->Add(dataLoc + Metal(0,1) + std::to_string(static_cast<long long>(q+1)) + fSimuExt + pionExt);
    //		fthrown = new TChain("thrown_pion");
    fthrown = new TChain("ntuple_thrown");
    for(Int_t q = 0; q < nSimuFiles; q++)
      fthrown->Add(dataLoc + Metal(0,1) + std::to_string(static_cast<long long>(q+1)) + fSimuExt + pionExt);
  }	

  fthrown->SetAlias("Xb","Q2/2./Nu/0.938272");
  faccept->SetAlias("Xb","Q2/2./Nu/0.938272");

  fntuple->Draw(">>list", cuts, "goff");
  fntuple->SetEventList((TEventList*) gDirectory->Get("list"));   
  

  faccept->Draw(">>list_acc",cuts_simul,"goff");
  faccept->SetEventList((TEventList*) gDirectory->Get("list_acc"));


  fthrown->Draw(">>list_thr", cuts_simul, "goff");
  fthrown->SetEventList((TEventList*)gDirectory->Get("list_thr"));

  fntuple->Draw((const char*)Form("PhiPQ>>htmp_data(%d,%f,%f)", N_PHI, PHI_MIN, PHI_MAX), "", "goff");
  TH1F *htmp_data = (TH1F*)gDirectory->GetList()->FindObject("htmp_data");
  htmp_data->Sumw2();

  faccept->Draw((const char*)Form("PhiPQ>>htmp_acc(%d,%f,%f)", N_PHI, PHI_MIN, PHI_MAX), "", "goff");
  TH1F *htmp_acc = (TH1F*)gDirectory->GetList()->FindObject("htmp_acc");
  htmp_acc->Sumw2();

  fthrown->Draw((const char*)Form("PhiPQ>>htmp_thr(%d,%f,%f)", N_PHI, PHI_MIN, PHI_MAX), "", "goff");
  TH1F *htmp_thr = (TH1F*)gDirectory->GetList()->FindObject("htmp_thr");
  htmp_thr->Sumw2();

  TH1F *htmp_acc_ratio = new TH1F("htmp_acc_ratio", "", N_PHI, PHI_MIN, PHI_MAX);
  TH1F *htmp_data_corr = new TH1F("htmp_data_corr", "", N_PHI, PHI_MIN, PHI_MAX);
  htmp_acc_ratio->Divide(htmp_acc,htmp_thr,1,1,"B");
  htmp_data_corr->Divide(htmp_data,htmp_acc_ratio,1,1);
  
  for(Int_t ii = 1; ii <= N_PHI; ii++){
    Float_t PhiVal = PHI_MIN + (ii-0.5)*delta_PHI;
    Float_t Val = htmp_data_corr->GetBinContent(ii);
    Float_t Err = htmp_data_corr->GetBinError(ii);
    fit_data->Fill((q2_min+q2_max)/2., (xb_min+xb_max)/2., (zh_min+zh_max)/2., (pt_min+pt_max)/2., PhiVal, Val, Err);
  }
  
  delete htmp_acc_ratio;
  delete htmp_data_corr;
  delete fntuple;
  delete faccept;
  delete fthrown;    
}
