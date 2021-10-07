#include "TSemiInclusiveModel.h"
#include "TMath.h"
#include "Partons.h"
#include "ConfigFile.h"
#include "haprad_constants.h"
#include "square_power.h"
#include <iostream>

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "THn.h"

namespace HapradUtils {

void SemiInclusiveModel(Double_t  q2,  Double_t  X,
                        Double_t  Y,   Double_t  Z,
                        Double_t  pt2, Double_t  mx2,
                        Double_t  pl,
                        Double_t& A, Double_t& Ac,
                        Double_t& Acc)
{
    using namespace TMath;

    Int_t bins[4] = {6,5,10,5};
    Double_t xmins[4] = {1.,0.1,0.,0.};
    Double_t xmaxs[4] = {4.,0.55,1.,1.};
    static THnD H12Hist("h12_hist", "H1+H2 histo", 4, bins, xmins, xmaxs);
    static THnD H3Hist("h3_hist", "H3 histo", 4, bins, xmins, xmaxs);
    static THnD H4Hist("h4_hist", "H4 histo", 4, bins, xmins, xmaxs);
    Double_t values[4];

    static bool initialized;

    gSystem->Load("libTree");

    if (!initialized){
	    initialized = true;
	    TFile file("newphihist.root");
	    TNtuple *tuple = (TNtuple*)file.Get("AAcAcc_data");
	    Float_t rXb, rQ2, rZh, rPt, rA, rAc, rAcc;
	    tuple->SetBranchAddress("Q2",&rQ2);
	    tuple->SetBranchAddress("Xb",&rXb);
	    tuple->SetBranchAddress("Zh",&rZh);
	    tuple->SetBranchAddress("Pt",&rPt);
	    tuple->SetBranchAddress("A",&rA);
	    tuple->SetBranchAddress("Ac",&rAc);
	    tuple->SetBranchAddress("Acc",&rAcc);
	    for (int i = 0; i < tuple->GetEntries(); i++){
		    tuple->GetEntry(i);
		    values[0] = rQ2; values[1] = rXb;
		    values[2] = rZh; values[3] = rPt;
		    H12Hist.Fill(values,rA);
		    H3Hist.Fill(values,rAc);
		    H4Hist.Fill(values,rAcc);

	    }
    }

    /*Check the kinematics*/
    if (X < 0 || X > 1) return;
    if (Z < 0 || Z > 1) return;

    //Double_t r  = Sqrt(1. + SQ(2 * kMassProton * X) / q2);

    if (mx2 < SQ(kMassProton + kMassPion)) return;

    Double_t bin[4];
    bin[0] = q2; bin[1] = X; bin[2] = Z; bin[3] = sqrt(pt2);

    for(Int_t i=0; i<4; i++){
      if(bin[i] > xmaxs[i]) bin[i] = xmaxs[i]-0.001;
      if(bin[i] < xmins[i]) bin[i] = xmins[i]+0.001;
    }
    
//    std::cout << "Vals: " <<  q2 << " " << X << " " << X << " " << H3Hist.GetBinContent(H4Hist.GetBin(bin)) << std::endl;
//    std::cout << q2 << " " << X << " " << X << " " << H4Hist.GetBinContent(H4Hist.GetBin(bin))  << std::endl;
    
    A = H12Hist.GetBinContent(H12Hist.GetBin(bin));   
    Ac = H3Hist.GetBinContent(H3Hist.GetBin(bin));
    Acc = H4Hist.GetBinContent(H4Hist.GetBin(bin));

    /*CHECK LATER
    // Check that cos(phi) and cos(2phi) are less than 1
    Double_t Ebeam, rt, rtz, cterm, m_cos_phi, m_cos_2phi;
    Ebeam   = q2 / (2 * kMassProton * X * Y);
    rt      = 1 - Y - kMassProton * X * Y / (2 * Ebeam);
    rtz     = Sqrt(rt / (1 + 2 * kMassProton * X / (Y * Ebeam)));
    cterm   = X * SQ(Y) * H1 + rt * H2;

    m_cos_phi = Sqrt(pt2 / q2) * (2 - Y) * rtz * H3m / (2 * cterm);

    if (Abs(4. * m_cos_phi) > 0.9) {
        H3m = 0.9 * Sign(1., m_cos_phi) * (2 * cterm) /
                (Sqrt(pt2 / q2) * (2 - Y) * rtz) / 4.;
    }

    m_cos_2phi = pt2 / q2 * SQ(rtz) * H4m / (2 * cterm);

    if (Abs(4 * m_cos_2phi) > 0.9)
        H4m = 0.9 * Sign(1., m_cos_2phi) * (2 * cterm) / (pt2 / q2 * SQ(rtz)) / 4.;
    */

    return;
}

} // End namespace HapradUtils
