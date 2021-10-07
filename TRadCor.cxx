#include "TRadCor.h"
#include "THapradConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "TDelta.h"
#include "TBorn.h"
#include "TQQTPhi.h"
#include "TRV2TR.h"
#include "THapradException.h"
#include "haprad_constants.h"
#include "square_power.h"
#include "Math/GSLIntegrator.h"
#include "Math/GSLMCIntegrator.h"
#include "ConfigFile.h"
#include <iostream>
#include <iomanip>


TRadCor::TRadCor()
: sigma_born(0.), sig_obs(0.)
{
    // Default constructor

    fConfig = new THapradConfig();
}



TRadCor::~TRadCor()
{
    // Default destructor

    delete fConfig;
}



void TRadCor::RegisteredLepton(Int_t type)
{
    // Set the registere lepton: 1 -- electron; 2 -- muon.
    //
    // Default is 1.

    fConfig->SetLepton(type);
}



void TRadCor::IntegratePhiRad(Int_t type)
{
    // Set whether to integrate over phi_{rad} (1) or to approximate (0).
    //
    // Default is 0;

    fConfig->SetIntegrationPhiRad(type);
}



void TRadCor::IntegratePhiHad(Int_t type)
{
    // Set whether to integrate over phi_{had} (1) or not (0).
    //
    // Default is 0;

    fConfig->SetIntegrationPhiHad(type);
}



void TRadCor::SetPolarization(Int_t type)
{
    // Set the type of the polarization:
    //
    //     1 -- long; 2 -- tran; 0 -- unpol.
    //
    // Default is 0;

    fConfig->SetPolarization(type);
}



void TRadCor::CalculateRCFactor(Double_t E, Double_t x, Double_t Q2,
                                Double_t z, Double_t p_t, Double_t phi,
                                Double_t maxMx2, Double_t targProp)
{
    // Calculate, without returning it, the radiative correction factor for
    // the given parameters. You can use then one of the three methods:
    // GetFactor1, GetFactor2 and GetFactor3 to get the value you want.
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount
    // of missing mass.


    //Set up the propostion for protons to neutrons Z/A
    NAZ = targProp;

    fKin    = new TKinematicalVariables(x,-Q2,z,p_t,phi/kRadianDeg,E);
    fInv    = new TLorentzInvariants(fConfig,fKin);
    fHadKin = new THadronKinematics(fConfig,fKin,fInv);

    Double_t S_x = - fKin->Y() / fKin->X();//HH: check if it is correct
    //Double_t Mx2 = SQ(kMassProton) + S_x * (1 - fKin->Z()) + fKin->T();//HH: check if T is the correct one, if the evaluate is called correctly
    Double_t Nu2 = SQ(S_x/2/kMassProton);
    Double_t PlQ2NU = TMath::Sqrt((Nu2 + S_x*x)*(SQ(z)*Nu2 - SQ(kMassPion) - SQ(p_t)));
    Double_t t_temp = SQ(kMassPion) - S_x*x + 2*(PlQ2NU - Nu2*z);
    Double_t Mx2 = SQ(kMassProton) + S_x * (1 - fKin->Z()) + t_temp;//HH: check if T is the correct one, if the evaluate is called correctly


    fKinError = false;
    fParametersError = false;

    try {
        if (Mx2 > maxMx2) {
            Initialization();
            SPhiH();
        }
    } catch (TKinematicException& wrongKin) {
        fKinError = true;
        std::cerr << wrongKin.what() << std::endl;
    } catch (ConfigFile::file_not_found) {
        fParametersError = true;
        std::cout << std::endl
                  << "File with parameters not found. Please create a "
                  << "file named 'parameters'" << std::endl
                  << "in the same directory of your executable "
                  << "with this format:" << std::endl
                  << std::endl
                  << "     parameter1 = value" << std::endl
                  << "     parameter2 = value" << std::endl
                  << "     parameter2 = value" << std::endl
                  << "               ..." << std::endl
                  << std::endl
                  << "The parameters you should put in the file are:" << std::endl
                  << std::endl
                  << "    par0, par1, par2, par3, par4," << std::endl
                  << "    A1, B1, C1, D1, E1" << std::endl
                  << std::endl;
    }

    delete fKin;
    delete fInv;
    delete fHadKin;
}



Double_t TRadCor::GetRCFactor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                              Double_t p_t, Double_t phi, Double_t maxMx2, Double_t targProp)
{
    // Calculate and return the default radiative correction factor for the
    // given parameters.
    //
    // You can get the same value first calling the CalculateRCFactor(...)
    // method and then GetFactor2().
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    CalculateRCFactor(E,x,Q2,z,p_t,phi,maxMx2,targProp);

    if (fKinError || fParametersError)
        return 0;
    else
        return (sig_obs + tai[0] + tai[1]) / sigma_born;
}



Double_t TRadCor::GetFactor1(void)
{
    // Return the radiative correction factor without the exclusive radiative
    // tail contribution. You need to call CalculateRCFactor(...) first.

    if (fKinError || fParametersError) {
        return 0;
    } else {
        Double_t sigma_obs_f1 = sig_obs + tai[0];

        return sigma_obs_f1 / sigma_born;
    }
}



Double_t TRadCor::GetFactor2(void)
{
    // Return the radiative correction factor with the exclusive radiative
    // tail contribution. You need to call CalculateRCFactor(...) first.
    //
    // This is the default value, returned by GetRCFactor(...).

    if (fKinError || fParametersError) {
        return 0;
    } else {
        Double_t sigma_obs_f2 = sig_obs + tai[0] + tai[1];

        return sigma_obs_f2 / sigma_born;
    }
}



Double_t TRadCor::GetFactor3(void)
{
    // Return the radiative correction factor with the exclusive radiative
    // tail contribution divided by 2. You need to call CalculateRCFactor(...)
    // first.

    if (fKinError || fParametersError) {
        return 0;
    } else {
        Double_t sigma_obs_f3 = sig_obs + tai[0] + tai[1] * NAZ;

        return sigma_obs_f3 / sigma_born;
    }
}



void TRadCor::Initialization(void)
{
    fInv->Evaluate();
    fHadKin->Evaluate();

    // Calculate normalization factor
    N = kBarn * (kPi * SQ(kAlpha) * fKin->Y() * fInv->Sx() * kMassProton) /
                (2 * fInv->SqrtLq());
    if (fKin->T() >= 0) {
      if (fHadKin->Ph() > fHadKin->Pt()) {
	N = N * fInv->SqrtLq() / 2 / kMassProton / fHadKin->Pl();
      } else {
	N = 0;
      }
    }
}



void TRadCor::SPhiH(void)
{
    TDelta fDeltas(this);
    fDeltas.Evaluate();

    // Calculate Born cross-section
    std::cout << "********** ita: " << 1 << " *********" << std::endl;

    TBorn fBornin(this);
    sigma_born = N * fBornin.Evaluate();
    std::cout << "sib: " << sigma_born << std::endl;

    // Calculate radiative tail contribution
    if (sigma_born == 0.0) {
        tai[0] = 0.;
    } else {
      tai[0] = N * RadiativeTail();
    }

    std::cout << "tai[" << 0 << "]\t"  << tai[0] << std::endl;

    // Calculate exclusive radiative tail calculation
    std::cout << "********** ita: " << 2 << " *********" << std::endl;

    tai[1] = N * ExclusiveRadiativeTail();
    std::cout << "tai[" << 1 << "]\t"  << tai[1] << std::endl;

    // Calculate observed cross-section without tail contribution
    Double_t extai1 = TMath::Exp(fDeltas.Inf());
    sig_obs = sigma_born * extai1 * (1 + fDeltas.VR() + fDeltas.Vac());

    // Print results (as original FORTRAN code does)
    std::cout << std::endl;
    std::cout.setf(std::ios::fixed);
    std::cout << "    sib    " << std::setw(10)
              << std::setprecision(7) << sigma_born * 1E-3
              << "    sig    " << std::setw(10)
              << std::setprecision(7) << (sig_obs + tai[0] + tai[1]) * 1E-3
              << "    tail   " << std::setw(10)
              << std::setprecision(7) << tai[1] * 1E-3 << std::endl;
}



Double_t TRadCor::RadiativeTail(void)
{
    Double_t result = 0;

    TQQTPhi qphi(this);

    if (fConfig->IntegratePhiRad() == 1) {
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE);
        ig.SetFunction(qphi);
        ig.SetRelTolerance(fConfig->EpsPhiR());
        ig.SetAbsTolerance(1E-18);
        result = ig.Integral(0, TMath::TwoPi());
        result = kAlpha * result / SQ(kPi) / 4. / fInv->SqrtLq();
    } else if (fConfig->IntegratePhiRad() == 0) {
        result = kAlpha / kPi * qphi(0.) / 2 / fInv->SqrtLq();
    }

    return result;
}



Double_t TRadCor::ExclusiveRadiativeTail(void)
{
    Double_t result = 0;

    Double_t tau_1, tau_2;
    Double_t tau[6];
    Double_t phi[4];

    tau_1 = - fInv->Q2() / fInv->S();
    tau_2 =   fInv->Q2() / fInv->X();

    phi[0] = 0.;
    phi[1] = 0.01 * kPi;
    phi[2] = 2. * kPi - 0.01 * kPi;
    phi[3] = 2. * kPi;

    double tau_max = (fInv->Sx() + fInv->SqrtLq()) / (2. * SQ(kMassProton));
    double tau_min = - fInv->Q2() / SQ(kMassProton) / tau_max;
    tau[0] = tau_min;
    tau[1] = tau_1 - 0.15 * (tau_1 - tau_min);
    tau[2] = tau_1 + 0.15 * (tau_2 - tau_1);
    tau[3] = tau_2 - 0.15 * (tau_2 - tau_1);
    tau[4] = tau_2 + 0.15 * (tau_max - tau_2);
    tau[5] = tau_max;

    ROOT::Math::GSLMCIntegrator ig(ROOT::Math::IntegrationMultiDim::kMISER,
                                   1E-6,
                                   1E-3,
                                   20000);
    TRV2TR rv2tr(this);
    ig.SetFunction(rv2tr);

    Double_t rere = 0.;
    for (Int_t iph = 0; iph < 3; iph++) {
        for (Int_t ico = 0; ico < 5; ico++) {
            Double_t lower_lim[2], upper_lim[2];

            lower_lim[0] = tau[ico];
            upper_lim[0] = tau[ico+1];
            lower_lim[1] = phi[iph];
            upper_lim[1] = phi[iph+1];

            if (lower_lim[1] > upper_lim[1])
                std::cout << lower_lim[1] << " < " << upper_lim[1] << std::endl;
            if (lower_lim[0] > upper_lim[0])
                std::cout << lower_lim[0] << " < " << upper_lim[1] << std::endl;

            Double_t re = ig.Integral(lower_lim,upper_lim);

            std::cout.setf(std::ios::fixed);
            std::cout << " tai: "
                        << std::setw(4)  << ico + 1 << "  "
                        << std::setw(4)  << iph + 1 << "  "
                        << std::setw(10) << std::setprecision(4) << re << "  "
                        << std::endl;

            rere = rere + re;
        }
    }

    Double_t d = 64 * TMath::Power(kPi,5.) * fInv->SqrtLq() * kMassProton;
    result = - kAlpha * rere / d;

    return result;
}
