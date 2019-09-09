#include <cmath>
#include "detectors.h"
#include "TRandom3.h"
#include "helpers.h"
#include "constants.h"

TRandom3 myRand(0);

bool HRS_hallA(TVector3 &v, double p_central, double phi_central)
{

  double theta_central = 0.5*M_PI;

  double mom_FWHM = 1e-4;
  double h_FWHM = 0.5e-3;
  double v_FWHM = 1e-3;
  double mom_res = mom_FWHM/(2*sqrt(2*log(2)));
  double h_res = h_FWHM/(2*sqrt(2*log(2)));
  double v_res = v_FWHM/(2*sqrt(2*log(2)));
  
  double p = v.Mag();
  double theta = v.Theta();
  double phi = v.Phi();

  // Acceptance Region

  if (fabs(p/p_central - 1) > 4.5e-2)
    return false;

  if (fabs(phi - phi_central) > 30e-3)
    return false;
  
  if (fabs(theta - theta_central) > 60e-3)
    return false;

  // Resolution Smearing
  
  p *= (1. + mom_res*myRand.Gaus());
  phi += h_res*myRand.Gaus();
  theta += v_res*myRand.Gaus();

  v.SetPtThetaPhi(p,theta,phi);
  
  return true;
}

bool BigBite(TVector3 &v, double phi_central)
{

  double theta_central = 0.5*M_PI;

  double mom_res = 1.5e-2;
  double h_res = 7e-3;
  double v_res = 13e-3;

  double p = v.Mag();
  double theta = v.Theta();
  double phi = v.Phi();

  // Acceptance Region

  if (fabs(phi - phi_central) > 80e-3)
    return false;
  
  if (fabs(theta - theta_central) > 300e-3)
    return false;

  // Resolution Smearing
  
  p *= (1. + mom_res*myRand.Gaus());
  phi += h_res*myRand.Gaus();
  theta += v_res*myRand.Gaus();
  
  v.SetPtThetaPhi(p,theta,phi);
  
  return true;
}

bool HAND(TVector3 &v, double phi_central)
{

  double theta_central = 0.5*M_PI;
  double d = 20.01; //ns

  double p = v.Mag();
  double theta = v.Theta();
  double phi = v.Phi();

  // Acceptance Region
  
  if (fabs(phi - phi_central) > 80e-3)
    return false;
  
  if (fabs(theta - theta_central) > 300e-3)
    return false;

  // Resolution Smearing
  
  double tof = d*sqrt(1.+sq(mN/p));
  tof += 1.5*myRand.Gaus();
  p = mN/sqrt(1-sq(tof/d));
    
  phi += 15e-3*myRand.Gaus();
  theta += 15e-3*myRand.Gaus();
  
  v.SetPtThetaPhi(p,theta,phi);
  
  return true;
}
