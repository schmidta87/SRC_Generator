#include "detectors.h"
#include "TRandom3.h"

TRandom3 myRand(0);

bool HRS_hallA(TVector3 &v, double p_central, double phi_central)
{

  double theta_central = 0.5*TMath::Pi();

  v.RotateY(0.5*TMath::Pi());
  
  double p = v.Mag();
  double theta = v.Theta();
  double phi = v.Phi();

  if (abs(p/p_central - 1) < 0.045)
    return false;

  if (abs(phi - phi_central) < 0.030)
    return false;
  
  if (abs(theta - theta_central) < 0.060)
    return false;
  
  p *= myRand.Gaus(1,1e-4/(2*sqrt(2*log(2))));
  phi += myRand.Gaus(0,0.5e-3/(2*sqrt(2*log(2))));
  theta += myRand.Gaus(0,1.0e-3/(2*sqrt(2*log(2))));

  v.RotateY(-0.5*TMath::Pi());
  
  return true;
}

bool BigBite(TVector3 &v, double p_central, double phi_central)
{

  double theta_central = 0.5*TMath::Pi();

  v.RotateY(0.5*TMath::Pi());
  
  double p = v.Mag();
  double theta = v.Theta();
  double phi = v.Phi();

  if (abs(p/p_central - 1) < 0.045)
    return false;

  if (abs(phi - phi_central) < 0.030)
    return false;
  
  if (abs(theta - theta_central) < 0.060)
    return false;
  
  p *= myRand.Gaus(1,1.6e-3);
  phi += myRand.Gaus(0,7e-3);
  theta += myRand.Gaus(0,13e-3);

  v.RotateY(-0.5*TMath::Pi());
  
  return true;
}
