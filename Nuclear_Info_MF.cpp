#include "Nuclear_Info_MF.h"

#include "constants.h"

#include <iostream>
#include <cstdlib>

#include <fstream>
#include <sstream>

using namespace std;

Nuclear_Info_MF::Nuclear_Info_MF(int thisA)
{
  A = thisA;
  
  fill_arrays();

  Estar = 0;
  sigmaE = 0;
  
  if (A==4)
    {
      mA=m_4He;
      mAmp=m_3H;
      mAmn=m_3He;
    }
  else
    {
      std::cerr << "You selected a nucleus with A=" << A << "\n"
	   << " which is not in the library. Aborting...\n";
      exit(-2);
    }
}

Nuclear_Info_MF::~Nuclear_Info_MF()
{
  
}

double Nuclear_Info_MF::get_mA()
{
  return mA;
}

double Nuclear_Info_MF::get_mAmp()
{
  return mAmp;
}

double Nuclear_Info_MF::get_mAmn()
{
  return mAmn;
}

double Nuclear_Info_MF::get_mAm1()
{
  return (mAmp+mAmn)/2;
}

double Nuclear_Info_MF::get_Pp(double k, double E)
{
  return get_P(Pp,k,E);
}

double Nuclear_Info_MF::get_Pn(double k, double E)
{
  return get_P(Pn,k,E);
}

double Nuclear_Info_MF::get_P(double PMF[200][1000], double k, double E)
{
  double bink = (k / GeVfm + 0.025) / 0.05;
  double binE = E / 0.001 + 0.5;

  if ((bink < 0.) or (binE < 0.) or (bink > 200.) or (binE > 1000.))
    return 0;
  if ((bink < 1.) and (binE < 1.))
    return bink * binE * PMF[0][0];

  int bk = bink;
  double xk = bink - bk;

  if (binE < 1.)
    return binE * (xk*PMF[bk][0] + (1.-xk)*PMF[bk-1][0]);
  
  int bE = binE;
  double xE = binE - bE;

  if (bink < 1.)
    return bink * (xE*PMF[0][bE] + (1.-xE)*PMF[0][bE-1]);

  return (xk*xE*PMF[bk][bE] + xk*(1.-xE)*PMF[bk][bE-1] + (1.-xk)*xE*PMF[bk-1][bE] + (1.-xk)*(1.-xE)*PMF[bk-1][bE-1]);
}

void Nuclear_Info_MF::fill_arrays()
{

  double temp_Pp[200][1000];
  
  double temp_Pn[200][1000];
  
  std::ifstream filep("../SRC_Generator/SpectralFunctions_MF/Pp_4He.csv");
  
  for (int row = 0; row < 200; ++row)
    {
      std::string line;
      std::getline(filep,line);
      if (!filep.good())
	break;
      
      std::stringstream iss(line);

      for (int col = 0; col < 1000; ++col)
	{
	  std::string val;
	  std::getline(iss,val,',');
	  if (!iss.good())
	    break;

	  std::stringstream converter(val);
	  converter >> temp_Pp[row][col];
	  
	  Pp[row][col] = temp_Pp[row][col] * 1000 / (GeVfm*GeVfm*GeVfm);
	}
    }

  std::ifstream filen("../SRC_Generator/SpectralFunctions_MF/Pn_4He.csv");

  for (int row = 0; row < 200; ++row)
    {
      std::string line;
      std::getline(filen,line);
      if (!filen.good())
	break;
      
      std::stringstream iss(line);

      for (int col = 0; col < 1000; ++col)
	{
	  std::string val;
	  std::getline(iss,val,',');
	  if (!iss.good())
	    break;

	  std::stringstream converter(val);
	  converter >> temp_Pn[row][col];
	  
	  Pn[row][col] = temp_Pn[row][col] * 1000 / (GeVfm*GeVfm*GeVfm);
	}
      
    }
  
}
