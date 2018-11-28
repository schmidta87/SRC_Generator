#ifndef __CROSS_SECTIONS_H__
#define __CROSS_SECTIONS_H__

#include "TVector3.h"

class Cross_Sections
{
 public:
  Cross_Sections();
  double sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  
 private:
  double Gdipole(double QSq);
  double sq(double x);
  
};

#endif
