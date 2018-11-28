#include "Cross_Sections.h"

#include "constants.h"

#include "TVector3.h"

#include <cstdlib>

Cross_Sections::Cross_Sections()
{
}

double Cross_Sections::sq(double x){ return x*x; };

double Cross_Sections::sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  TVector3 q = TVector3(0.,0.,Ebeam) - k;
  TVector3 pM = p-q;

  double QSq = q.Mag2() - sq(Ebeam - k.Mag());
  double E = sqrt(p.Mag2() + sq(mN));
  double Ebar = sqrt(pM.Mag2() + sq(mN));
  double omegabar = E-Ebar;
  double QSqbar = q.Mag2() - sq(omegabar);

  // Calculate form factors
  double GE = (isProton)? Gdipole(QSq) : 1.91 * QSq * Gdipole(QSq) / (4.*sq(mN) + 5.6 * QSq);
  double GM = (isProton)? 2.79*Gdipole(QSq) : -1.91*Gdipole(QSq);
  double F1 = 0.5 * (GE + QSq*GM/(4.*sq(mN)));
  double kF2 = (GM - GE)/(1. + QSq/(4.*sq(mN)));

  double wC = (sq(E+Ebar)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2)) - q.Mag2()*sq(F1 + kF2))/(4.*E*Ebar);
  double wT = QSqbar*sq(F1 + kF2)/(2.*Ebar*E);
  double wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
  double wI = -p.Mag()*sin(p.Angle(q))*(Ebar + E)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);

  double sigmaMott = cmSqGeVSq * 4. * sq(alpha) * k.Mag2() * sq(cos(k.Theta()/2.)) / sq(QSq);

  double phi = q.Cross(k).Angle( q.Cross(p) );
  return sigmaMott * ( sq(QSq)/q.Mag2() * wC +
                       (QSq/(2.*q.Mag2()) + sq(tan(k.Theta()/2.))) * wT +
                       QSq/q.Mag2() * sqrt(QSq/q.Mag2() + sq(tan(k.Theta()/2.))) * wI * cos(phi) +
                       (QSq/q.Mag2() * sq(cos(phi)) + sq(tan(k.Theta()/2.))) * wS
                       );
}

double Cross_Sections::Gdipole(double QSq){ return 1. / sq(1 + QSq/0.71); };
