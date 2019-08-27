#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// constants
const double me = 0.000511;
const double mN = 0.93892;
const double mU=0.931;
const double GeVfm=0.1973;
const double alpha=0.0072973525664;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;

// nuclear masses
const double m_1H = mN;
const double m_2H = 2.01410178 * mU - me;
const double m_3H = 3.01604928199 * mU - me;
const double m_3He = 3.0160293 * mU - 2*me;
const double m_4He = 4.00260325415 * mU - 2*me;
const double m_6Li = 6.015122795 * mU - 3*me;
const double m_8Be = 8.00530510 * mU - 4*me;
const double m_10Be = 10.013534 * mU - 4*me;
const double m_10B = 10.0129370 * mU - 5*me;
const double m_11B = 11.0093054 * mU - 5*me;
const double m_10C = 10.016853 * mU - 6*me;
const double m_12C = 12.*mU - 6*me;
const double m_14N = 14.0030740048*mU - 7*me;
const double m_16O = 15.99491461956*mU - 8*me;

// nucleon codes
const int pCode = 2212;
const int nCode = 2112;

#endif
