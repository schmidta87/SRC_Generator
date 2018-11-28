#ifndef __NUCLEAR_INFO_H__
#define __NUCLEAR_INFO_H__

class Nuclear_Info
{
 public:
  Nuclear_Info(int thisA);
  ~Nuclear_Info();
  double get_pp(double k_rel);
  double get_pn(double k_rel);
  double get_pn0(double k_rel);
  double get_pn1(double k_rel);
  double get_mA();
  double get_mAm2();
  double get_sigmaCM();
  double get_sigmaE();
  double get_Estar();
  double get_Cpp0();
  double get_Cpn0();
  double get_Cpn1();

  void setCustomValues(double newSigma, double newEstar, double newCpp0, double newCpn0, double newCpn1);
  void set_sigmaCM(double newSigma);
  void set_Estar(double newEstar);
  void set_Cpp0(double newCpp0);
  void set_Cpn0(double newCpn0);
  void set_Cpn1(double newCpn1);
  void set_sigmaE(double newSigE);

  
 private:
  int A;
  double mA;
  double mAm2;
  double sigmaCM;
  double phiSq_pp0[100];
  double phiSq_pn0[100];
  double phiSq_pn1[100];
  double Cpp0;
  double Cpn0;
  double Cpn1;
  double Estar;
  double sigmaE;
  
  double get_phiSq(double *phiPtr, double k_rel);
  void fill_arrays();
  
};

#endif
