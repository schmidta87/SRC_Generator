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
  double get_phiSq(double *phiPtr, double k_rel);
  void fill_arrays();
  
};

#endif
