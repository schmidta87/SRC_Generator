#ifndef __NUCLEAR_INFO_H__
#define __NUCLEAR_INFO_H__

enum NNModel {AV18, AV4Pc, N2LO_10, N2LO_12, N3LO_600};

class Nuclear_Info
{
 public:
  Nuclear_Info(int thisZ, int thisN, char* uType);
  Nuclear_Info(int thisZ, int thisN, NNModel uType);
  ~Nuclear_Info();
  double get_S(double krel, int l_type, int r_type);
  double get_pp(double k_rel);
  double get_nn(double k_rel);
  double get_pn(double k_rel);
  double get_pn0(double k_rel);
  double get_pn1(double k_rel);
  NNModel get_InteractionType();
  double get_mA();
  double get_mAm2();
  double get_sigmaCM();
  double get_sigmaE();
  double get_Estar();
  double get_Cnn0();
  double get_Cpp0();
  double get_Cpn0();
  double get_Cpn1();

  void set_Nucleus(int thisZ, int thisN);
  void setCustomValues(double newSigma, double newEstar, double newCpp0, double Cnn0, double newCpn0, double newCpn1);
  void set_Interaction(NNModel thisNNType);
  void set_Interaction(char* thisNNType);
  void set_sigmaCM(double newSigma);
  void set_Estar(double newEstar);
  void set_Cpp0(double newCpp0);
  void set_Cnn0(double newCnn0);
  void set_Cpn0(double newCpn0);
  void set_Cpn1(double newCpn1);
  void set_sigmaE(double newSigE);
  
 private:
  int Z;
  int N;
  int A;
  NNModel u;
  double mA;
  double mAm2;
  double sigmaCM;
  double phiSq_pp0[100];
  double phiSq_nn0[100];
  double phiSq_pn0[100];
  double phiSq_pn1[100];
  double Cpp0;
  double Cnn0;
  double Cpn0;
  double Cpn1;
  double Estar;
  double sigmaE;
  
  double get_phiSq(double *phiPtr, double k_rel);
  void fill_arrays_AV18();
  void fill_arrays_n2lo_local();
  void fill_arrays_n3lo_nonlocal();
  void fill_arrays_n2lo_12_local();
  void fill_arrays_AV4Pc();
  
};

#endif
