#ifndef __NUCLEAR_INFO_H__
#define __NUCLEAR_INFO_H__
#include <vector>

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
  double get_mAmpp();
  double get_mAmpn();
  double get_mAmnn();
  double get_sigmaCM();
  double get_sigmaE();
  double get_Estar();
  double get_Cnn0();
  double get_Cpp0();
  double get_Cpn0();
  double get_Cpn1();
  void randomize();
  
  void do_SCX(int &lead_type, int &rec_type, double r);
  std::vector<double> get_SCX_Ps();

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
  void set_SCX_Ps(double pPP2NP_new,
		 double pPP2PN_new,
		 double pPP2NN_new,
		 double pPN2NN_new,
		 double pPN2PP_new,
		 double pPN2NP_new,
		 double pNP2PP_new,
		 double pNP2NN_new,
		 double pNP2PN_new,
		 double pNN2PN_new,
		 double pNN2NP_new,
		 double pNN2PP_new);
  
 private:
  int Z;
  int N;
  int A;
  NNModel u;
  double mA;
  double mAmpp;
  double mAmpn;
  double mAmnn;
  double sigmaCM;
  double d_sigmaCM;
  double phiSq_pp0[100];
  double phiSq_nn0[100];
  double phiSq_pn0[100];
  double phiSq_pn1[100];
  double Cpp0;
  double d_Cpp0;
  double Cnn0;
  double d_Cnn0;
  double Cpn0;
  double d_Cpn0;
  double Cpn1;
  double d_Cpn1;
  double Estar;
  double Estar_max;
  double sigmaE;
  
  double get_phiSq(double *phiPtr, double k_rel);

  bool set_Contacts_SS_r();
  bool set_Contacts_SS_k();
  bool set_Contacts_EG2();
  
  void fill_arrays_AV18();
  void fill_arrays_n2lo_local();
  void fill_arrays_n3lo_nonlocal();
  void fill_arrays_n2lo_12_local();
  void fill_arrays_AV4Pc();

  double pPP2NP;
  double d_pPP2NP;
  double pPP2PN;
  double d_pPP2PN;
  double pPP2NN;
  double d_pPP2NN;

  double pPN2NN;
  double d_pPN2NN;
  double pPN2PP;
  double d_pPN2PP;
  double pPN2NP;
  double d_pPN2NP;

  double pNP2PP;
  double d_pNP2PP;
  double pNP2NN;
  double d_pNP2NN;
  double pNP2PN;
  double d_pNP2PN;

  double pNN2PN;
  double d_pNN2PN;
  double pNN2NP;
  double d_pNN2NP;
  double pNN2PP;
  double d_pNN2PP;
  
};

#endif
