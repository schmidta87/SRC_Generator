#ifndef __NUCLEAR_INFO_MF_H__
#define __NUCLEAR_INFO_MF_H__

class Nuclear_Info_MF
{
 public:
  Nuclear_Info_MF(int thisA);
  ~Nuclear_Info_MF();
  double get_mA();
  double get_mAmp();
  double get_mAmn();
  double get_Estar();
  double get_sigmaE();
  double get_Pp(double k, double E);
  double get_Pn(double k, double E);

  void set_Estar(double newEstar);
  void set_sigmaE(double newSigE);
  
 private:
  int A;
  double mA;
  double mAmp;
  double mAmn;
  double Estar;
  double sigmaE;
  double Pp[200][1000];
  double Pn[200][1000];

  double get_P(double PMF[200][1000], double k, double E);
  void fill_arrays();
  
};

#endif
