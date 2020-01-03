#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TVector3.h"

#include "Nuclear_Info.h"
#include "constants.h"
#include "helpers.h"

using namespace std;

double eta_pp = 0.85;

double Tp;
double Tpp;

int main(int argc, char ** argv)
{
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tShneor_HallA_sim /path/to/sim/file/(e,e'p) /path/to/sim/file/(e,e'pp) /path/to/out/file\n\n";
	return -1;
    }

  TFile * infile_full = new TFile(argv[1]);
  TFile * infile_sub = new TFile(argv[2]);

  TFile * outfile = new TFile(argv[3],"RECREATE");

  //Constants
  const double Ebeam = 4.627;
  const TVector3 vBeam(Ebeam,0.,0.);

  // Settings
  bool verbose = false;
  bool rand_flag = false;
  bool use_lc = false;
  TRandom3 myRand(0);

  int c;
  while ((c=getopt (argc-4, &argv[4], "vrCl")) != -1)
    switch(c)
      {
      case 'v':
	verbose = true;
	break;
      case 'r':
	rand_flag = true;
	break;
      case 'l':
	use_lc = false;
      case '?':
	return -1;
      default:
	abort();
      }

  // Create histograms
  vector<TH1*> h_list;
  
  TH1D * hp_setting = new TH1D("ep_setting","ep;pMiss [MeV];Counts",3,300,600);
  h_list.push_back(hp_setting);
  TH1D * hp_p_setting = new TH1D("ep_p_setting","ep(p);pMiss [MeV];Counts",3,300,600);
  h_list.push_back(hp_p_setting);
  TH1D * hpp_setting = new TH1D("epp_setting","epp;pMiss [MeV];Counts",3,300,600);
  h_list.push_back(hpp_setting);

  TH1D * hpp_cosgamma = new TH1D("epp_cosgamma","epp;cos gamma;Counts",8*4,-1.,-0.90);
  h_list.push_back(hpp_cosgamma);
  TH1D * hp_Em = new TH1D("ep_Em","ep;EMiss [GeV];Counts",150,0.0,0.35);
  h_list.push_back(hp_Em);

  TH1D * hpp_cosgamma_set[3];
  TH1D * hp_Em_set[3];
  TH1D * hp_omega_set[3];
  TH1D * hp_q_set[3];
  TH1D * hp_xB_set[3];
  TH1D * hp_thetaM_set[3];
  TH1D * hp_thetaRec_set[3];
  TH1D * hpp_thetaM_set[3];
  TH1D * hp_Pm_set[3];
  TH1D * hpp_Pm_set[3];
  TH2D * hp_Em_thetaM_set[3];
  TH2D * hp_phiM_thetaM_set[3];
  TH2D * hp_phiRec_thetaRec_set[3];
  TH2D * hp_omega_thetaM_set[3];
  TH2D * hp_q_thetaM_set[3];
  TH2D * hp_xB_thetaM_set[3];
  
  for (int i = 0; i<3; i++)
    {
      char temp[100];

      int set;
      switch(i)
	{
	  case 0:
	    set = 350;
	    break;
	  case 1:
	    set = 450;
	    break;
	  case 2:
	    set = 550;
	    break;
	  default:
	    abort();
	}

      sprintf(temp,"epp_cosgamma_%i",set);
      hpp_cosgamma_set[i] = new TH1D(temp,"epp;cos gamma;Counts",8*4,-1.,-.90);
      h_list.push_back(hpp_cosgamma_set[i]);
      
      sprintf(temp,"ep_Em_%i",set);
      hp_Em_set[i] = new TH1D(temp,"ep;Emiss [GeV];Counts",150,0.0,0.350);
      h_list.push_back(hp_Em_set[i]);

      sprintf(temp,"ep_omega_%i",set);
      hp_omega_set[i] = new TH1D(temp,"ep;omega [GeV];Counts",40,0.7,1.1);
      h_list.push_back(hp_omega_set[i]);

      sprintf(temp,"ep_q_%i",set);
      hp_q_set[i] = new TH1D(temp,"ep;q [GeV];Counts",40,1.5,1.8);
      h_list.push_back(hp_q_set[i]);

      sprintf(temp,"ep_xB_%i",set);
      hp_xB_set[i] = new TH1D(temp,"ep;xB;Counts",40,0.7,1.5);
      h_list.push_back(hp_xB_set[i]);

      sprintf(temp,"ep_thetaM_%i",set);
      hp_thetaM_set[i] = new TH1D(temp,"ep;thetaM [degrees];Counts",100,40.,130.);
      h_list.push_back(hp_thetaM_set[i]);

      sprintf(temp,"ep_thetaRec_%i",set);
      hp_thetaRec_set[i] = new TH1D(temp,"ep;thetaRec [degrees];Counts",100,40.,130.);
      h_list.push_back(hp_thetaRec_set[i]);

      sprintf(temp,"epp_thetaM_%i",set);
      hpp_thetaM_set[i] = new TH1D(temp,"epp;thetaM [degrees];Counts",100,40.,130.);
      h_list.push_back(hpp_thetaM_set[i]);

      sprintf(temp,"ep_Pm_%i",set);
      hp_Pm_set[i] = new TH1D(temp,"ep;Pm [GeV];Counts",60,0.1,0.7);
      h_list.push_back(hp_Pm_set[i]);

      sprintf(temp,"epp_Pm_%i",set);
      hpp_Pm_set[i] = new TH1D(temp,"epp;Pm [GeV];Counts",60,0.1,0.7);
      h_list.push_back(hpp_Pm_set[i]);

      sprintf(temp,"ep_Em_thetaM_%i",set);
      hp_Em_thetaM_set[i] = new TH2D(temp,"epp;Em [GeV]; thetaM [degrees]",150,0.0,0.350,100,40.,130.);
      h_list.push_back(hp_Em_thetaM_set[i]);      

      sprintf(temp,"ep_phiM_thetaM_%i",set);
      hp_phiM_thetaM_set[i] = new TH2D(temp,"epp;phiM [degrees]; thetaM [degrees]",40,-60,60,100,40.,130.);
      h_list.push_back(hp_phiM_thetaM_set[i]);      

      sprintf(temp,"ep_omega_thetaM_%i",set);
      hp_omega_thetaM_set[i] = new TH2D(temp,"epp; omega [GeV]; thetaM [degrees]",40,0.7,1.1,100,40.,130.);
      h_list.push_back(hp_omega_thetaM_set[i]);
      
      sprintf(temp,"ep_q_thetaM_%i",set);
      hp_q_thetaM_set[i] = new TH2D(temp,"epp; q [geV]; thetaM [degrees]",40,0.5,0.8,100,40.,130.);
      h_list.push_back(hp_q_thetaM_set[i]);
      
      sprintf(temp,"ep_xB_thetaM_%i",set);
      hp_xB_thetaM_set[i] = new TH2D(temp,"epp; xB; thetaM [degrees]",40,0.7,1.5,100,40.,130.);
      h_list.push_back(hp_xB_thetaM_set[i]);      

      sprintf(temp,"ep_phiRec_thetaRec_%i",set);
      hp_phiRec_thetaRec_set[i] = new TH2D(temp,"epp;phiRec [degrees]; thetaRec [degrees]",24,-60,60,24,39.,159.);
      h_list.push_back(hp_phiRec_thetaRec_set[i]);      

    }
  
  for (int i=0; i<h_list.size(); i++)
    h_list[i]->Sumw2();

  // Acceptance Factors:
  TGraphAsymmErrors * pp_acc = new TGraphAsymmErrors();
  pp_acc->SetName("pp_acc");
  pp_acc->SetTitle("pp_acc;p_miss [MeV]; pp acceptance");

  // Ratios:
  TGraphAsymmErrors * pp_to_p = new TGraphAsymmErrors();
  pp_to_p->SetName("pp_to_p");
  pp_to_p->SetTitle("pp_to_p;p_miss setting; pp/p");

  TGraphAsymmErrors * pp_to_p_cor = new TGraphAsymmErrors();
  pp_to_p_cor->SetName("pp_to_p_cor");
  pp_to_p_cor->SetTitle("pp_to_p_cor;p_miss setting; pp/p");

  // Tree Variable initialization
  Float_t Pe[3], Pp[2][3];
  Double_t weightp, weightpp;
  Double_t thetak_gen, phik_gen, thetaRel_gen, phiRel_gen, pRel_gen;
  Int_t setting, rec_code;
  
  // Input Tree (full)
  TTree * inTree_full = (TTree*)infile_full->Get("T");
  
  inTree_full->SetBranchAddress("Pe",Pe);
  inTree_full->SetBranchAddress("Pp",Pp);
  if (!use_lc)
    {
      inTree_full->SetBranchAddress("weightp",&weightp);
      inTree_full->SetBranchAddress("weightpp",&weightpp);
    }
  else
    {
      inTree_full->SetBranchAddress("lcweightp",&weightp);
      inTree_full->SetBranchAddress("lcweightpp",&weightpp);
    }
  inTree_full->SetBranchAddress("setting",&setting);
  inTree_full->SetBranchAddress("rec_code",&rec_code);
  
  const int nEvents_full = inTree_full->GetEntries();
  for (int event=0 ; event < nEvents_full ; event++)
    {

      inTree_full->GetEvent(event);

      int set_bin;
      switch(setting)
	{
	case 350:
	  set_bin = 0;
	  Tp = 0.79;
	  Tpp = Tp*0.84;
	  break;
	case 450:
	  set_bin = 1;
	  Tp = 0.77;
	  Tpp = Tp*0.82;
	  break;
	case 550:
	  set_bin = 2;
	  Tp = 0.76;
	  Tpp = Tp*0.81;
	  break;
	default:
	  abort();
	    }

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
      TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
      
      TVector3 vq = vBeam - ve;
      TVector3 vmiss = vq - vlead;
      double pmiss = vmiss.Mag();
      hp_setting->Fill(setting,weightp);
      
      // Missing Energy Definition from Thesis
      double omega = Ebeam - ve.Mag();
      double Em = mN - m_12C + sqrt(sq(omega + m_12C - sqrt(vlead.Mag2() + sq(mN))) - sq(pmiss));
      hp_Em_set[set_bin]->Fill(Em,weightp);

      double QSq = vq.Mag2() - sq(omega);
      double xB = QSq/(2*mN*omega);
	
      hp_omega_set[set_bin]->Fill(omega,weightp);
      hp_q_set[set_bin]->Fill(vq.Mag(),weightp);
      hp_xB_set[set_bin]->Fill(xB,weightp);
      hp_thetaM_set[set_bin]->Fill(vmiss.Phi()*180./M_PI,weightp);  
      hp_Pm_set[set_bin]->Fill(vmiss.Mag(),weightp);  

      hp_Em_thetaM_set[set_bin]->Fill(Em,vmiss.Phi()*180./M_PI,weightp);  

      hp_phiM_thetaM_set[set_bin]->Fill(vmiss.Theta()*180./M_PI-90.,vmiss.Phi()*180./M_PI,weightp);  

      hp_omega_thetaM_set[set_bin]->Fill(omega,vmiss.Phi()*180./M_PI,weightp);  
      hp_q_thetaM_set[set_bin]->Fill(vq.Mag(),vmiss.Phi()*180./M_PI,weightp);  
      hp_xB_thetaM_set[set_bin]->Fill(xB,vmiss.Phi()*180./M_PI,weightp);  

      if (rec_code == pCode)
	{
	  hp_p_setting->Fill(setting,weightp*Tpp/Tp*eta_pp);
	}
      
      hp_thetaRec_set[set_bin]->Fill(vrec.Phi()*180./M_PI,weightp);  

      hp_phiRec_thetaRec_set[set_bin]->Fill(vrec.Theta()*180./M_PI-90.,vrec.Phi()*180./M_PI,weightp);  

    }

  infile_full->Close();

  // Input Tree (sub)
  TTree * inTree_sub = (TTree*)infile_sub->Get("T");
  
  inTree_sub->SetBranchAddress("Pe",Pe);
  inTree_sub->SetBranchAddress("Pp",Pp);
  if (!use_lc)
    {
      inTree_sub->SetBranchAddress("weightp",&weightp);
      inTree_sub->SetBranchAddress("weightpp",&weightpp);
    }
  else
    {
      inTree_sub->SetBranchAddress("lcweightp",&weightp);
      inTree_sub->SetBranchAddress("lcweightpp",&weightpp);
    }
  inTree_sub->SetBranchAddress("setting",&setting);
  inTree_sub->SetBranchAddress("rec_code",&rec_code);
  
  const int nEvents_sub = inTree_sub->GetEntries();
  for (int event=0 ; event < nEvents_sub ; event++)
    {

      inTree_sub->GetEvent(event);

      if (weightpp <= 0.)
	continue;
      
      int set_bin;
      switch(setting)
	{
	case 350:
	  set_bin = 0;
	  Tp = 0.79;
	  Tpp = Tp*0.84;
	  break;
	case 450:
	  set_bin = 1;
	  Tp = 0.77;
	  Tpp = Tp*0.82;
	  break;
	case 550:
	  set_bin = 2;
	  Tp = 0.76;
	  Tpp = Tp*0.81;
	  break;
	default:
	  abort();
	    }

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
      TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);

      TVector3 vq = vBeam - ve;
      TVector3 vmiss = vq - vlead;
      TVector3 vcm = vmiss + vrec;

      double cosgamma = cos(vmiss.Angle(-vrec));

      if (rec_code == pCode)
	{
	  hpp_cosgamma_set[set_bin]->Fill(cosgamma,weightpp);
	  hpp_setting->Fill(setting,weightpp);
	  hpp_thetaM_set[set_bin]->Fill(vmiss.Phi()*180./M_PI,weightpp);  
	  hpp_Pm_set[set_bin]->Fill(vmiss.Mag(),weightpp);  
	}
	  
    }

  infile_sub->Close();

  double Np [] = {142000., 123000., 87000.};

  double normp;
  double normpp;

  for (int i = 0; i<3; i++)
    {
      normp = Np[i]/hp_Em_set[i]->Integral();
      normpp = normp;
      
      hp_Em_set[i]->Scale(normp);      
      hp_omega_set[i]->Scale(normp);      
      hp_q_set[i]->Scale(normp);      
      hp_xB_set[i]->Scale(normp);      
      hp_thetaM_set[i]->Scale(normp);      
      hp_thetaRec_set[i]->Scale(normp);      
      hp_Pm_set[i]->Scale(normp);      

      hpp_cosgamma_set[i]->Scale(normpp);
      hpp_thetaM_set[i]->Scale(normpp);      
      hpp_Pm_set[i]->Scale(normpp);      

    }

  for (int i = 1; i<3; i++)
    {
      hpp_cosgamma->Add(hpp_cosgamma_set[i]);
      hp_Em->Add(hp_Em_set[i]);
    }  

  outfile->cd();

  pp_acc->BayesDivide(hpp_setting,hp_p_setting);
  pp_acc->Write();

  pp_to_p->BayesDivide(hpp_setting,hp_setting);
  pp_to_p->Write();

  pp_to_p_cor->BayesDivide(hp_p_setting,hp_setting);
  //pp_to_p_cor->Scale(Tp/Tpp*1/eta_pp);
  pp_to_p_cor->Write();

  for (int i=0; i<h_list.size(); i++)
    h_list[i]->Write();
  
  outfile->Close();

  return 0;
}
