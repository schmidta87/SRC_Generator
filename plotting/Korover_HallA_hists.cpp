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
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TVector3.h"

#include "Nuclear_Info.h"
#include "constants.h"
#include "helpers.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tKorover_HallA_sim /path/to/sim/file/(e,e'p) /path/to/sim/file/(e,e'pN) /path/to/out/file\n\n";
	return -1;
    }

  TFile * infile_full = new TFile(argv[1]);
  TFile * infile_sub = new TFile(argv[2]);

  TFile * outfile = new TFile(argv[3],"RECREATE");

  //Constants
  const double Ebeam = 4.454;
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
  
  TH1D * hpn_cosgamma = new TH1D("epn_cosgamma","epn;cos gamma;Counts",25*4,-1.,-0.94);
  h_list.push_back(hpn_cosgamma);
  int mMiss_bins = 40*4;
  TH1D * hpn_mMiss = new TH1D("epn_Mmiss","epn;mMiss [GeV];Counts",mMiss_bins,1.7,2.2);
  h_list.push_back(hpn_mMiss);
  TH1D * hpp_mMiss = new TH1D("epp_Mmiss","epp;mMiss [GeV];Counts",mMiss_bins,1.7,2.2);
  h_list.push_back(hpp_mMiss);
  int Pm_bins = 40*4;
  TH1D * hp_Pm = new TH1D("ep_Pm","ep;pMiss [GeV];Counts",Pm_bins,0.2,1.);
  h_list.push_back(hp_Pm);
  TH1D * hpn_Pm = new TH1D("epn_Pm","epn;pMiss [GeV];Counts",Pm_bins,0.2,1.);
  h_list.push_back(hpn_Pm);
  TH1D * hpp_Pm = new TH1D("epp_Pm","epp;pMiss [GeV];Counts",Pm_bins,0.2,1.);
  h_list.push_back(hpp_Pm);

  TH1D * hp_setting = new TH1D("ep_setting","ep;pMiss [MeV];Counts",3,437.5,812.5);
  h_list.push_back(hp_setting);
  TH1D * hp_n_setting = new TH1D("ep_n_setting","ep(n);pMiss [MeV];Counts",3,437.5,812.5);
  h_list.push_back(hp_n_setting);
  TH1D * hp_p_setting = new TH1D("ep_p_setting","ep(p);pMiss [MeV];Counts",3,437.5,812.5);
  h_list.push_back(hp_p_setting);
  TH1D * hpn_setting = new TH1D("epn_setting","epn;pMiss [MeV];Counts",3,437.5,812.5);
  h_list.push_back(hpn_setting);
  TH1D * hpp_setting = new TH1D("epp_setting","epp;pMiss [MeV];Counts",3,437.5,812.5);
  h_list.push_back(hpp_setting);
  
  TH1D * hpn_cosgamma_set[3];
  TH1D * hpn_mMiss_set[3];
  TH1D * hpp_mMiss_set[3];
  TH1D * hp_Pm_set[3];
  TH1D * hpn_Pm_set[3];
  TH1D * hpp_Pm_set[3];

  for (int i = 0; i<3; i++)
    {
      char temp[100];

      int set;
      switch(i)
	{
	case 0:
	  set = 500;
	  break;
	case 1:
	  set = 625;
	  break;
	case 2:
	  set = 750;
	  break;
	default:
	  abort();
	    }

      sprintf(temp,"epn_cosgamma_%i",set);
      hpn_cosgamma_set[i] = new TH1D(temp,"epn;cos gamma;Counts",25*4,-1.,-0.94);
      h_list.push_back(hpn_cosgamma_set[i]);

      sprintf(temp,"epn_Mmiss_%i",set);
      hpn_mMiss_set[i] = new TH1D(temp,"epn;mMiss [GeV];Counts",mMiss_bins,1.7,2.2);
      h_list.push_back(hpn_mMiss_set[i]);

      sprintf(temp,"epp_Mmiss_%i",set);
      hpp_mMiss_set[i] = new TH1D(temp,"epp;mMiss [GeV];Counts",mMiss_bins,1.7,2.2);
      h_list.push_back(hpp_mMiss_set[i]);

      sprintf(temp,"ep_Pm_%i",set);
      hp_Pm_set[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",Pm_bins,0.2,1.);
      h_list.push_back(hp_Pm_set[i]);

      sprintf(temp,"epn_Pm_%i",set);
      hpn_Pm_set[i] = new TH1D(temp,"epn;pMiss [GeV];Counts",Pm_bins,0.2,1.);
      h_list.push_back(hpn_Pm_set[i]);

      sprintf(temp,"epp_Pm_%i",set);
      hpp_Pm_set[i] = new TH1D(temp,"epp;pMiss [GeV];Counts",Pm_bins,0.2,1.);
      h_list.push_back(hpp_Pm_set[i]);

    }
  
  for (int i=0; i<h_list.size(); i++)
    h_list[i]->Sumw2();

  // Acceptance Factors:
  TGraphAsymmErrors * pn_acc = new TGraphAsymmErrors();
  pn_acc->SetName("pn_acc");
  pn_acc->SetTitle("pn_acc;p_miss [MeV]; pn acceptance");
  TGraphAsymmErrors * pp_acc = new TGraphAsymmErrors();
  pp_acc->SetName("pp_acc");
  pp_acc->SetTitle("pp_acc;p_miss [MeV]; pp acceptance");
	  
  // Tree Variable initialization
  Float_t Pe[3], Pp[2][3];
  Double_t weightp, weightpN;
  Double_t thetak_gen, phik_gen, thetaRel_gen, phiRel_gen, pRel_gen;
  Int_t setting, rec_code;
  
  // Input Tree (full)
  TTree * inTree_full = (TTree*)infile_full->Get("T");
  
  inTree_full->SetBranchAddress("Pe",Pe);
  inTree_full->SetBranchAddress("Pp",Pp);
  if (!use_lc)
    {
      inTree_full->SetBranchAddress("weightp",&weightp);
      inTree_full->SetBranchAddress("weightpN",&weightpN);
    }
  else
    {
      inTree_full->SetBranchAddress("lcweightp",&weightp);
      inTree_full->SetBranchAddress("lcweightpN",&weightpN);
    }
  inTree_full->SetBranchAddress("thetak_gen",&thetak_gen);
  inTree_full->SetBranchAddress("phik_gen",&phik_gen);
  inTree_full->SetBranchAddress("thetaRel_gen",&thetaRel_gen);
  inTree_full->SetBranchAddress("phiRel_gen",&phiRel_gen);
  inTree_full->SetBranchAddress("pRel_gen",&pRel_gen);
  inTree_full->SetBranchAddress("setting",&setting);
  inTree_full->SetBranchAddress("rec_code",&rec_code);
  
  const int nEvents_full = inTree_full->GetEntries();
  for (int event=0 ; event < nEvents_full ; event++)
    {

      inTree_full->GetEvent(event);

      int set_bin;
      switch(setting)
	{
	case 500:
	  set_bin = 0;
	  break;
	case 625:
	  set_bin = 1;
	  break;
	case 750:
	  set_bin = 2;
	  break;
	default:
	  abort();
	    }

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
      
      TVector3 vq = vBeam - ve;
      TVector3 vmiss = vlead - vq;
      double pmiss = vmiss.Mag();
      hp_Pm_set[set_bin]->Fill(pmiss,weightp);
      hp_setting->Fill(setting,weightp);
      
      if (rec_code == pCode)
	{
	  hp_p_setting->Fill(setting,weightp);
	}
      else if (rec_code == nCode)
	{
	  hp_n_setting->Fill(setting,weightp);
	}
      else
	abort();
      
    }

  infile_full->Close();

  // Input Tree (sub)
  TTree * inTree_sub = (TTree*)infile_sub->Get("T");
  
  inTree_sub->SetBranchAddress("Pe",Pe);
  inTree_sub->SetBranchAddress("Pp",Pp);
  if (!use_lc)
    {
      inTree_sub->SetBranchAddress("weightp",&weightp);
      inTree_sub->SetBranchAddress("weightpN",&weightpN);
    }
  else
    {
      inTree_sub->SetBranchAddress("lcweightp",&weightp);
      inTree_sub->SetBranchAddress("lcweightpN",&weightpN);
    }
  inTree_sub->SetBranchAddress("thetak_gen",&thetak_gen);
  inTree_sub->SetBranchAddress("phik_gen",&phik_gen);
  inTree_sub->SetBranchAddress("thetaRel_gen",&thetaRel_gen);
  inTree_sub->SetBranchAddress("phiRel_gen",&phiRel_gen);
  inTree_sub->SetBranchAddress("pRel_gen",&pRel_gen);
  inTree_sub->SetBranchAddress("setting",&setting);
  inTree_sub->SetBranchAddress("rec_code",&rec_code);
  
  const int nEvents_sub = inTree_sub->GetEntries();
  for (int event=0 ; event < nEvents_sub ; event++)
    {

      inTree_sub->GetEvent(event);

      int set_bin;
      switch(setting)
	{
	case 500:
	  set_bin = 0;
	  break;
	case 625:
	  set_bin = 1;
	  break;
	case 750:
	  set_bin = 2;
	  break;
	default:
	  abort();
	    }

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
      TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
      
      TVector3 vq = vBeam - ve;
      TVector3 vmiss = vlead - vq;
      TVector3 vcm = vmiss + vrec;

      double cosgamma = cos(vmiss.Angle(vrec));

      double Elead = sqrt(vlead.Mag2() + sq(mN));
      double Erec = sqrt(vrec.Mag2() + sq(mN));
      double nu = Ebeam - ve.Mag();
      double m_miss = sqrt(sq(nu + m_4He - Elead - Erec) - vcm.Mag2());

      double pmiss = vmiss.Mag();

      if (rec_code == pCode)
	{
	  hpp_mMiss_set[set_bin]->Fill(m_miss,weightpN);
	  hpp_Pm_set[set_bin]->Fill(pmiss,weightpN);
	  hpp_setting->Fill(setting,weightpN);
	}
      else if (rec_code == nCode)
	{
	  hpn_cosgamma_set[set_bin]->Fill(cosgamma,weightpN);
	  hpn_mMiss_set[set_bin]->Fill(m_miss,weightpN);
	  hpn_Pm_set[set_bin]->Fill(pmiss,weightpN);
	  hpn_setting->Fill(setting,weightpN);
	}
      else
	abort();
	  
    }

  infile_sub->Close();

  double Np [] = {3100., 1619., 1228.};

  for (int i = 0; i<2; i++)
    {
      double normp = Np[i]/hp_Pm_set[i]->Integral();
      double normpN = normp;
      
      hpn_cosgamma_set[i]->Scale(normpN);
      hpn_mMiss_set[i]->Scale(normpN);
      hpp_mMiss_set[i]->Scale(normpN);
      hp_Pm_set[i]->Scale(normp);
      hpn_Pm_set[i]->Scale(normpN);
      hpp_Pm_set[i]->Scale(normpN);
      
      hpn_cosgamma->Add(hpn_cosgamma_set[i]);
      hpn_mMiss->Add(hpn_mMiss_set[i]);
      hpp_mMiss->Add(hpp_mMiss_set[i]);
      hp_Pm->Add(hp_Pm_set[i]);
      hpn_Pm->Add(hpn_Pm_set[i]);
      hpp_Pm->Add(hpp_Pm_set[i]);
    }
  
  outfile->cd();

  pn_acc->BayesDivide(hpn_setting,hp_n_setting);
  pp_acc->BayesDivide(hpp_setting,hp_p_setting);
  pn_acc->Write();
  pp_acc->Write();
  
  for (int i=0; i<h_list.size(); i++)
    h_list[i]->Write();
  
  outfile->Close();

  return 0;
}
