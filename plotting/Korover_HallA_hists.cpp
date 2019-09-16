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
  if (argc < 5)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tKorover_HallA_sim /path/to/sim/file/1 /path/to/sim/file/2 /path/to/sim/file/3 /path/to/out/file\n\n";
	return -1;
    }

  TFile * infile_1 = new TFile(argv[1]);
  TFile * infile_2 = new TFile(argv[2]);
  TFile * infile_3 = new TFile(argv[3]);

  TFile * outfile = new TFile(argv[4],"RECREATE");

  //Constants
  const double Ebeam = 4.454;
  const TVector3 vBeam(Ebeam,0.,0.);

  // Settings
  bool verbose = false;
  bool rand_flag = false;
  bool use_lc = false;
  TRandom3 myRand(0);

  int c;
  while ((c=getopt (argc-5, &argv[5], "vrCl")) != -1)
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
  TH1D * hpn_mMiss = new TH1D("epn_Mmiss","epn;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpn_mMiss);
  TH1D * hpp_mMiss = new TH1D("epp_Mmiss","epp;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpp_mMiss);

  TH1D * hpn_cosgamma_1 = new TH1D("epn_cosgamma_1","epn;cos gamma;Counts",25*4,-1.,-0.94);
  h_list.push_back(hpn_cosgamma_1);
  TH1D * hpn_mMiss_1 = new TH1D("epn_Mmiss_1","epn;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpn_mMiss_1);
  TH1D * hpp_mMiss_1 = new TH1D("epp_Mmiss_1","epp;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpp_mMiss_1);
  
  TH1D * hpn_cosgamma_2 = new TH1D("epn_cosgamma_2","epn;cos gamma;Counts",25*4,-1.,-0.94);
  h_list.push_back(hpn_cosgamma_2);
  TH1D * hpn_mMiss_2 = new TH1D("epn_Mmiss_2","epn;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpn_mMiss_2);
  TH1D * hpp_mMiss_2 = new TH1D("epp_Mmiss_2","epp;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpp_mMiss_2);
  
  TH1D * hpn_cosgamma_3 = new TH1D("epn_cosgamma_3","epn;cos gamma;Counts",25*4,-1.,-0.94);
  h_list.push_back(hpn_cosgamma_3);
  TH1D * hpn_mMiss_3 = new TH1D("epn_Mmiss_3","epn;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpn_mMiss_3);
  TH1D * hpp_mMiss_3 = new TH1D("epp_Mmiss_3","epp;mMiss [GeV];Counts",40*4,1.7,2.2);
  h_list.push_back(hpp_mMiss_3);

  for (int i=0; i<h_list.size(); i++)
    h_list[i]->Sumw2();

  // Tree Variable initialization
  Float_t Pe[3], Pp[2][3];
  Double_t weightpp, weightpn;
  Double_t thetak_gen, phik_gen, thetaRel_gen, phiRel_gen, pRel_gen;
  
  
  // Input Tree (1)
  TTree * inTree_1 = (TTree*)infile_1->Get("T");
  
  inTree_1->SetBranchAddress("Pe",Pe);
  inTree_1->SetBranchAddress("Pp",Pp);
  if (!use_lc)
    {
      inTree_1->SetBranchAddress("weightpp",&weightpp);
      inTree_1->SetBranchAddress("weightpn",&weightpn);
    }
  else
    {
      inTree_1->SetBranchAddress("lcweightpp",&weightpp);
      inTree_1->SetBranchAddress("lcweightpn",&weightpn);
    }
  inTree_1->SetBranchAddress("thetak_gen",&thetak_gen);
  inTree_1->SetBranchAddress("phik_gen",&phik_gen);
  inTree_1->SetBranchAddress("thetaRel_gen",&thetaRel_gen);
  inTree_1->SetBranchAddress("phiRel_gen",&phiRel_gen);
  inTree_1->SetBranchAddress("pRel_gen",&pRel_gen);

  const int nEvents_1 = inTree_1->GetEntries();
  for (int event=0 ; event < nEvents_1 ; event++)
    {

      inTree_1->GetEvent(event);

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
      TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
      
      TVector3 vq = vBeam - ve;
      TVector3 vmiss = vlead - vq;
      TVector3 vcm = vmiss + vrec;

      double cosgamma = cos(vmiss.Angle(vrec));
      hpn_cosgamma_1->Fill(cosgamma,weightpn);

      double Elead = sqrt(vlead.Mag2() + sq(mN));
      double Erec = sqrt(vrec.Mag2() + sq(mN));
      double nu = Ebeam - ve.Mag();
      double m_miss = sqrt(sq(nu + m_4He - Elead - Erec) - vcm.Mag2());
      hpn_mMiss_1->Fill(m_miss,weightpn);
      hpp_mMiss_1->Fill(m_miss,weightpp);
    }

  infile_1->Close();

  // Input Tree (2)
  TTree * inTree_2 = (TTree*)infile_2->Get("T");
  
  inTree_2->SetBranchAddress("Pe",Pe);
  inTree_2->SetBranchAddress("Pp",Pp);
  if (!use_lc)
    {
      inTree_2->SetBranchAddress("weightpp",&weightpp);
      inTree_2->SetBranchAddress("weightpn",&weightpn);
    }
  else
    {
      inTree_2->SetBranchAddress("lcweightpp",&weightpp);
      inTree_2->SetBranchAddress("lcweightpn",&weightpn);
    }
  inTree_2->SetBranchAddress("thetak_gen",&thetak_gen);
  inTree_2->SetBranchAddress("phik_gen",&phik_gen);
  inTree_2->SetBranchAddress("thetaRel_gen",&thetaRel_gen);
  inTree_2->SetBranchAddress("phiRel_gen",&phiRel_gen);
  inTree_2->SetBranchAddress("pRel_gen",&pRel_gen);

  const int nEvents_2 = inTree_2->GetEntries();
  for (int event=0 ; event < nEvents_2 ; event++)
    {

      inTree_2->GetEvent(event);

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
      TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
      
      TVector3 vq = vBeam - ve;
      TVector3 vmiss = vlead - vq;
      TVector3 vcm = vmiss + vrec;

      double cosgamma = cos(vmiss.Angle(vrec));
      hpn_cosgamma_2->Fill(cosgamma,weightpn);
	
      double Elead = sqrt(vlead.Mag2() + sq(mN));
      double Erec = sqrt(vrec.Mag2() + sq(mN));
      double nu = Ebeam - ve.Mag();
      double m_miss = sqrt(sq(nu + m_4He - Elead - Erec) - vcm.Mag2());
      hpn_mMiss_2->Fill(m_miss,weightpn);
      hpp_mMiss_2->Fill(m_miss,weightpp);
    }

  infile_2->Close();

  // Input Tree (3)
  TTree * inTree_3 = (TTree*)infile_3->Get("T");
  
  inTree_3->SetBranchAddress("Pe",Pe);
  inTree_3->SetBranchAddress("Pp",Pp);
  if (!use_lc)
    {
      inTree_3->SetBranchAddress("weightpp",&weightpp);
      inTree_3->SetBranchAddress("weightpn",&weightpn);
    }
  else
    {
      inTree_3->SetBranchAddress("lcweightpp",&weightpp);
      inTree_3->SetBranchAddress("lcweightpn",&weightpn);
    }
  inTree_3->SetBranchAddress("thetak_gen",&thetak_gen);
  inTree_3->SetBranchAddress("phik_gen",&phik_gen);
  inTree_3->SetBranchAddress("thetaRel_gen",&thetaRel_gen);
  inTree_3->SetBranchAddress("phiRel_gen",&phiRel_gen);
  inTree_3->SetBranchAddress("pRel_gen",&pRel_gen);

  const int nEvents_3 = inTree_3->GetEntries();
  for (int event=0 ; event < nEvents_3 ; event++)
    {

      inTree_3->GetEvent(event);

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
      TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
      
      TVector3 vq = vBeam - ve;
      TVector3 vmiss = vlead - vq;
      TVector3 vcm = vmiss + vrec;

      double cosgamma = cos(vmiss.Angle(vrec));
      hpn_cosgamma_3->Fill(cosgamma,weightpn);
	
      double Elead = sqrt(vlead.Mag2() + sq(mN));
      double Erec = sqrt(vrec.Mag2() + sq(mN));
      double nu = Ebeam - ve.Mag();
      double m_miss = sqrt(sq(nu + m_4He - Elead - Erec) - vcm.Mag2());
      hpn_mMiss_3->Fill(m_miss,weightpn);
      hpp_mMiss_3->Fill(m_miss,weightpp);
    }

  infile_3->Close();
  
  double N_pn_1 = 107.;
  double N_pn_2 = 66.;
  double N_pn_3 = 50.;
  double N_pp_1 = 7.5;
  double N_pp_2 = 20.;
  double N_pp_3 = 22.5;

  double N_pn = 113.;
  double N_pp = 37.;
  
  hpn_cosgamma_1->Scale(N_pn_1/hpn_cosgamma_1->Integral());
  hpn_cosgamma_2->Scale(N_pn_2/hpn_cosgamma_2->Integral());
  hpn_cosgamma_3->Scale(N_pn_3/hpn_cosgamma_3->Integral());
  
  hpn_mMiss_1->Scale(N_pn_1/hpn_mMiss_1->Integral());
  hpn_mMiss_2->Scale(N_pn_2/hpn_mMiss_2->Integral());
  hpn_mMiss_3->Scale(N_pn_3/hpn_mMiss_3->Integral());

  hpp_mMiss_1->Scale(N_pp_1/hpp_mMiss_1->Integral());
  hpp_mMiss_2->Scale(N_pp_2/hpp_mMiss_2->Integral());
  hpp_mMiss_3->Scale(N_pp_3/hpp_mMiss_3->Integral());

  hpn_cosgamma->Add(hpn_cosgamma_1);
  hpn_cosgamma->Add(hpn_cosgamma_2);
  hpn_cosgamma->Add(hpn_cosgamma_3);
  hpn_cosgamma->Scale(N_pn/hpn_cosgamma->Integral());

  hpn_mMiss->Add(hpn_mMiss_1);
  hpn_mMiss->Add(hpn_mMiss_2);
  hpn_mMiss->Add(hpn_mMiss_3);
  hpn_mMiss->Scale(N_pn/hpn_mMiss->Integral());
  
  hpp_mMiss->Add(hpp_mMiss_1);
  hpp_mMiss->Add(hpp_mMiss_2);
  hpp_mMiss->Add(hpp_mMiss_3);
  hpp_mMiss->Scale(N_pp/hpp_mMiss->Integral());
  
  outfile->cd();

  for (int i=0; i<h_list.size(); i++)
    h_list[i]->Write();
  
  outfile->Close();
  
  return 0;
}
