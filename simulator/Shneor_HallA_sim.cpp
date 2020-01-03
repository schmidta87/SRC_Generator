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

#include "Nuclear_Info.h"
#include "constants.h"
#include "detectors.h"
#include "helpers.h"

using namespace std;

bool thetaMcut = true;
bool xcut = false;

bool subcut = false;

int main(int argc, char ** argv)
{
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tShneor_HallA_sim /path/to/gen/file /path/to/out/file\n\n";
	return -1;
    }

  TFile * infile = new TFile(argv[1]);

  double theta_central = 0.5*M_PI;
  double pe_central = 3.724;
  double phie_central = -19.5*M_PI/180.;
  double phirec_central = 99.*M_PI/180.;
  
  double plead_central_350 = 1.45;
  double philead_central_350 = 40.1*M_PI/180.;    

  double plead_central_450 = 1.42;
  double philead_central_450 = 35.8*M_PI/180.;

  double plead_central_550 = 1.36;
  double philead_central_550 = 32.0*M_PI/180.;

  double mom_cut = 4.5e-2;
  double inplane_cut = 28.e-3;
  double outplane_cut = 60.e-3;
  
  double eta = 0.85;
  
  double Tp;
  double Tpp;
  double thetaMmin;
  
  const double Ebeam = 4.627;
  
  bool verbose = false;
  bool rand_flag = false;
  bool doFCuts = true;

  int c;
  while ((c=getopt (argc-2, &argv[2], "vrC")) != -1)
    switch(c)
      {
      case 'v':
	verbose = true;
	break;
      case 'r':
	rand_flag = true;
	break;
      case 'C':
	doFCuts = false;
	break;
      case '?':
	return -1;
      default:
	abort();
      }
  
  //Input Tree
  TTree * inTree = (TTree*)infile->Get("genT");

  Double_t gen_pe[3], gen_pLead[3], gen_pRec[3], gen_weight, gen_lcweight;
  Int_t lead_type, rec_type;
  inTree->SetBranchAddress("lead_type",&lead_type);
  inTree->SetBranchAddress("rec_type",&rec_type);
  inTree->SetBranchAddress("weight",&gen_weight);
  inTree->SetBranchAddress("lcweight",&gen_lcweight);
  inTree->SetBranchAddress("pe",gen_pe);
  inTree->SetBranchAddress("pLead",gen_pLead);
  inTree->SetBranchAddress("pRec",gen_pRec);
  
  // Other set up
  TRandom3 myRand(0);
    
  TFile * outfile = new TFile(argv[2],"RECREATE");
  
  // Output Tree
  TTree * outTree = new TTree("T","Simulated Data Tree");
  Float_t Pe[3], Pp[2][3];
  Double_t weightp, weightpp, lcweightp, lcweightpp;
  Double_t thetak_gen, phik_gen, thetaRel_gen, phiRel_gen, pRel_gen;
  Int_t setting, rec_code;
  
  outTree->Branch("Pe",Pe,"Pe[3]/F");
  outTree->Branch("Pp",Pp,"Pp[2][3]/F");
  outTree->Branch("weightp",&weightp,"weightp/D");
  outTree->Branch("weightpp",&weightpp,"weightpp/D");
  outTree->Branch("lcweightp",&lcweightp,"lcweightp/D");
  outTree->Branch("lcweightpp",&lcweightpp,"lcweightpp/D");
  outTree->Branch("thetak_gen",&thetak_gen,"thetak_gen/D");
  outTree->Branch("phik_gen",&phik_gen,"phik_gen/D");
  outTree->Branch("thetaRel_gen",&thetaRel_gen,"thetaRel_gen/D");
  outTree->Branch("phiRel_gen",&phiRel_gen,"phiRel_gen/D");
  outTree->Branch("pRel_gen",&pRel_gen,"pRel_gen/D");
  outTree->Branch("setting",&setting,"setting/I");
  outTree->Branch("rec_code",&rec_code,"rec_code/I");

  const int nEvents = inTree->GetEntries();
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %1000000==0 and verbose)
	cout << "Working on event " << event << " out of " << nEvents <<"\n";

      inTree->GetEvent(event);
      if (gen_weight == 0. and gen_lcweight == 0.)
	continue;

      setting = 0;
      
      // Require a leading proton
      if (lead_type != pCode)
	continue;

      // Create vectors for the particles; rotated for Hall-A geometry
      TVector3 ve(gen_pe[2],gen_pe[0],gen_pe[1]);
      TVector3 vlead(gen_pLead[2],gen_pLead[0],gen_pLead[1]);
      TVector3 vrec(gen_pRec[2],gen_pRec[0],gen_pRec[1]);

      // Electron Detection and Fiducial Cuts; common to all settings
      if (!HRS_hallA(ve, pe_central, phie_central))
	continue;
      
      if (fabs(ve.Mag()/pe_central - 1) > mom_cut)
	continue;
      
      if (fabs(ve.Phi() - phie_central) > inplane_cut)
	continue;
      
      if (fabs(ve.Theta() - theta_central) > outplane_cut)
	continue;

      // Lead Proton Detection and Fiducial Cuts; determines detection setting
      if (HRS_hallA(vlead, plead_central_350, philead_central_350))
	{
	  if (fabs(vlead.Mag()/plead_central_350 - 1) < mom_cut
	      and fabs(vlead.Phi() - philead_central_350) < inplane_cut
	      and fabs(vlead.Theta() - theta_central) < outplane_cut)
	    setting = 350;
	}
      else if (HRS_hallA(vlead, plead_central_450, philead_central_450))
	{
	  if (fabs(vlead.Mag()/plead_central_450 - 1) < mom_cut
	      and fabs(vlead.Phi() - philead_central_450) < inplane_cut
	      and fabs(vlead.Theta() - theta_central) < outplane_cut)
	    setting = 450;
	}
      else if (HRS_hallA(vlead, plead_central_550, philead_central_550))
	{
	  if (fabs(vlead.Mag()/plead_central_550 - 1) < mom_cut
	      and fabs(vlead.Phi() - philead_central_550) < inplane_cut
	      and fabs(vlead.Theta() - theta_central) < outplane_cut)
	    setting = 550;
	}
      
      if (setting == 0)
	continue;

      switch(setting)
	{
	case 350:
	  Tp = 0.79;
	  Tpp = Tp*0.84;
	  thetaMmin = 76.*M_PI/180.;
	  break;
	case 450:
	  Tp = 0.77;
	  Tpp = Tp*0.82;
	  thetaMmin = 84.*M_PI/180.;
	  break;
	case 550:
	  Tp = 0.76;
	  Tpp = Tp*0.81;
	  thetaMmin = 88.*M_PI/180.;
	  break;
	default:
	  abort();
	}

      weightp = gen_weight * 1.E33 * Tp;
      weightpp = gen_weight * 1.E33 * eta * Tpp;
      lcweightp = gen_lcweight * 1.E33 * Tp;
      lcweightpp = gen_lcweight * 1.E33 * eta * Tpp;
  
      if (weightp <= 0. and lcweightp <= 0.)
	continue;

      // Recoil Detection and Fiducial Cuts
      if (rec_type == nCode)
	{
	  weightpp = 0.;
	  lcweightpp = 0.;
	}
      else if (!BigBite(vrec,phirec_central))
	{
	  weightpp = 0.;
	  lcweightpp = 0.;
	}
            
      if (fabs(vrec.Phi() - phirec_central) > 5*M_PI/180.)
	{
	  weightpp *= 0;
	  lcweightpp *= 0;
	}
	
      if (fabs(vrec.Theta() - theta_central) > 15*M_PI/180.)
	{
	  weightpp *= 0;
	  lcweightpp *= 0;
	}
      
      if (vrec.Mag() > 0.9 or vrec.Mag() < 0.25)
	{
	  weightpp *= 0;
	  lcweightpp *= 0;
	}

      if (weightpp <= 0. and lcweightpp <= 0. and subcut)
	continue;
      
      // Derived vectors
      TVector3 vbeam(Ebeam,0.,0.);
      TVector3 vq = vbeam - ve;
      TVector3 vmiss=vq - vlead;
      TVector3 vcm=vmiss+vrec;
      TVector3 vrel=0.5*(vmiss-vrec);

      double omega = Ebeam - ve.Mag();
      
      // Cuts
      double Emiss = -m_12C + mN + sqrt( sq(omega + m_12C - sqrt(vlead.Mag2() + sq(mN))) - vmiss.Mag2());

      if (Emiss < 0.032)
	continue;
      
      if (thetaMcut)
	{
	  if (vmiss.Phi() < thetaMmin)
	    continue;
	}
      double QSq = vq.Mag2() - sq(omega);
      double xB = QSq/(2*mN*omega);

      if (xB < 1. and xcut)
	continue;
	
      // Load up tree
      rec_code = rec_type;
      for (int i=0 ; i<3 ; i++)
	{
	  Pe[i] = ve[i];
	  Pp[0][i] = vlead[i];
	  Pp[1][i] = vrec[i];
	}

      // Generated phase space tree for examining generating restrictions
      TVector3 ve_gen(gen_pe[0],gen_pe[1],gen_pe[2]);
      TVector3 vlead_gen(gen_pLead[0],gen_pLead[1],gen_pLead[2]);
      TVector3 vrec_gen(gen_pRec[0],gen_pRec[1],gen_pRec[2]);
      TVector3 vbeam_gen(0.,0.,Ebeam);
      TVector3 vq_gen = vbeam_gen - ve_gen;
      TVector3 vmiss_gen = vlead_gen - vq_gen;
      TVector3 vrel_gen = 0.5*(vmiss_gen - vrec_gen);

      pRel_gen = vrel_gen.Mag();
      thetaRel_gen = vrel_gen.Theta()*180./M_PI;
      phiRel_gen = vrel_gen.Phi()*180./M_PI;
      if (phiRel_gen < 0)
	phiRel_gen += 360.;
      
      thetak_gen = ve_gen.Theta()*180./M_PI;
      phik_gen = ve_gen.Phi()*180./M_PI;
      if (phik_gen < 0)
	phik_gen += 360.;

      outTree->Fill();
    }

  outfile->cd();
  outTree->Write();
  outfile->Close();

  infile->Close();
    
  return 0;
}
