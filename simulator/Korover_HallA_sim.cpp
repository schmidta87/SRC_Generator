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

int main(int argc, char ** argv)
{
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tKorover_HallA_sim /path/to/gen/file /path/to/out/file\n\n";
	return -1;
    }

  TFile * infile = new TFile(argv[1]);

  double theta_central = 0.5*M_PI;
  double pe_central = 3.602;
  double phie_central = -20.3*M_PI/180.;
  double phirec_central;
  
  double plead_central_500 = 1.38;
  double philead_central_500 = 33.5*M_PI/180.;
  double phirec_central_500 = 97.*M_PI/180.;      

  double plead_central_625 = 1.3;
  double philead_central_625 = 29.*M_PI/180.;

  double plead_central_750 = 1.19;
  double philead_central_750 = 24.5*M_PI/180.;
  double phirec_central_750 = 92.*M_PI/180.;

  double eta_pp = 0.73;
  double eta_pn = 0.40;
  double sig_eta_pp = 0.01;
  double sig_eta_pn = 0.014;
  double eta;
  
  double TL;
  double TR;
      
  const double Ebeam = 4.454;
  
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

  if (rand_flag)
    {
      eta_pp += myRand.Gaus(0,sig_eta_pp);
      eta_pn += myRand.Gaus(0,sig_eta_pn);
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
  Double_t weightp, weightpN, lcweightp, lcweightpN;
  Double_t thetak_gen, phik_gen, thetaRel_gen, phiRel_gen, pRel_gen;
  Int_t setting, rec_code;
  
  outTree->Branch("Pe",Pe,"Pe[3]/F");
  outTree->Branch("Pp",Pp,"Pp[2][3]/F");
  outTree->Branch("weightp",&weightp,"weightp/D");
  outTree->Branch("weightpN",&weightpN,"weightpN/D");
  outTree->Branch("lcweightp",&lcweightp,"lcweightp/D");
  outTree->Branch("lcweightpN",&lcweightpN,"lcweightpN/D");
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
      
      if (fabs(ve.Mag()/pe_central - 1) > 4.5e-2)
	continue;
      
      if (fabs(ve.Phi() - phie_central) > 27e-3)
	continue;
      
      if (fabs(ve.Theta() - theta_central) > 58e-3)
	continue;

      // Lead Proton Detection and Fiducial Cuts; determines detection setting
      if (HRS_hallA(vlead, plead_central_500, philead_central_500))
	{
	  if (fabs(vlead.Mag()/plead_central_500 - 1) < 4.5e-2
	      and fabs(vlead.Phi() - philead_central_500) < 27e-3
	      and fabs(vlead.Theta() - theta_central) < 58e-3)
	    setting = 500;
	}
      else if (HRS_hallA(vlead, plead_central_625, philead_central_625))
	{
	  if (fabs(vlead.Mag()/plead_central_625 - 1) < 4.5e-2
	      and fabs(vlead.Phi() - philead_central_625) < 27e-3
	      and fabs(vlead.Theta() - theta_central) < 58e-3)
	    setting = 625;
	}
      else if (HRS_hallA(vlead, plead_central_750, philead_central_750))
	{
	  if (fabs(vlead.Mag()/plead_central_750 - 1) < 4.5e-2
	      and fabs(vlead.Phi() - philead_central_750) < 27e-3
	      and fabs(vlead.Theta() - theta_central) < 58e-3)
	    setting = 750;
	}

      if (setting == 0)
	continue;

      switch(rec_type)
	{
	case pCode:
	  eta = eta_pp;
	  break;
	case nCode:
	  eta = eta_pn;
	  break;
	default:
	  abort();
	}

      switch(setting)
	{
	case 500:
	  TL = 0.75;
	  TR = 0.66;
	case 625:
	  TL = 0.752;
	  TR = 0.7;
	  phirec_central = phirec_central_500;
	  break;
	case 750:
	  TL = 0.756;
	  TR = 0.734;
	  phirec_central = phirec_central_750;
	  break;
	default:
	  abort();
	}

      weightp = gen_weight * 1.E33 * TL;
      weightpN = gen_weight * 1.E33 * eta * TL * TR;
      lcweightp = gen_lcweight * 1.E33 * TL;
      lcweightpN = gen_lcweight * 1.E33 * eta * TL * TR;
  
      if (weightp <= 0. and lcweightp <= 0.)
	continue;

      // Recoil Detection and Fiducial Cuts
      if (rec_type == pCode)
	{
	  if (!BigBite(vrec,phirec_central))
	    {
	      weightpN = 0.;
	      lcweightpN = 0.;
	    }
	}
      else
	{
	  if (!HAND(vrec,phirec_central))
	    {
	      weightpN = 0.;
	      lcweightpN = 0.;
	    }
	}
            
      if (fabs(vrec.Phi() - phirec_central) > 4*M_PI/180.)
	{
	  weightpN *= 0;
	  lcweightpN *= 0;
	}
	
      if (fabs(vrec.Theta() - theta_central) > 14*M_PI/180.)
	{
	  weightpN *= 0;
	  lcweightpN *= 0;
	}
      
      if (vrec.Mag() > 0.9 or vrec.Mag() < 0.3)
	{
	  weightpN *= 0;
	  lcweightpN *= 0;
	}
      
      // Derived vectors
      TVector3 vbeam(Ebeam,0.,0.);
      TVector3 vq = vbeam - ve;
      TVector3 vmiss=vlead-vq;
      TVector3 vcm=vmiss+vrec;
      TVector3 vrel=0.5*(vmiss-vrec);

      double gen_Nu = Ebeam - ve.Mag();

      // Delta Cut
      double W = sqrt(sq(m_4He + gen_Nu) - vq.Mag2());
      double Lambda = 0.5*(sq(m_3H) - sq(mN) + sq(W));
      double y = ((m_4He + gen_Nu)*sqrt(sq(Lambda) - sq(m_3H*W)) - vq.Mag()*Lambda)/sq(W);

      if (gen_Nu > -1.28*y + (0.901))
	continue;

      // (e,e'p) Missing Mass Cut
      double Elead = sqrt(vlead.Mag2() + sq(mN));
      double Erec = sqrt(vrec.Mag2() + sq(mN));
      double m_miss = sqrt(sq(gen_Nu + 2*mN - Elead) - vmiss.Mag2());

      if (m_miss > 1.)
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
