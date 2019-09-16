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
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tKorover_HallA_sim /path/to/gen/file setting /path/to/out/file\n\n";
	return -1;
    }

  TFile * infile = new TFile(argv[1]);
  int setting = atoi(argv[2]);

  double theta_central = 0.5*M_PI;
  double pe_central = 3.602;
  double phie_central = -20.3*M_PI/180.;
  double plead_central, philead_central, phirec_central;
  switch(setting)
    {
    case 1:
      plead_central = 1.38;
      philead_central = 33.5*M_PI/180.;
      phirec_central = 97.*M_PI/180.;
      break;
    case 2:
      plead_central = 1.3;
      philead_central = 29.*M_PI/180.;
      phirec_central = 97.*M_PI/180.;
      break;
    case 3:
      plead_central = 1.19;
      philead_central = 24.5*M_PI/180.;
      phirec_central = 92.*M_PI/180.;
      break;
    case '?':
      return -1;
    default:
      abort();
    }
    
  const double Ebeam = 4.454;
  
  bool verbose = false;
  bool rand_flag = false;
  bool doFCuts = true;

  int c;
  while ((c=getopt (argc-3, &argv[3], "vrC")) != -1)
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
  
  TFile * outfile = new TFile(argv[3],"RECREATE");
  
  // Output Tree
  TTree * outTree = new TTree("T","Simulated Data Tree");
  Float_t Pe[3], Pp[2][3];
  Double_t weightpp, weightpn, lcweightpp, lcweightpn;
  Double_t thetak_gen, phik_gen, thetaRel_gen, phiRel_gen, pRel_gen;
  
  outTree->Branch("Pe",Pe,"Pe[3]/F");
  outTree->Branch("Pp",Pp,"Pp[2][3]/F");
  outTree->Branch("weightpp",&weightpp,"weightpp/D");
  outTree->Branch("weightpn",&weightpn,"weightpn/D");
  outTree->Branch("lcweightpp",&lcweightpp,"lcweightpp/D");
  outTree->Branch("lcweightpn",&lcweightpn,"lcweightpn/D");
  outTree->Branch("thetak_gen",&thetak_gen,"thetak_gen/D");
  outTree->Branch("phik_gen",&phik_gen,"phik_gen/D");
  outTree->Branch("thetaRel_gen",&thetaRel_gen,"thetaRel_gen/D");
  outTree->Branch("phiRel_gen",&phiRel_gen,"phiRel_gen/D");
  outTree->Branch("pRel_gen",&pRel_gen,"pRel_gen/D");

  const int nEvents = inTree->GetEntries();
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %1000000==0 and verbose)
	cout << "Working on event " << event << " out of " << nEvents <<"\n";

      inTree->GetEvent(event);

      // Require a leading proton
      if (lead_type != pCode)
	continue;

      // Create vectors for the particles; rotated for Hall-A geometry
      TVector3 ve(gen_pe[2],gen_pe[0],gen_pe[1]);
      TVector3 vlead(gen_pLead[2],gen_pLead[0],gen_pLead[1]);
      TVector3 vrec(gen_pRec[2],gen_pRec[0],gen_pRec[1]);

      // Electron Detection and Fiducial Cuts
      if (!HRS_hallA(ve, pe_central, phie_central))
	continue;
      
      if (fabs(ve.Mag()/pe_central - 1) > 4.5e-2)
	continue;
      
      if (fabs(ve.Phi() - phie_central) > 27e-3)
	continue;
      
      if (fabs(ve.Theta() - theta_central) > 58e-3)
	continue;

      // Lead Proton Detection and Fiducial Cuts
      if (!HRS_hallA(vlead, plead_central, philead_central))
	continue;
      
      if (fabs(vlead.Mag()/plead_central - 1) > 4.5e-2)
	continue;
      
      if (fabs(vlead.Phi() - philead_central) > 27e-3)
	continue;
	
      if (fabs(vlead.Theta() - theta_central) > 58e-3)
	continue;
      

      if (gen_weight == 0. and gen_lcweight == 0.)
	continue;
      
      weightpp = gen_weight * 1.E33;
      weightpn = gen_weight * 1.E33;
      lcweightpp = gen_lcweight * 1.E33;
      lcweightpn = gen_lcweight * 1.E33;
      
      // Recoil Detection and Fiducial Cuts
      if (rec_type == pCode)
	{
	  weightpn = 0.;
	  lcweightpn = 0.;
	  if (!BigBite(vrec,phirec_central))
	    {
	      weightpp = 0.;
	      lcweightpp = 0.;
	    }
	}
      else
	{
	  weightpp = 0.;
	  lcweightpp = 0.;
	  if (!HAND(vrec,phirec_central))
	    {
	      weightpn = 0.;
	      lcweightpn = 0.;
	    }
	}
            
      if (fabs(vrec.Phi() - phirec_central) > 4*M_PI/180.)
	continue;
	
      if (fabs(vrec.Theta() - theta_central) > 14*M_PI/180.)
	continue;
      
  
      if (weightpp <= 0. and weightpn <= 0. and lcweightpp <= 0. and lcweightpn <= 0)
	continue;
      
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

  infile->Close();
  
  outfile->cd();
  outTree->Write();
  outfile->Close();
  
  return 0;
}
