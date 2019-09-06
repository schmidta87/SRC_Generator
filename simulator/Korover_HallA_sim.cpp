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

  double theta_central = 0.5*TMath::Pi();
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
  
  // Output Tree (1p)
  TTree * outTree = new TTree("T","Simulated Data Tree");
  Float_t Q2, Xb, Pe[3], Pe_size, theta_e, phi_e, Pp[2][3], Pp_size[2], pq_angle[2], Ep[2], theta_p[2], phi_p[2], Nu, q[3];
  Float_t Pmiss_q_angle[2], Pmiss_size[2], Pmiss[2][3];
  Float_t z = 0.;
  Double_t weightp, weightpp, weightpn, lcweightp, lcweightpp, lcweightpn;
  Double_t Q2_gen, Xb_gen, phik_gen, phirec_gen, thetarec_gen;
  
  outTree->Branch("Q2",&Q2,"Q2/F");
  outTree->Branch("Xb",&Xb,"Xb/F");
  outTree->Branch("Nu",&Nu,"Nu/F");
  outTree->Branch("q",q,"q[3]/F");
  outTree->Branch("Pe",Pe,"Pe[3]/F");
  outTree->Branch("Pe_size",&Pe_size,"Pe_size/F");
  outTree->Branch("theta_e",&theta_e,"theta_e/F");
  outTree->Branch("phi_e",&phi_e,"phi_e/F");
  outTree->Branch("Pp",Pp,"Pp[1][3]/F");
  outTree->Branch("Pp_size",Pp_size,"Pp_size[1]/F");
  outTree->Branch("pq_angle",pq_angle,"pq_angle[1]/F");
  outTree->Branch("Ep",Ep,"Ep[1]/F");
  outTree->Branch("theta_p",theta_p,"theta_p[1]/F");
  outTree->Branch("phi_p",phi_p,"phi_p[1]/F");
  outTree->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[1]/F");
  outTree->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[1]/F");
  outTree->Branch("Pmiss",Pmiss,"Pmiss[1][3]/F");
  outTree->Branch("weightp",&weightp,"weightp/D");
  outTree->Branch("weightpp",&weightpp,"weightpp/D");
  outTree->Branch("weightpn",&weightpn,"weightpn/D");
  outTree->Branch("lcweightp",&lcweightp,"lcweightp/D");
  outTree->Branch("lcweightpp",&lcweightpp,"lcweightpp/D");
  outTree->Branch("lcweightpn",&lcweightpn,"lcweightpn/D");
  outTree->Branch("Q2_gen",&Q2_gen,"Q2_gen/D");
  outTree->Branch("Xb_gen",&Xb_gen,"Xb_gen/D");
  outTree->Branch("phik_gen",&phik_gen,"phik_gen/D");
  outTree->Branch("phirec_gen",&phirec_gen,"phirec_gen/D");
  outTree->Branch("thetarec_gen",&thetarec_gen,"thetarec_gen/D");

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

      TVector3 ve_gen(gen_pe[0],gen_pe[1],gen_pe[2]);
      TVector3 vrec_gen(gen_pRec[0],gen_pRec[1],gen_pRec[2]);

      // Electron Detection and Fiducial Cuts
      if (!HRS_hallA(ve, pe_central, phie_central))
	continue;
      
      if (abs(ve.Mag()/pe_central - 1) > 4.5e-2)
	continue;
      
      if (abs(ve.Phi() - phie_central) > 30e-3)
	continue;
      
      if (abs(ve.Theta() - theta_central) > 60e-3)
	continue;

      // Lead Proton Detection and Fiducial Cuts
      if (!HRS_hallA(vlead, plead_central, philead_central))
	continue;
      
      if (abs(vlead.Mag()/plead_central - 1) > 4.5e-2)
	continue;
      
      if (abs(vlead.Phi() - philead_central) > 30e-3)
	continue;
	
      if (abs(vlead.Theta() - theta_central) > 60e-3)
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
      
      if (abs(vrec.Phi() - phirec_central) > 4*M_PI/180.)
	continue;
	
      if (abs(vrec.Theta() - theta_central) > 14*M_PI/180.)
	continue;
      
  
      if (weightpp == 0. and weightpn == 0. and lcweightpp == 0. and lcweightpn == 0)
	continue;
      
      // Derived vectors
      TVector3 vbeam = TVector3(Ebeam,0.,0.);
      TVector3 vq = vbeam - ve;
      TVector3 vmiss=vlead-vq;
      TVector3 vcm=vmiss+vrec;
      TVector3 vrel=0.5*(vmiss-vrec);

      double gen_pMiss_Mag = vmiss.Mag();
      double gen_pe_Mag = ve.Mag();
      double gen_Nu = Ebeam - ve.Mag();
      double gen_QSq = vq.Mag2() - sq(gen_Nu);
      double gen_xB = gen_QSq/(2.*mN*gen_Nu);
      double gen_q_Mag = vq.Mag();
      double gen_pLead_Mag = vlead.Mag();
      double gen_pRec_Mag = vrec.Mag();

      // Delta Cut
      double w = sqrt(sq(m_12C + gen_Nu) - vq.Mag2());
      double Lambda = 0.5*(sq(m_11B) - sq(mN) + sq(w));
      double y = ((m_12C + gen_Nu)*sqrt(sq(Lambda) - sq(m_11B*w)) - vq.Mag()*Lambda)/sq(w);

      if (gen_Nu > -0.22*y + (0.95))
	continue;
      
      // Load up tree
      Q2 = gen_QSq;
      Xb = gen_xB;
      Nu = gen_Nu;
      Pe_size = gen_pe_Mag;
      theta_e = ve.Theta() * 180./M_PI;
      phi_e = ve.Phi()*180./M_PI;
      Pp_size[0] = gen_pLead_Mag;
      Pp_size[1] = gen_pRec_Mag;
      pq_angle[0] = vq.Angle(vlead)*180./M_PI;
      pq_angle[1] = vq.Angle(vrec)*180./M_PI;
      Ep[0] = sqrt(gen_pLead_Mag*gen_pLead_Mag + mN*mN);
      Ep[1] = sqrt(gen_pRec_Mag*gen_pRec_Mag + mN*mN);
      theta_p[0] = vlead.Theta()*180./M_PI;
      theta_p[1] = vrec.Theta()*180./M_PI;
      phi_p[0] = vlead.Phi()*180./M_PI;
      phi_p[1] = vrec.Phi()*180./M_PI;
      for (int i=0 ; i<3 ; i++)
	{
	  q[i] = vq[i];
	  Pe[i] = gen_pe[i];
	  Pp[0][i] = vlead[i];
	  Pp[1][i] = vrec[i];
	  Pmiss[0][i] = vlead[i] - vq[i];
	  Pmiss[1][i] = vrec[i] - vq[i];
	}

      Pmiss_q_angle[0] = (vlead - vq).Angle(vq) * 180./M_PI;
      Pmiss_q_angle[1] = (vrec - vq).Angle(vq) * 180./M_PI;
      
      Pmiss_size[0] = (vlead - vq).Mag();
      Pmiss_size[1] = (vrec - vq).Mag();

      // Gen Phase Space
      TVector3 vbeam_gen = TVector3(0.,0.,Ebeam);
      TVector3 vq_gen = vbeam_gen - ve_gen;
      double Nu_gen = Ebeam - ve_gen.Mag();
      Q2_gen = vq_gen.Mag2() - sq(Nu_gen);
      Xb_gen = Q2_gen/(2.*mN*Nu_gen);
      phik_gen = ve_gen.Phi()*180./M_PI;
      if (phik_gen < 0)
	phik_gen += 360.;
      
      phirec_gen = vrec_gen.Phi()*180./M_PI;
      thetarec_gen = vrec_gen.Theta()*180./M_PI;
      
      outTree->Fill();
    }

  infile->Close();
  
  outfile->cd();
  outTree->Write();
  outfile->Close();
  
  return 0;
}
