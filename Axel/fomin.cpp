#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"

#include "constants.h"

using namespace std;

const double Ebeam=5.766;
const double deg_margin=1.;

int main(int argc, char ** argv)
{						
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tfomin /path/to/input /path/to/output\n\n";
      return -1;
    }

  TFile * inF = new TFile(argv[1]);
  TTree * inT = (TTree*) inF->Get("T");
  const int nEvents = inT->GetEntries();
  Double_t pe[3], weight;
  inT->SetBranchAddress("pe",pe);
  inT->SetBranchAddress("weight",&weight);
  TVector3 vbeam(0.,0.,Ebeam);

  TFile * outF = new TFile(argv[2],"RECREATE");
  TH2D * h2_kin = new TH2D("kin","Inclusive kinematics;xB;QSq;Counts",20,1.,2.,25,0.,5.);
  h2_kin->Sumw2();
  TH1D * h_xB_18 = new TH1D("xB_18","Inclusive theta=18;xB;Counts",10,1.,2.);
  h_xB_18->Sumw2();
  TH1D * h_xB_22 = new TH1D("xB_22","Inclusive theta=22;xB;Counts",10,1.,2.);
  h_xB_22->Sumw2();
  TH1D * h_xB_26 = new TH1D("xB_26","Inclusive theta=26;xB;Counts",10,1.,2.);
  h_xB_26->Sumw2();

  for (int event=0; event < nEvents ; event++)
    {
      if ((event % 100000)==0)
	cerr << "Working on event " << event << " out of " << nEvents << "\n";

      inT->GetEvent(event);

      // Move the weight into nb to make macroscopic
      weight *= 1.E33;
      if (weight <= 0.)
	continue;

      // Determine kinematics
      TVector3 ve(pe[0],pe[1],pe[2]);
      TVector3 vq = vbeam - ve;
      double omega = Ebeam - ve.Mag();
      double QSq = vq.Mag2() - omega*omega;
      double xB = QSq / (2. * mN * omega);      
      double thetaDeg = ve.Theta() * 180./M_PI;

      h2_kin->Fill(xB,QSq,weight);

      if ( fabs(thetaDeg - 18.) < deg_margin)
	h_xB_18->Fill(xB, weight);

      if ( fabs(thetaDeg - 22.) < deg_margin)
	h_xB_22->Fill(xB, weight);

      if ( fabs(thetaDeg - 26.) < deg_margin)
	h_xB_26->Fill(xB, weight);
    }

  inF->Close();
  
  outF->cd();
  h2_kin->Write();
  h_xB_18->Write();
  h_xB_22->Write();
  h_xB_26->Write();
  outF->Close();

  return 0;
}
