#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "helpers.h"
#include "constants.h"
#include "Nuclear_Info_MF.h"
#include "Cross_Sections.h"

using namespace std;

// Probability windows
const double Qmin=1.;
const double Qmax=5.;
const double Xmin=1.;
const double Xmax=2.;

double gaussian(double x, double sig){return 1/sqrt(2*M_PI*sq(sig)) * exp(sq(x/sig)/2);};

int main(int argc, char ** argv)
{
  if ((argc !=2)&&(argc !=5))
    {
	      cerr << "Wrong number of arguments. Insteady try\n\t"
	   << "gen_weight [A] [Beam energy (GeV)] /path/to/output/file [# of events]\n\n";
      return -1;
    }

  double Ebeam;
  int nEvents;
  int Anum;
  int arn;

  // Read in the arguments
  if (argc == 2){
    std::string Anumstr;
    cout << "Nuclear mass: ";
    getline(std::cin, Anumstr);
    Anum = atoi(Anumstr.c_str());
    arn = 1;
  } else {
    Anum = atoi(argv[1]);
    arn = 3;
  }

  Cross_Sections myCS;
  
  Nuclear_Info_MF myInfo(Anum);
  TFile * outfile = new TFile(argv[arn],"RECREATE");
   
  if (argc == 2){
    std::string Ebeamstr;
    cout << "Beam energy [GeV]: ";
    getline(std::cin, Ebeamstr);
    Ebeam = atof(Ebeamstr.c_str());
    
    std::string nEventsstr;
    cout << "Number of events: ";
    getline(std::cin, nEventsstr);
    nEvents = atoi(nEventsstr.c_str());
    
  } else {
    Ebeam = atof(argv[2]);
    nEvents = atoi(argv[4]);

  }
  
  const TVector3 v1(0.,0.,Ebeam);
    
  // Set up the tree
  TTree * outtree = new TTree("T","Generator Tree");
  Double_t pe[3], q[3], pLead[3], pMiss[3];
  Double_t QSq, xB, nu, pe_Mag, q_Mag, pLead_Mag, pMiss_Mag, theta_pmq, weight, Erem, Estar;
  Int_t nucleon_type;
  outtree->Branch("nucleon_type",&nucleon_type,"nucleon_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("q",q,"q[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pMiss",pMiss,"pMiss[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  outtree->Branch("QSq",&QSq,"QSq/D");
  outtree->Branch("xB",&xB,"xB/D");
  outtree->Branch("nu",&nu,"nu/D");
  outtree->Branch("pe_Mag",&pe_Mag,"pe_Mag/D");
  outtree->Branch("q_Mag",&q_Mag,"q_Mag/D");
  outtree->Branch("pLead_Mag",&pLead_Mag,"pLead_Mag/D");
  outtree->Branch("pMiss_Mag",&pMiss_Mag,"pMiss_Mag/D");
  outtree->Branch("theta_pmq",&theta_pmq,"theta_pmq/D");
  outtree->Branch("Erem",&Erem,"Erem/D");
  outtree->Branch("Estar",&Estar,"Estar/D");

  // Masses and sigma of CM momentum
  TRandom3 myRand(0);
  const double mAm1 = myInfo.get_mAm1();
  const double mA = myInfo.get_mA();
  
  // Loop over events
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %10000 ==0)
	cerr << "Working on event " << event << "\n";

      // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0
      weight =1.;

      // Decide what kind of proton or neutron we are dealing with
      //2212 (proton)
      //2112 (neutron)
      nucleon_type = (myRand.Rndm() > 0.5) ? 2212:2112;
      weight *= 2.;

      // Pick random x, QSq to set up the electron side
      QSq = Qmin + (Qmax-Qmin)*myRand.Rndm();
      xB = Xmin + (Xmax - Xmin)*sqrt(myRand.Rndm());
      nu = QSq/(2.*mN*xB);
      pe_Mag = Ebeam - nu;
      double cosTheta3 = 1. - QSq/(2.*Ebeam*pe_Mag);
      double phi3 = 2.*M_PI*myRand.Rndm();
      TVector3 v3;
      v3.SetMagThetaPhi(pe_Mag,acos(cosTheta3),phi3);
      pe[0]=v3.X();
      pe[1]=v3.Y();
      pe[2]=v3.Z();
      TVector3 vq = v1 - v3;
      q[0]=vq.X();
      q[1]=vq.Y();
      q[2]=vq.Z();
      q_Mag = vq.Mag();

      // Pick random initial motion
      double sig0 = 0.25/sqrt(3.0);
      pMiss[0] = myRand.Gaus(0.,sig0);
      pMiss[1] = myRand.Gaus(0.,sig0);
      pMiss[2] = myRand.Gaus(0.,sig0);
      TVector3 vMiss(pMiss[0],pMiss[1],pMiss[2]);
      TVector3 vLead = vMiss + vq;
      double ELead = sqrt(vLead.Mag2()+sq(mN));
      TVector3 vAm1 = -vMiss;
      double EAm1 = nu + mA - ELead;
      Erem = vAm1.Mag2()/(2*mAm1) + nu - (ELead - mN);
      double Msq = sq(EAm1) - vAm1.Mag2();
      if ((Msq < 0.) or (Erem < 0.)){
	weight=0.;
      }
      else
	{
	  Estar = sqrt(Msq) - mAm1;
	  if (Estar < 0.)
	    {
	      weight = 0.;
	    }
	  else
	    {
	      // Define the other vectors in the tree
	      pLead[0] = pMiss[0] + q[0];
	      pLead[1] = pMiss[1] + q[1];
	      pLead[2] = pMiss[2] + q[2];
	      pLead_Mag = vLead.Mag();
	      pMiss_Mag = vMiss.Mag();
	      theta_pmq = acos((pMiss[0]*q[0] + pMiss[1]*q[1] + pMiss[2]*q[2])/pMiss_Mag /q_Mag);

	      // Calculate the weight
	      weight *= myCS.sigma_eN(Ebeam, v3, vLead, (nucleon_type==2122)) // eN cross section
		* nu/(2.*xB*Ebeam*pe_Mag) // Jacobian for QSq,xB from electron angle and momentum
		* (Qmax-Qmin) * ((sq(Xmax-Xmin))/(2*(xB-Xmin))) // Normalization over electron range
		* 1./(4*sq(M_PI)) // Angular terms
		* 1/(gaussian(pMiss[0], sig0) * gaussian(pMiss[1], sig0) * gaussian(pMiss[2], sig0)) // Normalization over proton range
		* ((nucleon_type==2122) ? myInfo.get_Pp(pMiss_Mag, Erem) : myInfo.get_Pn(pMiss_Mag, Erem) ); // Spectral function sampling.

	    }
	}
      // Fill the tree
      outtree->Fill();      
    } 	  

  // Clean up
  outtree->Write();
  outfile->Close();
  return 0;
}
