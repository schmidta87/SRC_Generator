#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "constants.h"
#include "helpers.h"
#include "Nuclear_Info.h"
#include "Cross_Sections.h"

using namespace std;

void print_help()
{
  cerr << "Usage: ./gen_weight [A] [Beam energy (GeV)] [path/to/output.root] [# of events] [optional flags]\n"
       << "Optional flags:\n"
       << "-h: Help\n"
       << "-S: Silent\n"
       << "-T: Print full output tree\n"
       << "-z: Print zero-weight events\n"
       << "-s <Sigma_CM [GeV]>\n"
       << "-E <E* [GeV]>==<0>\n"
       << "-k <kRel cutoff [GeV]==0.25>\n"
       << "-x <minimum xB>==1\n"
       << "-X <maximum xB>==2\n"
       << "-q <minimum Q2>==1\n"
       << "-Q <minimum Q2>==5\n"
       << "-c <Cross section method>==<cc1>\n"
       << "-f <Form Factor model>==<kelly>\n\n\n";
}

int main(int argc, char ** argv)
{
  if (argc<5)
    {
      print_help();
      return -1;
    }

  // Read in the arguments
  Nuclear_Info myInfo(atoi(argv[1]));
  const double Ebeam=atof(argv[2]);
  const TVector3 v1(0.,0.,Ebeam);
  TFile * outfile = new TFile(argv[3],"RECREATE");
  int nEvents = atoi(argv[4]);

  // Custom settings
  double pRel_cut = 0.25;
  csMethod csMeth=cc1;
  ffModel ffMod=kelly;
  bool quiet=false;
  bool print_full_tree=false;
  bool print_zeros = false;
  // Probability windows
  double Qmin=1.;
  double Qmax=5.;
  double Xmin=1.;
  double Xmax=2.;

  int c;
  while ((c=getopt (argc-4, &argv[4], "hSTzs:E:k:c:f:q:Q:x:X:")) != -1) // First five arguments are not optional flags.
    switch(c)
      {
      case 'h':
        print_help();
        return -1;
      case 'S':
        quiet = true;
        break;
      case 'T':
	print_full_tree=true;
	break;
      case 'z':
	print_zeros=true;
	break;
      case 's':
	myInfo.set_sigmaCM(atof(optarg));
        break;
      case 'E':
	myInfo.set_Estar(atof(optarg));
        break;
      case 'k':
        pRel_cut = atof(optarg);
        break;
      case 'c':
        if (strcmp(optarg,"onshell")==0)
          csMeth=onshell;
        else if (strcmp(optarg,"cc1")==0)
          csMeth=cc1;
        else if (strcmp(optarg,"cc2")==0)
	  csMeth=cc2;
        else if (atoi(optarg)==1)
          csMeth=cc1;
        else if (atoi(optarg)==2)
          csMeth=cc2;
        else
          {
            cerr << "Invalid cross section designation. Allowed values are onshell, cc1 and cc2. Aborting...\n";
            return -1;
          }
        break;
      case 'f':
        if (strcmp(optarg,"kelly")==0)
          ffMod=kelly;
        else if (strcmp(optarg,"dipole")==0)
          ffMod=dipole;
        else{
          cerr << "Invalid form factor model provided. Allowed options are kelly and dipole. Aborting...\n";
          return -1;
	}
	break;
      case 'x':
	Xmin=atof(optarg);
	break;
      case 'X':
	Xmax=atof(optarg);
	break;
      case 'q':
	Qmin=atof(optarg);
	break;
      case 'Q':
	Qmax=atof(optarg);
	break;
      case '?':
	return -1;
      default:
	abort();
  }

  // Adapt cross section to custom arguments
  Cross_Sections myCS(csMeth,ffMod);

  // Set up the tree
  TTree * outtree = new TTree("T","Generator Tree");
  Double_t pe[3], q[3], pLead[3], pRec[3], pMiss[3], pCM[3], pRel[3];
  Double_t  QSq, xB, nu, pe_Mag, q_Mag, pLead_Mag, pRec_Mag, pMiss_Mag, pCM_Mag, pRel_Mag, theta_pmq, theta_prq, weight;
  Int_t lead_type, rec_type;
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  if (print_full_tree)
    {
      outtree->Branch("q",q,"q[3]/D");
      outtree->Branch("pMiss",pMiss,"pMiss[3]/D");
      outtree->Branch("pCM",pCM,"pCM[3]/D");
      outtree->Branch("pRel",pRel,"pRel[3]/D");
      outtree->Branch("QSq",&QSq,"QSq/D");
      outtree->Branch("xB",&xB,"xB/D");
      outtree->Branch("nu",&nu,"nu/D");
      outtree->Branch("pe_Mag",&pe_Mag,"pe_Mag/D");
      outtree->Branch("q_Mag",&q_Mag,"q_Mag/D");
      outtree->Branch("pLead_Mag",&pLead_Mag,"pLead_Mag/D");
      outtree->Branch("pRec_Mag",&pRec_Mag,"pRec_Mag/D");
      outtree->Branch("pMiss_Mag",&pMiss_Mag,"pMiss_Mag/D");
      outtree->Branch("pCM_Mag",&pCM_Mag,"pCM_Mag/D");
      outtree->Branch("pRel_Mag",&pRel_Mag,"pRel_Mag/D");
      outtree->Branch("theta_pmq",&theta_pmq,"theta_pmq/D");
      outtree->Branch("theta_prq",&theta_prq,"theta_prq/D");
    }

  // Masses and sigma of CM momentum
  TRandom3 myRand(0);
  const double mAm2 = myInfo.get_mAm2(); // this includes the effect of Estar
  const double mA = myInfo.get_mA();
  const double sigCM =myInfo.get_sigmaCM();

  // Loop over events
  for (int event=0 ; event < nEvents ; event++)
    {
      if ((event %100000==0) && (!quiet))
	cerr << "Working on event " << event << "\n";

      // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0
      weight =1.;

      // Decide what kind of proton or neutron pair we are dealing with
      lead_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      rec_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      weight *= 4.;

      // Pick random x, QSq to set up the electron side
      QSq = Qmin + (Qmax-Qmin)*myRand.Rndm();
      xB = Xmin + (Xmax - Xmin)*sqrt(myRand.Rndm()); // This draws xB in a triangular distribution, favoring x-> 2
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

      // Pick random CM motion
      pCM[0] = myRand.Gaus(0.,sigCM);
      pCM[1] = myRand.Gaus(0.,sigCM);
      pCM[2] = myRand.Gaus(0.,sigCM);
      TVector3 vCM(pCM[0],pCM[1],pCM[2]);
      pCM_Mag = vCM.Mag();
      TVector3 vAm2 = -vCM;
      double EAm2 = sqrt(vCM.Mag2() + sq(mAm2));
      TVector3 vZ = vCM + vq; // 3 momentum of the pair, useful for calculating kinematics
      double X = mA + nu - EAm2;
      double Z = vZ.Mag();
      double YSq = sq(X) - sq(Z);

      // Pick random recoil angles
      double phiRec = myRand.Rndm()*2.*M_PI;
      double cosThetaRec = 2.*(myRand.Rndm()-0.5); // Maybe in the future figure out min and max angles, and restrict
      double thetaRec = acos(cosThetaRec);

      // Determine other angles
      double thetaCM = vCM.Theta();
      double phiCM = vCM.Phi();
      double cosThetaCMRec = sin(thetaCM)*sin(thetaRec) * (cos(phiCM)*cos(phiRec) + sin(phiCM)*sin(phiRec)) + cos(thetaCM)*cosThetaRec;
      double thetaQ = vq.Theta();
      double phiQ = vq.Phi();
      double cosThetaQRec = sin(thetaQ)*sin(thetaRec) * (cos(phiQ)*cos(phiRec) + sin(phiQ)*sin(phiRec)) + cos(thetaQ)*cosThetaRec;
      double thetaZ = vZ.Theta();
      double phiZ = vZ.Phi();
      double cosThetaZRec = sin(thetaZ)*sin(thetaRec) * (cos(phiZ)*cos(phiRec) + sin(phiZ)*sin(phiRec)) + cos(thetaZ)*cosThetaRec;

      // Determine recoil momentum
      double D = sq(X)*(sq(YSq) + sq(2.*mN)*(vZ.Mag2()*sq(cosThetaZRec) - sq(X)));

      if (D<0)
	{
	  weight=0.;
	}
      else
	{
	  double momRec1 = 0.5*(YSq*Z*cosThetaZRec + sqrt(D))/(sq(X) - vZ.Mag2()*sq(cosThetaZRec));
	  double momRec2 = 0.5*(YSq*Z*cosThetaZRec - sqrt(D))/(sq(X) - vZ.Mag2()*sq(cosThetaZRec));

          bool momRec1Valid = (momRec1>=0.) && (X - sqrt(sq(momRec1) + sq(mN)) >=0.) && (sq(X) - sq(Z) + 2.*Z*momRec1*cosThetaZRec > 0.);
          bool momRec2Valid = (momRec2>=0.) && (X - sqrt(sq(momRec2) + sq(mN)) >=0.) && (sq(X) - sq(Z) + 2.*Z*momRec2*cosThetaZRec >= 0.);

          if ( (!momRec1Valid) && (!momRec2Valid))
            {
              weight=0.;
            }
          else
            {
              if (!momRec1Valid)
                pRec_Mag = momRec2;
              else if (!momRec2Valid)
                pRec_Mag = momRec1;
              else
		{
		  pRec_Mag = (gRandom->Rndm()>0.5)? momRec1 : momRec2;
		  weight*=2.; // because the solution we picked is half as likely
		}

	      // Define the recoil vector
	      TVector3 vRec;
	      vRec.SetMagThetaPhi(pRec_Mag,acos(cosThetaRec),phiRec);
	      pRec[0]=vRec.X();
	      pRec[1]=vRec.Y();
	      pRec[2]=vRec.Z();

	      // Define the other vectors in the tree
	      pMiss[0] = pCM[0] - pRec[0];
	      pMiss[1] = pCM[1] - pRec[1];
	      pMiss[2] = pCM[2] - pRec[2];
	      pLead[0] = pMiss[0] + q[0];
	      pLead[1] = pMiss[1] + q[1];
	      pLead[2] = pMiss[2] + q[2];
	      pRel[0] = 0.5*(pMiss[0] - pRec[0]);
	      pRel[1] = 0.5*(pMiss[1] - pRec[1]);
	      pRel[2] = 0.5*(pMiss[2] - pRec[2]);
	      TVector3 vLead(pLead[0],pLead[1],pLead[2]);
	      TVector3 vRel(pRel[0],pRel[1],pRel[2]);
	      TVector3 vMiss(pMiss[0],pMiss[1],pMiss[2]);
	      pLead_Mag = vLead.Mag();
	      pMiss_Mag = vMiss.Mag();
	      pRel_Mag = vRel.Mag();
	      theta_pmq = acos((pMiss[0]*q[0] + pMiss[1]*q[1] + pMiss[2]*q[2])/pMiss_Mag /q_Mag);
	      theta_prq = acos((pRec[0]*q[0] + pRec[1]*q[1] + pRec[2]*q[2])/pRec_Mag /q_Mag);

	      double Elead = sqrt(sq(mN) + vLead.Mag2());
	      double Erec = sqrt(sq(mN) + vRec.Mag2());

	      // Do a safeguard cut
	      if (pRel_Mag < pRel_cut)
		weight=0.;

	      // Calculate the weight
	      weight *= myCS.sigma_eN(Ebeam, v3, vLead, (lead_type==pCode)) // eN cross section
		* nu/(2.*xB*Ebeam*pe_Mag) // Jacobian for QSq,xB from electron angle and momentum
		* ((sq(Xmax-Xmin))/(2*(xB-Xmin))) // Normalization over range
		* 1./(4.*sq(M_PI)) // Angular terms
		* ((lead_type==rec_type) ? myInfo.get_pp(pRel_Mag) : myInfo.get_pn(pRel_Mag)) // Relative pair probability (from contacts)
		* vRec.Mag2() * Erec * Elead / fabs(Erec*(pRec_Mag - Z*cosThetaZRec) + Elead*pRec_Mag); // Jacobian for delta fnc.
	    }
	}

      // Fill the tree
      if ((weight > 0.) || print_zeros)
	outtree->Fill();

    } 	  

  // Clean up
  outtree->Write();
  outfile->Close();
  return 0;
}
