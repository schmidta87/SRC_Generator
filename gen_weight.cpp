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
  cerr << "Usage: ./gen_weight [Z] [N] [Beam energy (GeV)] [path/to/output.root] [# of events] [optional flags]\n"
       << "Optional flags:\n"
       << "-h: Help\n"
       << "-v: Verbose\n"
       << "-T: Print full output tree\n"
       << "-z: Print zero-weight events\n"
       << "-s <Sigma_CM [GeV]>\n"
       << "-E <E* [GeV]>==<0>\n"
       << "-u <NN interaction>==<AV18>\n"
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
  if (argc<6)
    {
      print_help();
      return -1;
    }
  
    // Read in the arguments
  int Z = atoi(argv[1]);
  int Anum = atoi(argv[2]); 
  const double Ebeam=atof(argv[3]);
  TFile * outfile = new TFile(argv[4],"RECREATE");
  int nEvents = atoi(argv[5]);

  // Custom settings
  bool verbose=false;
  char* u = "AV18";
  bool do_sCM = false;
  double sCM;
  double Estar = 0.;
  double pRel_cut = 0.3;
  csMethod csMeth=cc1;
  ffModel ffMod=kelly;
  bool rand_flag = false;
  bool doRad = true;
  bool print_full_tree=false;
  bool print_zeros = false;
  // Probability windows
  double Qmin=1.;
  double Qmax=5.;
  double Xmin=1.;
  double Xmax=2.;

  double deltaHard(double QSq);
  
  int c;
  while ((c=getopt (argc-4, &argv[4], "hvTzs:E:u:k:c:f:q:Q:x:X:rR")) != -1) // First five arguments are not optional flags.
    switch(c)
      {
      case 'h':
        print_help();
        return -1;
      case 'v':
        verbose = true;
        break;
      case 'T':
	print_full_tree=true;
	break;
      case 'z':
	print_zeros=true;
	break;
      case 's':
	do_sCM = true;
	sCM = atof(optarg);
        break;
      case 'E':
	Estar = atof(optarg);
        break;
      case 'u':
	u = optarg;
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
      case 'r':
	rand_flag = true;
	break;
      case 'R':
	doRad = false;
	break;
      case '?':
	return -1;
      default:
	abort();
  }

  // Initialize Nucleus
  Nuclear_Info myInfo(Z,Anum,u);
  if (do_sCM)
    myInfo.set_sigmaCM(sCM);
  myInfo.set_Estar(Estar);
  
  // Adapt cross section to custom arguments
  Cross_Sections myCS(csMeth,ffMod);

  const TVector3 v1(0.,0.,Ebeam);
  const double lambda_ei = alpha/M_PI * (log( 4.*Ebeam*Ebeam/(me*me)) - 1.);
  
  // Set up the tree
  TTree * outtree = new TTree("genT","Generator Tree");
  Double_t pe[3], q[3], pLead[3], pRec[3], pMiss[3], pCM[3], pRel[3];
  Double_t  QSq, xB, nu, pe_Mag, q_Mag, pLead_Mag, pRec_Mag, pMiss_Mag, pCM_Mag, pRel_Mag, theta_pmq, theta_prq, weight, lcweight;
  Int_t lead_type, rec_type;
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  outtree->Branch("lcweight",&lcweight,"lcweight/D");
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
  const double mA = myInfo.get_mA();
  const double mbar = mA/Anum;
  const double mAmpp = myInfo.get_mAmpp(); // this includes the effect of Estar
  const double mAmpn = myInfo.get_mAmpn();
  const double mAmnn = myInfo.get_mAmnn();
  const double sigCM =myInfo.get_sigmaCM();

  // Loop over events
  for (int event=0 ; event < nEvents ; event++)
    {
      if ((event %100000==0) && (verbose))
	cerr << "Working on event " << event << "\n";

      // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0
      weight =1.;

      // Decide what kind of proton or neutron pair we are dealing with
      lead_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      rec_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      weight *= 4.;

      // Determine mass of A-2 system
      double mAm2;
      if (lead_type == pCode and rec_type == pCode)
	mAm2 = mAmpp;
      else if (lead_type == nCode and rec_type == nCode)
	mAm2 = mAmnn;
      else
	mAm2 = mAmpn;

      // Start with the radiation off the incoming electron
      double DeltaEi = doRad ? pow(myRand.Rndm(),1./lambda_ei) * Ebeam : 0.;
      double Ebeam_eff = Ebeam - DeltaEi;

      // Pick a random xB, QSq prior to FSR
      double xB_eff = Xmin + (Xmax - Xmin)*myRand.Rndm();
      double QSq_eff = Qmin + (Qmax-Qmin)*myRand.Rndm();
      double nu_eff = QSq_eff/(2.*mN*xB_eff);
      double pe_Mag_eff = Ebeam_eff - nu_eff;

      if (pe_Mag_eff < 0.)
	{
	  weight=0.;
	  lcweight=0.;
	}
      else
	{
      
      // The outgoing electron angle won't change in the peaking approximation
      double cosTheta3 = 1. - QSq_eff/(2.*Ebeam_eff*pe_Mag_eff);
      if (fabs(cosTheta3) > 1.)
	{
	  weight=0.;
	  lcweight=0.;
	}
      else
      {
	
      double phi3 = 2.*M_PI*myRand.Rndm();
      double theta3 = acos(cosTheta3);
      TVector3 v3_eff;
      v3_eff.SetMagThetaPhi(pe_Mag_eff, theta3, phi3);
      TVector3 vq_eff = TVector3(0.,0.,Ebeam_eff) - v3_eff;
      TVector3 vqhat_eff = vq_eff.Unit();
	  
      // Sample radiation off the outgoing electron
      double lambda_ef = alpha/M_PI * (log( 4.*pe_Mag_eff*pe_Mag_eff/(me*me)) - 1.);
      double DeltaEf = doRad? pow(myRand.Rndm(),1./lambda_ef) * pe_Mag_eff : 0.;
      pe_Mag = pe_Mag_eff - DeltaEf;

      // This will allow us to calculate apparent quantities
      QSq = 2. * Ebeam * pe_Mag * (1.-cosTheta3);
      nu = Ebeam - pe_Mag;
      xB = QSq/(2.*mN*nu);

      // Fill into vectors
      TVector3 v3;
      v3.SetMagThetaPhi(pe_Mag,theta3,phi3);
      pe[0]=v3.X();
      pe[1]=v3.Y();
      pe[2]=v3.Z();
      TVector3 vq = v1 - v3;
      q[0]=vq.X();
      q[1]=vq.Y();
      q[2]=vq.Z();
      q_Mag = vq.Mag();

      // Pick random CM motion
      TVector3 vCM_eff(myRand.Gaus(0.,sigCM),myRand.Gaus(0.,sigCM),myRand.Gaus(0.,sigCM));
      TVector3 vAm2 = -vCM_eff;
      double EAm2 = sqrt(vCM_eff.Mag2() + sq(mAm2));
      TVector3 vZ = vCM_eff + vq_eff; // 3 momentum of the pair, useful for calculating kinematics
      double X = mA + nu_eff - EAm2;
      double Z = vZ.Mag();
      double YSq = sq(X) - sq(Z);

      // Pick random recoil angles
      double phiRec = myRand.Rndm()*2.*M_PI;
      double cosThetaRec = 2.*(myRand.Rndm()-0.5); // Maybe in the future figure out min and max angles, and restrict
      double thetaRec = acos(cosThetaRec);

      // Determine other angles
      double thetaCM = vCM_eff.Theta();
      double phiCM = vCM_eff.Phi();
      double cosThetaCMRec = sin(thetaCM)*sin(thetaRec) * (cos(phiCM)*cos(phiRec) + sin(phiCM)*sin(phiRec)) + cos(thetaCM)*cosThetaRec;
      double thetaQ = vq_eff.Theta();
      double phiQ = vq_eff.Phi();
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
	      lcweight=0.;
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
		  lcweight*=2.; 
		}

	      // Define the recoil vector
	      TVector3 vRec;
	      vRec.SetMagThetaPhi(pRec_Mag,acos(cosThetaRec),phiRec);
	      pRec[0]=vRec.X();
	      pRec[1]=vRec.Y();
	      pRec[2]=vRec.Z();

	      // Define the other vectors in the tree
	      TVector3 vMiss_eff = vCM_eff - vRec; // This is the true pmiss
	      TVector3 vLead = vMiss_eff + vq_eff; // This is the true plead
	      pLead[0] = vLead.X();
	      pLead[1] = vLead.Y();
	      pLead[2] = vLead.Z();
	      pLead_Mag = vLead.Mag();
	      TVector3 vMiss = vLead - vq; // This is the apparent pmiss
	      pMiss[0] = vMiss.X();
	      pMiss[1] = vMiss.Y();
	      pMiss[2] = vMiss.Z();
	      pMiss_Mag = vMiss.Mag();
	      TVector3 vRel_eff = 0.5*(vMiss_eff - vRec); // This is the true pRel
	      double pRel_eff_Mag = vRel_eff.Mag();

	      double Elead = sqrt(sq(mN) + vLead.Mag2()); // True values
	      double Erec = sqrt(sq(mN) + vRec.Mag2());

	      // Calculate some lightcone quantities
	      double alpha2 = (Erec - vRec.Dot(vqhat_eff))/mbar;
	      double alpha1 = (Elead - nu_eff - vMiss_eff.Dot(vqhat_eff))/mbar;
	      double alphaCM = alpha1 + alpha2;
	      double alpharel = 2*alpha2/alphaCM;
	      double alphaAm2 = Anum - alphaCM;

	      // Perpendicular components are also EFFECTIVE quantities
	      TVector3 vMiss_eff_perp = vMiss_eff - vMiss_eff.Dot(vqhat_eff)*vqhat_eff;
	      TVector3 vRec_perp = vRec - vRec.Dot(vqhat_eff)*vqhat_eff;
	      TVector3 k_perp = alpha1/alphaCM*vRec_perp - alpha2/alphaCM*vMiss_eff_perp;
	      double kSq = (sq(mN) + k_perp.Mag2())/(alpharel*(2.-alpharel)) - sq(mN);
	      double k = sqrt(kSq);

	      
	      // Do a safeguard cut
	      if (pRel_eff_Mag < pRel_cut)
		  weight=0.;
	      if (k < pRel_cut)
		  lcweight=0.;

	      pRel[0] = 0.5*(pMiss[0] - pRec[0]); // This is the apparent pRel;
	      pRel[1] = 0.5*(pMiss[1] - pRec[1]);
	      pRel[2] = 0.5*(pMiss[2] - pRec[2]);
	      TVector3 vRel(pRel[0],pRel[1],pRel[2]);
	      pRel_Mag = vRel.Mag();

	      pCM[0] = pMiss[0] + pRec[0]; // Apparent pCM
	      pCM[1] = pMiss[1] + pRec[1];
	      pCM[2] = pMiss[2] + pRec[2];
	      TVector3 vCM(pCM[0],pCM[1],pCM[2]);
	      pCM_Mag = vCM.Mag();

	      // These are apparent angles
	      theta_pmq = acos((pMiss[0]*q[0] + pMiss[1]*q[1] + pMiss[2]*q[2])/pMiss_Mag /q_Mag);
	      theta_prq = acos((pRec[0]*q[0] + pRec[1]*q[1] + pRec[2]*q[2])/pRec_Mag /q_Mag);

	      // Calculate the weight
	      weight *= myCS.sigma_eN(Ebeam, v3, vLead, (lead_type==pCode)) // eN cross section
		* nu/(2.*xB*Ebeam*pe_Mag) // Jacobian for QSq,xB from electron angle and momentum
		* ((sq(Xmax-Xmin))/(2*(xB-Xmin))) // Normalization over range
		* 1./(4.*sq(M_PI)) // Angular terms
		* myInfo.get_S(pRel_Mag,lead_type,rec_type)// Relative pair probability (from contacts)
		* vRec.Mag2() * Erec * Elead / fabs(Erec*(pRec_Mag - Z*cosThetaZRec) + Elead*pRec_Mag); // Jacobian for delta fnc.

	      if (kSq < 0)
		lcweight = 0;
	      else
		{
		  // Calculate the lightcone weight
		  lcweight *= myCS.sigma_eN(Ebeam_eff, v3_eff, vLead, (lead_type==pCode))/alpha1 // eN cross section
		    * nu_eff/(2.*xB_eff*Ebeam_eff*pe_Mag_eff) * (Qmax-Qmin) * (Xmax-Xmin) // Jacobian for QSq,xB
		    * (doRad ? (1. - deltaHard(QSq_eff)) * pow(Ebeam/sqrt(Ebeam*pe_Mag),lambda_ei) * pow(pe_Mag_eff/sqrt(Ebeam*pe_Mag),lambda_ef) : 1.) // Radiative weights
		    * 1./(4.*sq(M_PI)) // Angular terms
		    * sqrt(mN*mN + kSq)/Erec * 1./(2.-alpharel) * ((lead_type==rec_type) ? myInfo.get_pp(k) : myInfo.get_pn(k)) // Contacts
		    * vRec.Mag2() * Erec * Elead / fabs(Erec*(pRec_Mag - Z*cosThetaZRec) + Elead*pRec_Mag) // Jacobian for delta fnc.
		    * mbar*((Anum>2)?(alphaAm2/EAm2 * exp((sq(vCM_eff.Dot(vqhat_eff))-sq(mbar*(2.-alphaCM)))/(2.*sq(sigCM)))):1.); // Change in center-of-mass motion in lightcone picture
		}
	    }
	}
	}
	}

      // Fill the tree
      if ((weight > 0.) || (lcweight > 0.) || print_zeros)
	outtree->Fill();

    } 	  

  // Clean up
  outtree->Write();
  outfile->Close();
  return 0;
}

double deltaHard(double QSq)
{
  return 2.*alpha/M_PI * ( -13./12.*log(QSq/(me*me)) + 8./3.);
}
