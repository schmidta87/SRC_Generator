#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "helpers.h"
#include "Nuclear_Info.h"
#include "Cross_Sections.h"

using namespace std;

void do_SCX(int $lead_type, int &rec_type, double r);

void print_help()
{
  cerr << "Usage: ./gen_weight [Z] [N] [Beam energy (GeV)] [path/to/output.root] [# of events] [optional flags]\n"
       << "Optional flags:\n"
       << "-h: Help\n"
       << "-v: Verbose\n"
       << "-T: Print full output tree\n"
       << "-z: Print zero-weight events\n"
       << "-P: Use text file to generate phase space\n"
       << "-O: Turn off radiative effects\n"
       << "-o: Turn on Single-Charge Exchange\n"
       << "-s <Sigma_CM [GeV]>\n"
       << "-E <E* [GeV]>==<0>\n"
       << "-u <NN interaction>==<AV18>\n"
       << "-k <kRel cutoff [GeV]==0.25>\n"
       << "-c <Cross section method>==<cc1>\n"
       << "-f <Form Factor model>==<kelly>\n"
       << "-r: Randomize uncertain parameters\n"
       << "-R: Use same randomized parameters as given file\n\n\n";
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
  int N = atoi(argv[2]);
  int Anum = Z + N;
  const double Ebeam=atof(argv[3]);
  TFile * outfile = new TFile(argv[4],"RECREATE");
  int nEvents = atoi(argv[5]);

  // Custom settings
  bool verbose=false;
  char* u = "AV18";
  bool do_sCM = false;
  double sCM;
  bool do_Estar = false;
  double Estar;
  double pRel_cut = 0.30;
  double pRel_range = 0.10;
  csMethod csMeth=cc1;
  ffModel ffMod=kelly;
  bool rand_flag = false;
  bool doRad = true;
  bool doSCX = false;
  bool print_full_tree=false;
  bool print_zeros = false;
  bool use_input_params = false;
  TFile * param_file;
  TVectorT<double> params(19);
  // Probability windows
  bool custom_ps = false;
  char* phase_space;
  double phi3min=0.;
  double phi3max=2*M_PI;
  double theta3min=0.;
  double theta3max=M_PI/2.;
  double pRelmin = 0.25;
  double pRelmax = 1.05;
  double phiRelmin=0.;
  double phiRelmax=2*M_PI;
  double thetaRelmin=0.;
  double thetaRelmax=M_PI;  
  
  double deltaHard(double QSq);
  
  int c;
  while ((c=getopt (argc-5, &argv[5], "hvTzs:E:u:k:c:f:rP:OoR:")) != -1) // First five arguments are not optional flags.
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
	do_Estar = true;
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
      case 'r':
	rand_flag = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
	break;
      case 'O':
	doRad = false;
	break;
      case 'o':
	doSCX = true;
	break;
      case 'R':
	use_input_params = true;
	param_file = new TFile(optarg,"READ");
	break;
      case '?':
	return -1;
      default:
	abort();
  }

  if (custom_ps)
    {
      ifstream ps_file(phase_space);
      string param;
      double low, high;
      while (ps_file >> param >> low >> high)
	{
	  if (param == "phi3" or param == "phik")
	    {
	      phi3min = low;
	      phi3max = high;
	    }
	  else if (param == "phi3_deg" or param == "phik_deg")
	    {
	      phi3min = low*M_PI/180.;
	      phi3max = high*M_PI/180.;
	    }
	  else if (param == "theta3" or param == "thetak")
	    {
	      theta3min = low;
	      theta3max = high;
	    }
	  else if (param == "theta3_deg" or param == "thetak_deg")
	    {
	      theta3min = low*M_PI/180.;
	      theta3max = high*M_PI/180.;
	    }
	  else if (param == "pRel" or param == "prel")
	    {
	      pRelmin = low;
	      pRelmax = high;
	    }
	  else if (param == "phiRel" or param == "phik")
	    {
	      phi3min = low;
	      phi3max = high;
	    }
	  else if (param == "phiRel_deg" or param == "phirel_deg")
	    {
	      phiRelmin = low*M_PI/180.;
	      phiRelmax = high*M_PI/180.;
	    }
	  else if (param == "thetaRel" or param == "thetarel")
	    {
	      thetaRelmin = low;
	      thetaRelmax = high;
	    }
	  else if (param == "thetaRel_deg" or param == "thetarel_deg")
	    {
	      thetaRelmin = low*M_PI/180.;
	      thetaRelmax = high*M_PI/180.;
	    }
	  else
	    {
	      cerr << "Invalid phase space parameter provided. Aborting...\n";
	      return -1;
	    }
	}
    }

  double cosTheta3min = cos(theta3max);
  double cosTheta3max = cos(theta3min);
  double cosThetaRelmin = cos(thetaRelmax);
  double cosThetaRelmax = cos(thetaRelmin);

  // Initialize Nucleus
  Nuclear_Info myInfo(Z,N,u);
  if (do_sCM)
    myInfo.set_sigmaCM(sCM);
  if (do_Estar)
    myInfo.set_Estar(Estar);

  // Randomization
  TRandom3 myRand(0);
  if (rand_flag)
    {
      myInfo.randomize();
      pRel_cut = pRel_cut + (myRand.Uniform() - 0.5)*pRel_range;
    }
  else if (use_input_params)
    {
      params = *(TVectorT<double>*)param_file->Get("parameters");
      param_file->Close();
      myInfo.set_sigmaCM(params[0]);
      myInfo.set_Cpp0(params[1]);
      myInfo.set_Cpn0(params[2]);
      myInfo.set_Cnn0(params[3]);
      myInfo.set_Cpn1(params[4]);
      myInfo.set_Estar(params[5]);
      
      std::vector<double> Pscx = myInfo.get_SCX_Ps();
      myInfo.set_SCX_Ps(Pscx[0],
			Pscx[1],
			Pscx[2],
			Pscx[3],
			Pscx[4],
			Pscx[5],
			Pscx[6],
			Pscx[7],
			Pscx[8],
			Pscx[9],
			Pscx[10],
			Pscx[11]);
    }
  
  // Adapt cross section to custom arguments
  Cross_Sections myCS(csMeth,ffMod);

  const TVector3 v1(0.,0.,Ebeam);
  const double lambda_ei = alpha/M_PI * (log( 4.*Ebeam*Ebeam/(me*me)) - 1.);
  
  // Set up the tree
  outfile->cd();
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
  outtree->Branch("QSq",&QSq,"QSq/D");
  outtree->Branch("xB",&xB,"xB/D");
  if (print_full_tree)
    {
      outtree->Branch("q",q,"q[3]/D");
      outtree->Branch("pMiss",pMiss,"pMiss[3]/D");
      outtree->Branch("pCM",pCM,"pCM[3]/D");
      outtree->Branch("pRel",pRel,"pRel[3]/D");
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
  const double mA = myInfo.get_mA();
  const double mbar = mA/Anum;
  const double mAmpp = myInfo.get_mAmpp(); // this includes the effect of Estar
  const double mAmpn = myInfo.get_mAmpn();
  const double mAmnn = myInfo.get_mAmnn();
  const double sigCM = myInfo.get_sigmaCM();

  // Prepare vector of parameters to be output
  params[0] = sigCM;
  params[1] = myInfo.get_Cpp0();
  params[2] = myInfo.get_Cpn0();
  params[3] = myInfo.get_Cnn0();
  params[4] = myInfo.get_Cpn1();
  params[5] = myInfo.get_Estar();
  std::vector<double> Ps = myInfo.get_SCX_Ps();
  params[6] = Ps[0];
  params[7] = Ps[1];
  params[8] = Ps[2];
  params[9] = Ps[3];
  params[10] = Ps[4];
  params[11] = Ps[5];
  params[12] = Ps[6];
  params[13] = Ps[7];
  params[14] = Ps[8];
  params[15] = Ps[9];
  params[16] = Ps[10];
  params[17] = Ps[11];
  params[18] = pRel_cut;
  
  // Loop over events
  for (int event=0 ; event < nEvents ; event++)
    {
      if ((event %100000==0) && (verbose))
	cout << "Working on event " << event << "\n";

      // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0
      weight = 1.;
      lcweight = 1.;

      // Decide what kind of proton or neutron pair we are dealing with
      lead_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      rec_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      weight *= 4.;
      lcweight *= 4.;

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

      // Pick random electron angles
      double phi3 = phi3min + (phi3max-phi3min)*myRand.Rndm();
      double cosTheta3 = cosTheta3min + (cosTheta3max - cosTheta3min)*myRand.Rndm();
      double theta3 = acos(cosTheta3);
      TVector3 v3hat_eff;
      v3hat_eff.SetMagThetaPhi(1,theta3,phi3);
      
      // Pick random CM motion
      TVector3 vCM_eff(myRand.Gaus(0.,sigCM),myRand.Gaus(0.,sigCM),myRand.Gaus(0.,sigCM));
      TVector3 vAm2 = -vCM_eff;
      double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));

      // Pick random relative motion
      double phiRel = phiRelmin + (phiRelmax-phiRelmin)*myRand.Rndm();
      double cosThetaRel = cosThetaRelmin + (cosThetaRelmax - cosThetaRelmin)*myRand.Rndm();
      double thetaRel = acos(cosThetaRel);
      double pRel_Mag_eff = pRelmin + (pRelmax - pRelmin)*myRand.Rndm();
      TVector3 vRel_eff;
      vRel_eff.SetMagThetaPhi(pRel_Mag_eff,thetaRel,phiRel);

      // Determine initial nucleon momenta
      TVector3 vMiss_eff = 0.5*vCM_eff + vRel_eff;
      TVector3 vRec = 0.5*vCM_eff - vRel_eff;
      pRec[0]=vRec.X();
      pRec[1]=vRec.Y();
      pRec[2]=vRec.Z();
      double Erec = sqrt(sq(mN) + vRec.Mag2());	      
      
      // Calculate scattered electron energy
      TVector3 veN = TVector3(0.,0.,Ebeam_eff) + vMiss_eff;
      double pe_Mag_eff = 0.5*(sq(mA + Ebeam_eff - Erec - EAm2) - veN.Mag2() - sq(mN))/(mA + Ebeam_eff - Erec - EAm2 - veN.Dot(v3hat_eff));

      if ((pe_Mag_eff < 0.) or (pe_Mag_eff > Ebeam))
	{
	  weight=0.;
	  lcweight=0.;
	}
      else
	{
	  TVector3 v3_eff;
	  v3_eff.SetMagThetaPhi(pe_Mag_eff,theta3,phi3);
	  TVector3 vq_eff = TVector3(0.,0.,Ebeam_eff) - v3_eff;
	  TVector3 vqhat_eff = vq_eff.Unit();
	  double QSq_eff = 2.*Ebeam_eff * pe_Mag_eff * (1 - cosTheta3); 
	  double nu_eff = Ebeam_eff - pe_Mag_eff;
	  
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

	  // Determine other angles
	  //double thetaCM = vCM_eff.Theta();
	  //double phiCM = vCM_eff.Phi();
	  //double cosThetaCMRec = sin(thetaCM)*sin(thetaRec) * (cos(phiCM)*cos(phiRec) + sin(phiCM)*sin(phiRec)) + cos(thetaCM)*cosThetaRec;
	  //double thetaQ = vq_eff.Theta();
	  //double phiQ = vq_eff.Phi();
	  //double cosThetaQRec = sin(thetaQ)*sin(thetaRec) * (cos(phiQ)*cos(phiRec) + sin(phiQ)*sin(phiRec)) + cos(thetaQ)*cosThetaRec;

	  // Define the other vectors in the tree
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
	  
	  double Elead = sqrt(sq(mN) + vLead.Mag2()); // True values

	  // Calculate some lightcone quantities
	  double alpha2 = (Erec - vRec.Dot(vqhat_eff))/mbar;
	  double alpha1 = (Elead - nu_eff - vMiss_eff.Dot(vqhat_eff))/mbar;
	  double alphaCM = alpha1 + alpha2;
	  double alpharel = 2*alpha2/alphaCM;
	  double alphaAm2 = Anum - alphaCM;
	  
	  TVector3 vMiss_eff_perp = vMiss_eff - vMiss_eff.Dot(vqhat_eff)*vqhat_eff;
	  TVector3 vRec_perp = vRec - vRec.Dot(vqhat_eff)*vqhat_eff;
	  TVector3 k_perp = alpha1/alphaCM*vRec_perp - alpha2/alphaCM*vMiss_eff_perp;
	  double kSq = (sq(mN) + k_perp.Mag2())/(alpharel*(2.-alpharel)) - sq(mN);
	  double k = sqrt(kSq);

	  // Do a safeguard cut
	  if (pRel_Mag_eff < pRel_cut)
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
	  weight *= myCS.sigma_eN(Ebeam_eff, v3_eff, vLead, (lead_type==pCode)) // eN cross section
	    * (doRad ? (1. - deltaHard(QSq_eff)) * pow(Ebeam/sqrt(Ebeam*pe_Mag),lambda_ei) * pow(pe_Mag_eff/sqrt(Ebeam*pe_Mag),lambda_ef) : 1.) // Radiative weights
	    * 1./(2.*sq(M_PI)) * (phi3max - phi3min)/(2*M_PI) * (cosTheta3max - cosTheta3min)/(2) * (phiRelmax - phiRelmin)/(2*M_PI) * (cosThetaRelmax - cosThetaRelmin)/(2) // Angular terms
	    * myInfo.get_S(pRel_Mag_eff,lead_type,rec_type) * (pRelmax - pRelmin) // Contacts
	    * sq(pRel_Mag_eff) / fabs(1 - vLead.Dot(v3hat_eff)/Elead); // Jacobian for delta fnc.
	    
	  if (kSq < 0)
	    lcweight = 0;
	  else
	    {
	      // Calculate the lightcone weight
	      lcweight *= myCS.sigma_eN(Ebeam_eff, v3_eff, vLead, (lead_type==pCode))/alpha1 // eN cross section
		* (doRad ? (1. - deltaHard(QSq_eff)) * pow(Ebeam/sqrt(Ebeam*pe_Mag),lambda_ei) * pow(pe_Mag_eff/sqrt(Ebeam*pe_Mag),lambda_ef) : 1.) // Radiative weights
		* 1./(2.*sq(M_PI)) * (phi3max - phi3min)/(2*M_PI) * (cosTheta3max - cosTheta3min)/(2) * (phiRelmax - phiRelmin)/(2*M_PI) * (cosThetaRelmax - cosThetaRelmin)/(2) // Angular terms
		* sqrt(mN*mN + kSq)/Erec * 1./(2.-alpharel) * myInfo.get_S(k,lead_type,rec_type) * (pRelmax - pRelmin) // Contacts
		* sq(pRel_Mag_eff) / fabs(1 - vLead.Dot(v3hat_eff)/Elead) // Jacobian for delta fnc.
		* mbar * ((Anum>2)?(alphaAm2/EAm2 * exp((sq(vCM_eff.Dot(vqhat_eff))-sq(mbar*(2.-alphaCM)))/(2.*sq(sigCM)))):1.); // Change in center-of-mass motion in lightcone picture
	    }
	  
	}

      // Here is where we do single charge exchange
      if (doSCX)
	{
	  myInfo.do_SCX(lead_type, rec_type, gRandom->Rndm());
	}
      
      // Fill the tree
      if ((weight > 0.) || (lcweight > 0.) || print_zeros)
	outtree->Fill();

    } 	  

  // Clean up
  outtree->Write();
  params.Write("parameters");
  outfile->Close();
  return 0;
}

double deltaHard(double QSq)
{
  return 2.*alpha/M_PI * ( -13./12.*log(QSq/(me*me)) + 8./3.);
}
