#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector3.h"


using namespace std;

double square(double x){ return x*x; };
int main(int argc, char ** argv){

  if( argc != 6){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "pMissReader /path/to/input/nucleus/file /path/to/input/Electron/Map/file /path/to/input/Proton/Map/file /path/to/output/root/file [sigma pe] [sigma pLead] [sigma pRec] \n";

    return -1;
    
  }

  //Get A trees. Open output file.
  TFile * inData = new TFile(argv[1]);
  TFile * EMData = new TFile(argv[2]);
  TFile * PMData = new TFile(argv[3]);
  cerr<<"The program has begun\n";
  TFile * smearData = new TFile(argv[4],"RECREATE");
  Double_t sigmaPe = atof(argv[5]);
  Double_t sigmaPLead = atof(argv[6]);
  Double_t sigmaPRec = atof(argv[7]);
  

  cerr<<"Nucleus file has been opened from: "<<argv[1]<<"\n";
  

  //Make get tree and Map histos and Random number generator
  TTree * TreeH = (TTree*)inData->Get("T");
  TH3D * EGen = (TH3D*)EMData->Get("Generated Particles");
  TH3D * EAcc = (TH3D*)EMData->Get("Accepted Particles");
  TH3D * PGen = (TH3D*)EMData->Get("Generated Particles");
  TH3D * PAcc = (TH3D*)EMData->Get("Accepted Particles");
  TTree * outtree = new TTree("T","Smear Data Tree");
  TRandom3 myRand(0);  

  cerr<<"Trees successfully created\n";

  //Nucleon Masses
  const double Mp = 0.9383;
  const double Mn = 0.9396;
  Double_t m_Lead;
  Double_t m_Rec;

  //Set up variables for maps
  TVector3 smome;
  TVector3 smomLead;
  Double_t scostE,sphiE,scostLead,sphiLead,genWE,accWE,genWP,accWP;
  Int_t spe_MagBin,scostEBin,sphiEBin,spLead_MagBin,scostLeadBin,sphiLeadBin;
  
  //Set up inTree  
  //Define variables needed for inTree
  Double_t pe[3], q[3], pLead[3], pRec[3], pMiss[3], pCM[3], pRel[3];
  Double_t QSq, xB, nu, pe_Mag, q_Mag, pLead_Mag, pRec_Mag, pMiss_Mag, pCM_Mag, pRel_Mag, theta_pmq, theta_prq, weight;
  Int_t lead_type, rec_type;
  //Set int addresses for inTree
  TreeH->SetBranchAddress("lead_type",&lead_type);
  TreeH->SetBranchAddress("rec_type",&rec_type);
  //Set double addresses for inTree
  TreeH->SetBranchAddress("weight",&weight);
  //TreeH->SetBranchAddress("QSq",&QSq);
  //TreeH->SetBranchAddress("xB",&xB);
  //TreeH->SetBranchAddress("nu",&nu);
  TreeH->SetBranchAddress("pe_Mag",&pe_Mag);
  TreeH->SetBranchAddress("q_Mag",&q_Mag);  
  TreeH->SetBranchAddress("pLead_Mag",&pLead_Mag);
  TreeH->SetBranchAddress("pRec_Mag",&pRec_Mag);
  //TreeH->SetBranchAddress("pMiss_Mag",&pMiss_Mag);
  //TreeH->SetBranchAddress("pCM_Mag",&pCM_Mag);
  //TreeH->SetBranchAddress("pRel_Mag",&pRel_Mag);
  //TreeH->SetBranchAddress("theta_pmq",&theta_pmq);
  //TreeH->SetBranchAddress("theta_prq",&theta_prq);
  //Set vector addresses for inTree
  TreeH->SetBranchAddress("pe",pe);
  TreeH->SetBranchAddress("q",q);
  TreeH->SetBranchAddress("pLead",pLead);
  TreeH->SetBranchAddress("pRec",pRec);
  //TreeH->SetBranchAddress("pMiss",pMiss);
  //TreeH->SetBranchAddress("pCM",pCM);
  //TreeH->SetBranchAddress("pRel",pRel);


  //Set up smearTree
  //Define variables needed for inTree
  Double_t spe[3], sq[3], spLead[3], spRec[3], spMiss[3], spCM[3], spRel[3];
  Double_t sQSq, sxB, snu, spe_Mag, sq_Mag, spLead_Mag, spRec_Mag, spMiss_Mag, spCM_Mag, spRel_Mag, stheta_pmq, stheta_prq, sweight;
  Int_t slead_type, srec_type;
  //Set int addresses for inTree
  outtree->Branch("lead_type",&slead_type,"slead_type/I");
  outtree->Branch("rec_type",&srec_type,"srec_type/I");
  //Set double addresses for inTree
  outtree->Branch("weight",&sweight,"sweight/D");
  outtree->Branch("QSq",&sQSq,"sQSq/D");
  outtree->Branch("xB",&sxB,"sxB/D");
  outtree->Branch("nu",&snu,"snu/D");
  outtree->Branch("pe_Mag",&spe_Mag,"spe_Mag/D");
  outtree->Branch("q_Mag",&sq_Mag,"sq_Mag/D");
  outtree->Branch("pLead_Mag",&spLead_Mag,"spLead_Mag/D");
  outtree->Branch("pRec_Mag",&spRec_Mag,"spRec_Mag/D");
  outtree->Branch("pMiss_Mag",&spMiss_Mag,"spMiss_Mag/D");
  outtree->Branch("pCM_Mag",&spCM_Mag,"spCM_Mag/D");
  outtree->Branch("pRel_Mag",&spRel_Mag,"spRel_Mag/D");
  outtree->Branch("theta_pmq",&stheta_pmq,"stheta_pmq/D");
  outtree->Branch("theta_prq",&stheta_prq,"stheta_prq/D");
  //Set vector addresses for inTree
  outtree->Branch("pe",spe,"spe[3]/D");
  outtree->Branch("q",sq,"sq[3]/D");
  outtree->Branch("pLead",spLead,"spLead[3]/D");
  outtree->Branch("pRec",spRec,"spRec[3]/D");
  outtree->Branch("pMiss",spMiss,"spMiss[3]/D");
  outtree->Branch("pCM",spCM,"spCM[3]/D");
  outtree->Branch("pRel",spRel,"spRel[3]/D");

  cerr << "Finished setting addresses \n";


  int fin = TreeH->GetEntries();
  TreeH->GetEntry(0);
  const Double_t Ebeam[3] = {0.,0.,(q[2]+pe[2])};
  
  //Loop over TTree
  for(int i = 0; i < fin; i++){


    //Display completed
    if(i%100000 == 0){
     cerr << i*100./fin <<"% complete \n";}
    
    TreeH->GetEntry(i);

    //Keep these values the same
    slead_type = lead_type;
    srec_type = rec_type;

    
    if(slead_type == 2212) m_Lead = Mp;
    else m_Lead = Mn;
    if(srec_type == 2212) m_Rec = Mp;
    else m_Rec = Mn;


    //Directly smear the following three values
    spe_Mag = pe_Mag*(1 + myRand.Gaus(0.,sigmaPe));
    spLead_Mag = pLead_Mag*(1 + myRand.Gaus(0.,sigmaPLead));
    spRec_Mag = pRec_Mag*(1 + myRand.Gaus(0.,sigmaPRec));


    //Calculate the 7 momenta
    for(int j = 0; j < 3; j++){
      spe[j] = (spe_Mag/pe_Mag)*pe[j];
      spLead[j] = (spLead_Mag/pLead_Mag) * pLead[j];
      spRec[j] = (spRec_Mag/pRec_Mag) * pRec[j];
      sq[j] = Ebeam[j] - spe[j];
      spMiss[j] = spLead[j] - sq[j];
      spCM[j] = spMiss[j] + spRec[j];
      spRel[j] = 0.5*( spMiss[j] - spRec[j] );      
    };

    //calculate the 4 other momenta magnitudes
    sq_Mag = sqrt(square(sq[0])+square(sq[1])+square(sq[2]));
    spMiss_Mag = sqrt(square(spMiss[0])+square(spMiss[1])+square(spMiss[2]));
    spCM_Mag = sqrt(square(spCM[0])+square(spCM[1])+square(spCM[2]));
    spRel_Mag = sqrt(square(spRel[0])+square(spRel[1])+square(spRel[2]));
    
    //Determine smeared electron values
    snu = Ebeam[2] - spe_Mag;
    sQSq = -(square(snu)-square(sq_Mag));
    sxB = sQSq/(2*m_Lead*snu);

    //Determine smeared angles
    stheta_pmq = acos((spMiss[0]*q[0] + spMiss[1]*q[1] + spMiss[2]*q[2])/spMiss_Mag /sq_Mag);
    stheta_prq = acos((spRec[0]*q[0] + spRec[1]*q[1] + spRec[2]*q[2])/spRec_Mag /sq_Mag);

    //Apply electron maps
    //Get variables for maps
    smome.SetX(spe[0]);
    smome.SetY(spe[1]);
    smome.SetZ(spe[2]);
    scostE = smome.CosTheta();
    sphiE = smome.Phi();
    sphiE = (180/M_PI())*sphiE;
    if(sphiE<(-30)){
      sphiE+=360;
    }

    //Get bins and weights from maps
    spe_MagBin = EGen->GetXaxis()->FindBin(spe_Mag);
    scostEBin = EGen->GetYaxis()->FindBin(scostE);
    sphiEBin = EGen->GetZaxis()->FindBin(sphiE);
    genWE = EGen->GetBinContent(spe_MagBin,scostEBin,sphiEBin);
    accWE =  EAcc->GetBinContent(spe_MagBin,scostEBin,sphiEBin);

    //Apply proton maps
    //Get variables for maps
    smomLead.SetX(spLead[0]);
    smomLead.SetY(spLead[1]);
    smomLead.SetZ(spLead[2]);
    scostLead = smomLead.CosTheta();
    sphiLead = smomLead.Phi();
    sphiLead = (180/M_PI())*sphiLead;
    if(sphiLead<(-30)){
      sphiLead+=360;
    }

    //Get bins and weight from maps
    spLead_MagBin = PGen->GetXaxis()->FindBin(spLead_Mag);
    scostLeadBin = PGen->GetYaxis()->FindBin(scostLead);
    sphiLeadBin = PGen->GetZaxis()->FindBin(sphiLead);
    genWP = PGen->GetBinContent(spe_MagBin,scostEBin,sphiEBin);
    accWP =  PAcc->GetBinContent(spe_MagBin,scostEBin,sphiEBin);
				 
    
    sweight = weight*(accWE/genWE)*(accWP/genWP);
    outtree->Fill();
    
  }
    cerr<<"Finished filling histogram for Nucleus\n";



    outtree->Write();
    smearData->Close();

  cerr<< argv[2]<<" has been completed. \n\n\n";
  
  return 0;
}
