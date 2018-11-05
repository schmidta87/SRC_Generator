#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"


using namespace std;

double sq(double x){ return x*x; };
int main(int argc, char ** argv){

  if( argc != 7){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "pMissReader /path/to/input/nucleus/file /path/to/output/text/file /path/to/output/root/file [A] [Mass Min] [Mass Max]\n";

    return -1;
    
  }
  
  //Get A trees. Open output file.
  TFile * DataH = new TFile(argv[1]);
  ofstream file;
  file.open(argv[2]);
  TFile * outfile = new TFile(argv[3],"RECREATE");
  double A = atof(argv[4]);
  Double_t nMassMin = atof(argv[5]);
  Double_t nMassMax = atof(argv[6]);
  

  cerr<<"Nucleus file has been opened from: "<<argv[1]<<"\n";
  
  //Number of bins
  Double_t nPBins = 50;
  Double_t nPMin = 0.2;
  Double_t nPMax = 0.6;
  Double_t nPDelta = (nPMax-nPMin)/nPBins;
  Double_t nMassBins = 150;

  //Make trees and histograms for the nuclei
  TTree * TreeH = (TTree*)DataH->Get("T");
  TH2D * hH = new TH2D("Hist2D","A;pMiss_Mag;mMiss;Counts",nPBins,nPMin,nPMax,nMassBins,nMassMin,nMassMax);
  hH ->Sumw2();

  cerr<<"Histograms and Trees successfully created\n";

  
  //Define variables needed for histograms
  Int_t lead_type, rec_type;
  Double_t weight,QSq,xB,nu,pLead_Mag,q_Mag;
  Double_t q[3], pLead[3];
  
  //Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

  //Nucleon Masses
  const double Mp = 0.9383;
  const double Mn = 0.9396;


  //Set int addresses for TreeH
  TreeH->SetBranchAddress("lead_type",&lead_type);
  TreeH->SetBranchAddress("rec_type",&rec_type);

  //Set double addresses for TreeH
  TreeH->SetBranchAddress("weight",&weight);
  TreeH->SetBranchAddress("QSq",&QSq);
  TreeH->SetBranchAddress("xB",&xB);
  TreeH->SetBranchAddress("nu",&nu);
  TreeH->SetBranchAddress("pLead_Mag",&pLead_Mag);
  TreeH->SetBranchAddress("q_Mag",&q_Mag);

  //Set vector addresses for TreeH
  TreeH->SetBranchAddress("q",q);
  TreeH->SetBranchAddress("pLead",pLead);

  cerr << "Finished setting addresses \n";

  int fin = TreeH->GetEntries();

  //Loop over TTree
  for(int i = 0; i < fin; i++){


    //Display completed
    if((i*100)%fin == 0){
    cerr << (i*100)/fin <<"% complete \n";
    }
    
    TreeH->GetEntry(i);

    //Apply first cuts
    if(weight <= 0.) continue;
    if(lead_type == 2112) continue;    
    if(xB < 1.15) continue;
    if(((pLead_Mag/q_Mag) < 0.62) || ((pLead_Mag/q_Mag) > 0.96)) continue;

    Double_t cosThetaplq = ((q[0]*pLead[0]) + (q[1]*pLead[1]) + (q[2]*pLead[2]))/(pLead_Mag * q_Mag);
    Double_t thetaplq = (180/3.14159265) * (TMath::ACos(cosThetaplq));
    if(thetaplq > 25) continue;
    
    //Calculate desired variables
    Double_t Ep = sqrt(sq(Mp) + sq(pLead_Mag));
    Double_t pMiss_Mag = sqrt( sq(q[0]-pLead[0]) + sq(q[1]-pLead[1]) + sq(q[2]-pLead[2]) );
    Double_t mMissSq = sq(nu + (Mp+Mn) - Ep) - sq(pMiss_Mag);

    //Get positive mass squared
    if(mMissSq < 0 ) continue;
    Double_t mMiss = sqrt(mMissSq);
    
    //Apply cut to pMiss given in paper
    if( (pMiss_Mag < nPMin) || (pMiss_Mag > nPMax) ) continue;

    //Fill 2D histogram
    hH->Fill(pMiss_Mag,mMiss,weight);


    
  }
    cerr<<"Finished filling histogram for Nucleus\n";


    Double_t pAct = nPMin - (nPDelta * (1/2) );



    file << "# [Column 1: Missing Momentum] [Column 2: Mean Mass] [Column 3: Mean Mass Error] [Column 4: Mass Sigma] [Column 5: Mass Sigma Error] \n";

    
    for(int j = 0; j < nPBins; j++){ 

      pAct += nPDelta;
      TH1D * hMass = hH->ProjectionY("OnePMiss",j,(j+1));
      cerr<<"The value for P acting and the mean is: "<< pAct <<", "<<(hMass->GetEntries())<<  "\n" <<"";
       if(hMass->GetEntries() == 0) continue;
      
      hMass->Fit("gaus","q","",nMassMin,nMassMax);
      Double_t mean = hMass->GetFunction("gaus")->GetParameter(1);
      Double_t mean_error = hMass->GetFunction("gaus")->GetParError(1);
      Double_t sigma = hMass->GetFunction("gaus")->GetParameter(2);
      Double_t sigma_error = hMass->GetFunction("gaus")->GetParError(2);
      

      file<<pAct<<" "<<mean<<" "<<mean_error<<" "<<sigma<<" "<<sigma_error<<"\n";

      
      hMass->Write();

    }
  


    hH->Write();
    outfile->Close();
  
  file.close();

  cerr<< argv[2]<<" has been completed. \n\n\n";
  
  return 0;
}
