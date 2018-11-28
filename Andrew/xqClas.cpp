#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cstring> 

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"


using namespace std;

double sq(double x){ return x*x; };
int main(int argc, char ** argv){

  if( argc != 4){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "xqClas /path/to/input/data/file /path/to/output/text/file /path/to/output/root/file \n";

    return -1;
    
  }
  
  //Get A trees. Open output file.
  TFile * DataH = new TFile(argv[1]);
  ofstream file;
  file.open(argv[2]);
  TFile * outfile = new TFile(argv[3],"RECREATE");

  //Bins
  Double_t nPBins = 50;
  Double_t nPMin = 0.05;
  Double_t nPMax = 0.85;
  Double_t nPDelta = (nPMax-nPMin)/nPBins;
  Double_t nXBins = 50;
  Double_t nXMin = 1;
  Double_t nXMax = 2;
  Double_t nQBins = 50;
  Double_t nQMin = 0;
  Double_t nQMax = 6;
  
  
  //Make trees and histograms for data
  TTree * TreeH = (TTree*)DataH->Get("T");
  TH2D * hX = new TH2D("x_B","A;pMiss_Mag;x_B;Counts",nPBins,nPMin,nPMax,nXBins,nXMin,nXMax);
  TH2D * hQ = new TH2D("Q^2","A;pMiss_Mag;Q^2;Counts",nPBins,nPMin,nPMax,nQBins,nQMin,nQMax);

  //Define varialbes needed
  Int_t npar;
  Double_t Qsq,xB,nu;
  Double_t px[9],py[9],pz[9];
  const Double_t Ebeam = 4.461;
  const double Mp = 0.9383;
  const double Mn = 0.9396;
  


  cerr<<"Nucleus file has been opened from: "<<argv[1]<<"\n";


  //Set addresses for TTree
  TreeH->SetBranchAddress("nParticles",&npar);
  TreeH->SetBranchAddress("Q2",&Qsq);
  TreeH->SetBranchAddress("Xb",&xB);
  TreeH->SetBranchAddress("Nu",&nu);
  TreeH->SetBranchAddress("mom_x",px);
  TreeH->SetBranchAddress("mom_y",py);
  TreeH->SetBranchAddress("mom_z",pz);

  cerr<<"Nucleus file has been opened from: "<<argv[1]<<"\n";

 int fin = TreeH->GetEntries();
 
  //Loop over TTree
  for(int i = 0; i < fin; i++){

    //Display completed
    if((i%100000) == 0){
    cerr << (i*100.)/fin <<"% complete \n";
    file.flush();
    }
    
    TreeH->GetEntry(i);
    
    //get q
    Double_t q[3] = {-px[0],-py[0],Ebeam-pz[0]};
    Double_t q_Mag = sqrt(sq(q[0])+sq(q[1])+sq(q[2]));
    Double_t pMiss_Mag = 7.6;

    int j = 1;
    bool pass = false;

        while((j<npar)&&(!pass)){

	  Double_t pLead[3] = {px[j],py[j],pz[j]};
	  j++;
	  Double_t pLead_Mag = sqrt(sq(pLead[0])+sq(pLead[1])+sq(pLead[2]));
      
	  //Apply first cuts   
	  if(xB < 1.15) continue;
	  if(((pLead_Mag/q_Mag) < 0.62) || ((pLead_Mag/q_Mag) > 0.96)) continue;

	  Double_t cosThetaplq = ((q[0]*pLead[0]) + (q[1]*pLead[1]) + (q[2]*pLead[2]))/(pLead_Mag * q_Mag);
	  Double_t thetaplq = (180/3.14159265) * (TMath::ACos(cosThetaplq));
	  if(thetaplq > 25) continue;

	  //Calculate desired variables
	  Double_t Ep = sqrt(sq(Mp) + sq(pLead_Mag));
	  pMiss_Mag = sqrt( sq(q[0]-pLead[0]) + sq(q[1]-pLead[1]) + sq(q[2]-pLead[2]) );

	  //Apply cut to pMiss given in paper
	  if( (pMiss_Mag < nPMin) || (pMiss_Mag > nPMax) ) continue;

	  pass = true;
      
    }
   
 
    if(pass){
      //Fill 2D histogram
      hX->Fill(pMiss_Mag,xB);
      hQ->Fill(pMiss_Mag,Qsq);
      file<<pMiss_Mag<<" "<<xB<<" "<<Qsq<<endl;
    }


    memset(px, 0, sizeof(px));
    memset(py, 0, sizeof(py));
    memset(pz, 0, sizeof(pz));
        
    
    
  }
    cerr<<"Finished filling histogram for Nucleus\n";

  
  hX->Write();
  hQ->Write();
  outfile->Close();
  file.close();


  
  return 0;
}
