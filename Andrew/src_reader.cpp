#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"

#include "constants.h"
#include "helpers.h"
#include "Nuclear_Info.h"
#include "Cross_Sections.h"


using namespace std;


int main(int argc, char ** argv){

  if( argc != 3){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "src_reader /path/to/input/nucleus/file /path/to/output/file \n";

    return -1;
    
  }
  
  //Get D and A trees. Open output file.
  TFile * DataH = new TFile(argv[1]);
  TFile * fo = new TFile(argv[2],"RECREATE");
  
  cerr<<"Nucleus and Deuterium files have been opened from: "<<argv[1]<<"\n";
  
  //Make trees and histograms for the nuclei
  TTree * TreeH = (TTree*)DataH->Get("genT");
  TH1D * h1p_QSq = new TH1D("ep_QSq","ep;QSq [GeV^2];Counts",36,1.,5.);
  TH1D * h1p_xB =  new TH1D("ep_xB" ,"ep;xB;Counts",40,1.2,2.);
  TH1D * h1p_phiRec = new TH1D("ep_phiRec","ep;phiRec;Counts",36,0,2*M_PI);
  TH1D * h1p_cosThetaRec =  new TH1D("epw_cosThetaRec" ,"ep;cosTheta;Counts",40,-1.0,1.0);
  
  h1p_QSq ->Sumw2();
  h1p_xB ->Sumw2();
  h1p_phiRec ->Sumw2();
  h1p_cosThetaRec ->Sumw2();

  cerr<<"Histograms and Trees successfully created\n";
  
  //Define variables needed for histograms
  Double_t specWeight,weight,QSq,xB,phiRec,cosThetaRec,phi3,pCM[3];
  Int_t lead_type, rec_type;
  
  //Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

  //Set addresses for D
  TreeH->SetBranchAddress("weight",&weight);
  TreeH->SetBranchAddress("specWeight",&specWeight);
  TreeH->SetBranchAddress("QSq",&QSq);
  TreeH->SetBranchAddress("xB",&xB);
  TreeH->SetBranchAddress("phiRec",&phiRec);
  TreeH->SetBranchAddress("cosThetaRec",&cosThetaRec);
  TreeH->SetBranchAddress("phi3",&phi3);
  TreeH->SetBranchAddress("pCM",&pCM);
  TreeH->SetBranchAddress("lead_type",&lead_type);
  TreeH->SetBranchAddress("rec_type",&rec_type);

  //Loop over TTree
  for(int i = 0; i < TreeH->GetEntries(); i++){

    if (i %10000 ==0){
      //  cerr << "specWeight " << specWeight << "\n";
    }
    TreeH->GetEntry(i); 
    h1p_QSq->Fill(QSq,specWeight);
    h1p_xB->Fill(xB,specWeight);
    h1p_phiRec->Fill(phiRec,specWeight);
    h1p_cosThetaRec->Fill(cosThetaRec,specWeight);

  }
    cerr<<"Finished filling histogram for Deuterium\n";


    DataH->Close();
    h1p_QSq->Write();
    h1p_xB->Write();
    h1p_phiRec->Write();
    h1p_cosThetaRec->Write();
    fo->Close();
  cerr<< argv[3]<<" has been completed. \n\n\n";
  
  return 0;
}
