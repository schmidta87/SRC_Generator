
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

using namespace std;

// given by Andrew on Friday April 12th 2019; modified after

int main(int argc, char ** argv){

  if( argc != 3){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "./his_yay [input file = ../DATA/gen_file.root] [output file = ../DATA/his_file.root] \n";

    return -1;
    
  }
  
  //Get D and A trees. Open output file.
  TFile * DataH = new TFile(argv[1]);
  TFile * fo = new TFile(argv[2],"RECREATE");
  
  cerr<<"File has been open from: "<<argv[1]<<"\n";
  
//Make trees and histograms for the nuclei
  TTree * TreeH = (TTree*)DataH->Get("genT");
  TH1D * his_Q_2 = new TH1D("Q_2","Q_2 [GeV^2];Counts",36,1.,5.);
  TH1D * his_X_b =  new TH1D("Xb" ,"X_b [unitless]",40,0.5,1.5);
   

// Sum of squares (error)
  his_Q_2 ->Sumw2();
  his_X_b ->Sumw2();

  cerr<<"Histograms and Trees successfully created\n";

//Define variables needed for histograms
  Double_t Q_2, X_b,  weight;
  
//Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

//Set addresses for D
  TreeH->SetBranchAddress("QSq",&Q_2);
  TreeH->SetBranchAddress("xB",&X_b);
  TreeH->SetBranchAddress("weight",&weight);

    
//Loop over TTree
  for(int i = 0; i < TreeH->GetEntries(); i++){
    TreeH->GetEntry(i);
    his_Q_2->Fill(Q_2, weight);
    his_X_b->Fill(X_b, weight);
  }
    cerr<<"Finished filling histogram\n";


    DataH->Close();
    his_Q_2->Write();
    his_X_b->Write();
    fo->Close();
  cerr<< argv[3]<<" has been completed. \n\n\n";
  
  return 0;
}
