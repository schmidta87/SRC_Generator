#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"

#include <iostream>
#include <cstdlib>

using namespace std;


int main(int argc, char ** argv){

  if( argc != 5){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "src_reader /path/to/input/nucleus/file /path/to/input/deuterium/file /path/to/output/file [A] \n";
    
  }
  
  //Get Data from the three nuclei and define new file for Histogram
  TFile * DataH = new TFile(argv[1]);
  TFile * DataD = new TFile(argv[2]);
  TFile * histfile = new TFile(argv[3],"RECREATE");
  double A = atof(argv[4]);
  
  //Make trees and histograms for the nuclei
  TTree * TreeH = (TTree*)DataH->Get("T");
  TTree * TreeD = (TTree*)DataD->Get("T");
  TH1D * hH = new TH1D("Helium","Helium;xB;Counts",100,1,2);
  TH1D * hD = new TH1D("Deuterium","Deuterium;xB;Counts",100,1,2);
  hH ->Sumw2();
  hD ->Sumw2();
  
  //Define variables needed for histograms
  Double_t pe[3],weight,QSq,xB;
  Int_t lead_type, rec_type;
  //Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

  //Set addresses for D
  TreeD->SetBranchAddress("weight",&weight);
  TreeD->SetBranchAddress("QSq",&QSq);
  TreeD->SetBranchAddress("xB",&xB);
  TreeD->SetBranchAddress("lead_type",&lead_type);
  TreeD->SetBranchAddress("rec_type",&rec_type);

  //Loop over TTree
  for(int i = 0; i < TreeD->GetEntries(); i++){
    TreeD->GetEntry(i);
    if(QSq < 2) continue;
    hD->Fill(xB,weight);

    
  }

  //Set addresses for H
  TreeH->SetBranchAddress("pe",pe);
  TreeH->SetBranchAddress("weight",&weight);
  TreeH->SetBranchAddress("QSq",&QSq);
  TreeH->SetBranchAddress("xB",&xB);
  TreeH->SetBranchAddress("lead_type",&lead_type);
  TreeH->SetBranchAddress("rec_type",&rec_type);

  //Loop over TTree
  for(int i = 0; i < TreeH->GetEntries(); i++){
    TreeH->GetEntry(i);
    if(QSq < 2) continue;
    if((lead_type == neunum) && (rec_type == neunum)) continue;
    if((lead_type == pronum) && (rec_type == pronum)) continue;
    hH->Fill(xB,weight);

  }


  //Make histogram for ratio and divide it
  TH1D * Hratio = (TH1D*)hH->Clone("ratioHelium");
  Hratio->Divide(hD);
  Hratio->Scale(2/A);
  
  //Fit it to a function and find a2
  TF1 * p0H = new TF1("Hfit","[0]");
  Hratio->Fit("Hfit","r","",1.4,1.7);
  //double a2_H = 2*(p0H->GetParameter(0))/4;
  double a2_H = (p0H->GetParameter(0));
  cout<<"\n"<<"Helium a2: "<<a2_H<<endl;
  
  DataD->Close();
  histfile->cd();
  hH->Write();
  Hratio->Write();

  return 0;
}
