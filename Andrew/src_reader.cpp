#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"

#include <iostream>

using namespace std;


int main(int argc, char ** argv){

  //Get Data from the three nuclei and define new file for Histogram
  TFile * DataH = new TFile("srcProject/He_2md_sig000_E030_P1.root");
  TFile * DataD = new TFile("srcProject/D_P1.root");
  TFile * histfile = new TFile("srcProject/Hist_2md_sig000_E030_p1.root","RECREATE");

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
  Hratio->Scale(0.5);
  
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
