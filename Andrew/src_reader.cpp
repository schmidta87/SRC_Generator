#include <iostream>
#include <fstream>
#include <cstdlib>


#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"


using namespace std;


int main(int argc, char ** argv){

  if( argc != 5){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "src_reader /path/to/input/nucleus/file /path/to/input/deuterium/file /path/to/output/file [A] \n";

    return -1;
    
  }
  
  //Get Data from the three nuclei and define new file for Histogram
  TFile * DataH = new TFile(argv[1]);
  TFile * DataD = new TFile(argv[2]);
  //TFile * histfile = new TFile(argv[3],"RECREATE");
  ofstream file;
  file.open(argv[3]);
  double A = atof(argv[4]);

  cerr<<"Nucleus and Deuterium files have been opened at:"<<DataH<<" "<<DataD<<"\n";
  
  //Make trees and histograms for the nuclei
  //TTree * TreeH = NULL;
  //TTree * TreeD = NULL;
  TTree * TreeH = (TTree*)DataH->Get("T");
  TTree * TreeD = (TTree*)DataD->Get("T");
  TH1D * hH = new TH1D("Helium","Helium;xB;Counts",40,1,2);
  TH1D * hD = new TH1D("Deuterium","Deuterium;xB;Counts",40,1,2);
  hH ->Sumw2();
  hD ->Sumw2();

  cerr<<"Histograms and Trees successfully created\n";
  
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

    cerr<<"Finished filling histogram for Deuterium\n";

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
    hH->Fill(xB,weight);

  }
    cerr<<"Finished filling histogram for Nucleus\n";

  //Make histogram for ratio and divide it
  TH1D * Hratio = (TH1D*)hH->Clone("ratioHelium");
  Hratio->Divide(hD);
  Hratio->Scale(2/A);
  
  //Fit it to a function and find a2
  //TF1 * p0H = new TF1("Hfit","[0]");
  //Hratio->Fit("Hfit","r","",1.4,1.7);
  //double a2_H = 2*(p0H->GetParameter(0))/4;
  //double a2_H = (p0H->GetParameter(0));
  //cout<<"\n"<<"Helium a2: "<<a2_H<<endl;

  file << "# [Column 1: x] [Column 2: He counts] [Column 3: He error] [Column 4: d counts] [Column 5: d error] [Column 6: Hratio counts] [Column 7: Hratio error] \n";
  for(int j = 1; j <= (hH->GetXaxis()->GetNbins()); j++){

    file<<hH->GetBinCenter(j)<<" "<<hH->GetBinContent(j)<<" "<<hH->GetBinError(j)<<" "<<hD->GetBinContent(j)<<" "<<hD->GetBinError(j)<<" "<<Hratio->GetBinContent(j)<<" "<<Hratio->GetBinError(j)<<"\n";;

  }
  
  DataD->Close();
  //histfile->cd();
  //hH->Write();
  //Hratio->Write();
  file.close();
  
  return 0;
}
