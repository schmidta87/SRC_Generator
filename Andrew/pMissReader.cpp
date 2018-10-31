#include <iostream>
#include <fstream>
#include <cstdlib>


#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"


using namespace std;


int main(int argc, char ** argv){

  if( argc != 4){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "src_reader /path/to/input/nucleus/file /path/to/output/file [A] \n";

    return -1;
    
  }
  
  //Get D and A trees. Open output file.
  TFile * DataH = new TFile(argv[1]);
  ofstream file;
  file.open(argv[2]);
  double A = atof(argv[3]);

  cerr<<"Nucleus file has been opened from: "<<argv[1]<<"\n";
  
  //Make trees and histograms for the nuclei
  TTree * TreeH = (TTree*)DataH->Get("T");
  TH1D * hH = new TH1D("Helium","Helium;xB;Counts",40,1,2);
  hH ->Sumw2();

  cerr<<"Histograms and Trees successfully created\n";
  
  //Define variables needed for histograms
  Double_t weight,QSq,xB;
  Int_t lead_type, rec_type;
  
  //Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

  //Nucleon Mass
  const double Mn = 0.938;

    cerr<<"Finished filling histogram for Deuterium\n";

  //Set addresses for H
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
  

  file << "# [Column 1: x] [Column 2: He counts] [Column 3: He error] [Column 4: d counts] [Column 5: d error] [Column 6: Hratio counts] [Column 7: Hratio error] \n";
  for(int j = 1; j <= (hH->GetXaxis()->GetNbins()); j++){

    file<<hH->GetBinCenter(j)<<" "<<hH->GetBinContent(j)<<" "<<hH->GetBinError(j)<<" "<<hD->GetBinContent(j)<<" "<<hD->GetBinError(j)<<" "<<Hratio->GetBinContent(j)<<" "<<Hratio->GetBinError(j)<<"\n";;

  }

  
  DataD->Close();
  file.close();

  cerr<< argv[3]<<" has been completed. \n\n\n";
  
  return 0;
}