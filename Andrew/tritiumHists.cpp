#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TVector3.h"

#include "constants.h"
#include "helpers.h"
#include "Nuclear_Info.h"
#include "Cross_Sections.h"

using namespace std;

double centThetaE = 20.88;
double centThetaPLow = 48.82;
double centThetaPHigh = 58.50;
double dTheta = (180/M_PI) * 0.0275;
double centMomE = 3.543;
double centMomPLow = 1.481;
double centMomPHigh = 1.246;
double dMom = 0.04;

bool checkSpot(int lh, TVector3 vq, TVector3 ve, TVector3 vLead, TVector3 vRRec);

int main(int argc, char ** argv){

  if( argc != 5){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "tritiumHists /path/to/input/Tritium/file /path/to/input/Helium3/file /path/to/output/file Ebeam\n";

    return -1;
    
  }
  
  TFile * ftri = new TFile(argv[1]);
  TFile * fhe3 = new TFile(argv[2]);
  TFile * fout = new TFile(argv[3],"RECREATE");
  const TVector3 vEbeam(0,0,atoi(argv[4]));
  
  cerr<<"Files have been opened from: "<<argv[1]<<", "<<argv[2]<<"\n";

  vector<TH1*> histList;

  //Create Histograms
  //First element corresponds to Tritium (0) and Helium(1). The second corresponds to low(0) and high(1) pmiss kinematics
  /*  TH1D * H_Emiss[2][2];
  TH1D * H_pmiss[2][2];
  TH2D * H_QSq[2][2];

  for( int i = 0; i < 2; i++){
    for( int j = 0; j < 2; j++){
      char temp[100];

      sprintf(temp,"Emiss_%d_%d",i,j);
      H_Emiss[i][j] = new TH1D(temp,"Emiss;Emiss [GeV];Counts",40,0,0.1);
      histList.push_back(h1p_Emiss_split[i][j]);
		
    }
    }*/
  
  TH1D * tri_Emiss_low = new TH1D("tri_Emiss_low","Triton;Emiss_low [GeV];Counts",40,0,1);
  histList.push_back(tri_Emiss_low);
  TH1D * tri_Emiss_high = new TH1D("tri_Emiss_high","Triton;Emiss_high [GeV];Counts",40,0,1);
  histList.push_back(tri_Emiss_high);
  TH1D * tri_pmiss_low = new TH1D("tri_pmiss_low","Triton;pmiss_low [GeV];Counts",40,0,1);
  histList.push_back(tri_pmiss_low);
  TH1D * tri_pmiss_high = new TH1D("tri_pmiss_high","Triton;pmiss_high [GeV];Counts",40,0,1);
  histList.push_back(tri_pmiss_high);  
  TH1D * tri_QSq_low = new TH1D("tri_QSq_low","Triton;QSq_low [GeV^2];Counts",40,1.3,2.4);  
  histList.push_back(tri_QSq_low);
  TH1D * tri_QSq_high = new TH1D("tri_QSq_high","Triton;QSq_high [GeV^2];Counts",40,1.3,2.4);  
  histList.push_back(tri_QSq_high);

  TH1D * he3_Emiss_low = new TH1D("he3_Emiss_low","Helium-3;Emiss_low [GeV];Counts",40,0,0.1);
  histList.push_back(he3_Emiss_low);
  TH1D * he3_Emiss_high = new TH1D("he3_Emiss_high","Helium-3;Emiss_high [GeV];Counts",40,0,0.1);
  histList.push_back(he3_Emiss_high);
  TH1D * he3_pmiss_low = new TH1D("he3_pmiss_low","Helium-3;pmiss_low [GeV];Counts",40,0,0.1);
  histList.push_back(he3_pmiss_low);
  TH1D * he3_pmiss_high = new TH1D("he3_pmiss_high","Helium-3;pmiss_high [GeV];Counts",40,0,0.1);
  histList.push_back(he3_pmiss_high);  
  TH1D * he3_QSq_low = new TH1D("he3_QSq_low","Helium-3;QSq_low [GeV^2];Counts",40,1.3,2.4);  
  histList.push_back(he3_QSq_low);
  TH1D * he3_QSq_high = new TH1D("he3_QSq_high","Helium-3;QSq_high [GeV^2];Counts",40,1.3,2.4);  
  histList.push_back(he3_QSq_high);

  //Add weights
  for( int i=0; i <histList.size() ; i++) histList[i]->Sumw2();

  //Get tri to he3 ratio
  TGraphAsymmErrors * tri_to_he3_low = new TGraphAsymmErrors();
  tri_to_he3_low->SetName("tri_he3_ratio_low");
  tri_to_he3_low->SetTitle("tri/he3_low;p_miss;Cross Section Ratio");

  TGraphAsymmErrors * tri_to_he3_high = new TGraphAsymmErrors();
  tri_to_he3_high->SetName("tri_he3_ratio_high");
  tri_to_he3_high->SetTitle("tri/he3_high;p_miss;Cross Section Ratio");

  cerr<<"Histograms and TGraphAsymmErrors have been made\n";
  
  //Make trees and histograms for the nuclei
  TTree * triTree = (TTree*)ftri->Get("genT");
  TTree * he3Tree = (TTree*)fhe3->Get("genT");

  cerr<<"triTree successfully identified\n";
  
  //Define variables needed for histograms of Helium
  Double_t pe[3],pLead[3],pRec[3],weight;
  Int_t lead_type, rec_type;
    
  //Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

  //Set addresses for triTree
  triTree->SetBranchAddress("lead_type",&lead_type);
  cerr<<"triTree variable lead type successfully read\n";
  triTree->SetBranchAddress("rec_type",&rec_type);
  triTree->SetBranchAddress("pe",&pe);
  triTree->SetBranchAddress("pLead",&pLead);
  triTree->SetBranchAddress("pRec",&pRec);
  triTree->SetBranchAddress("weight",&weight);

  cerr<<"triTree successfully read\n";

  //Loop over triTree
  for(int i = 0; i < triTree->GetEntries(); i++){

    if (i %100000 ==0){
      cerr<<"Working on event "<<i<<" of Triton Tree\n";
    }
    triTree->GetEntry(i);

    bool isLow = false;
    bool isHigh = false;
    
    TVector3 ve(pe[0],pe[1],pe[2]);
    TVector3 vLead(pLead[0],pLead[1],pLead[2]);
    TVector3 vq = vEbeam-ve;
    TVector3 vMiss = vLead - vq;
    TVector3 vRRec = -vMiss;

    double ELead = sqrt(vLead.Mag2()+mN*mN);
    double omega = vEbeam.Mag() - ve.Mag();
    double pmiss = vMiss.Mag();
    double Emiss = vMiss.Mag();
    double QSq = vq.Mag2() - omega*omega;
    double xB = QSq/(2*mN*omega);

    isLow=checkSpot(0,vq,ve,vLead,vRRec);
    isHigh=checkSpot(1,vq,ve,vLead,vRRec);
    
    if(isLow){
      tri_Emiss_low->Fill(Emiss,weight);
      tri_pmiss_low->Fill(pmiss,weight);
      tri_QSq_low->Fill(QSq,weight);     
    }
    else if(isHigh){
      tri_Emiss_high->Fill(Emiss,weight);
      tri_pmiss_high->Fill(pmiss,weight);
      tri_QSq_high->Fill(QSq,weight);     

    }

  }
      cerr<<"Finished filling histograms for Triton\n";


  //Set addresses for he3Tree
  he3Tree->SetBranchAddress("lead_type",&lead_type);
  he3Tree->SetBranchAddress("rec_type",&rec_type);
  he3Tree->SetBranchAddress("pe",&pe);
  he3Tree->SetBranchAddress("pLead",&pLead);
  he3Tree->SetBranchAddress("pRec",&pRec);
  he3Tree->SetBranchAddress("weight",&weight);

  cerr<<"triTree successfully read\n";

  //Loop over he3Tree
  for(int i = 0; i < he3Tree->GetEntries(); i++){

    if (i %100000 ==0){
      cerr<<"Working on event "<<i<<" of Helium-3 Tree\n";
    }
    he3Tree->GetEntry(i);

    bool isLow = false;
    bool isHigh = false;
    
    TVector3 ve(pe[0],pe[1],pe[2]);
    TVector3 vLead(pLead[0],pLead[1],pLead[2]);
    TVector3 vq = vEbeam-ve;
    TVector3 vMiss = vLead - vq;
    TVector3 vRRec = -vMiss;

    double ELead = sqrt(vLead.Mag2()+mN*mN);
    double omega = vEbeam.Mag() - ve.Mag();
    double pmiss = vMiss.Mag();
    double Emiss = vMiss.Mag();
    double QSq = vq.Mag2() - omega*omega;
    double xB = QSq/(2*mN*omega);

    //Get Angles for Cuts
    double theta_Recq = (180/M_PI) * vRRec.Angle(vq);
    double thetaE = (180/M_PI) * ve.Theta();
    double thetaP = (180/M_PI) * vLead.Theta();

    if(theta_Recq > 37.5) continue;
    if(xB < 1.3) continue;
    if( (thetaE < (centThetaE-dTheta)) || (thetaE > (centThetaE+dTheta)) ) continue;
    if( (ve.Mag() < (centMomE*(1-dMom))) || (ve.Mag() > (centMomE*(1+dMom))) ) continue;
    
    if( (thetaP > (centThetaPLow - dTheta)) && (thetaP < (centThetaPLow + dTheta)) ){
      if( (vLead.Mag() > (centMomPLow*(1-dMom))) && (vLead.Mag() < (centMomPLow*(1+dMom))) ){
	isLow = true;
      }
    } 

    if( (thetaP > (centThetaPHigh - dTheta)) && (thetaP < (centThetaPHigh + dTheta)) ){
      if( (vLead.Mag() > (centMomPHigh*(1-dMom))) && (vLead.Mag() < (centMomPHigh*(1+dMom))) ){
	isHigh = true;
      }
    } 

    
    if(isLow){
      he3_Emiss_low->Fill(Emiss,weight);
      he3_pmiss_low->Fill(pmiss,weight);
      he3_QSq_low->Fill(QSq,weight);     
    }
    else if(isHigh){
      he3_Emiss_high->Fill(Emiss,weight);
      he3_pmiss_high->Fill(pmiss,weight);
      he3_QSq_high->Fill(QSq,weight);     

    }

  }

  cerr<<"Finished filling histograms for Helium-3\n";


  tri_to_he3_low->BayesDivide(tri_pmiss_low,he3_pmiss_low);
  tri_to_he3_high->BayesDivide(tri_pmiss_high,he3_pmiss_high);

  cerr<<"Finished filling BayesDivide\n";

  //Write Out
  fout->cd();
  tri_to_he3_low->Write();
  tri_to_he3_high->Write();

   
  for( int i = 0; i < histList.size(); i++){
       histList[i]->Write();
  }

   cerr<<"Finished writing Histos and TGAEs\n";

  fout->Close();

  return 0;
}


bool checkSpot(int lh, TVector3 vq, TVector3 ve, TVector3 vLead, TVector3 vRRec){

  //Get Angles for Cuts
  double theta_Recq = (180/M_PI) * vRRec.Angle(vq);
  double thetaE = (180/M_PI) * ve.Theta();
  double thetaP = (180/M_PI) * vLead.Theta();
  
  //if(theta_Recq > 37.5) return false;
  //if(xB < 1.3) return false;
  //if( (thetaE < (centThetaE-dTheta)) || (thetaE > (centThetaE+dTheta)) ) return false;
  //if( (ve.Mag() < (centMomE*(1-dMom))) || (ve.Mag() > (centMomE*(1+dMom))) ) return false;

  if(lh = 0){
    if( (thetaP > (centThetaPLow - dTheta)) && (thetaP < (centThetaPLow + dTheta)) ){
      if( (vLead.Mag() > (centMomPLow*(1-dMom))) && (vLead.Mag() < (centMomPLow*(1+dMom))) ){
	return true;
      }
    } 
  }
  else if(lh = 1){
    if( (thetaP > (centThetaPHigh - dTheta)) && (thetaP < (centThetaPHigh + dTheta)) ){
      if( (vLead.Mag() > (centMomPHigh*(1-dMom))) && (vLead.Mag() < (centMomPHigh*(1+dMom))) ){
	return true;
      }
    } 
  }
  return false;
}

