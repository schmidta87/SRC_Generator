
#include <iostream>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "spectral_interpolation.h"
#include "Cross_Sections.h"
#include "constants.h"


// USER SPECIFIC: Randomly generated variabe's range (Freely change as neccessary)
  const double prob_neutron = 0.;  // probability that ejected nucleon is a neutron
  const double xB_max = 1.6;       // Maximum in X_b; unitless
  const double xB_min = 1.;        // Minimum in X_b; unitless
  const double QSq_max = 4.1;      // Maximum in Q^2; GeV^2
  const double QSq_min = 0.5;      // Minimum in Q^2; Don't go too low (must be > 0.2 for sure); GeV^2
  const double P1_xyz_mean = 0.;   // mean of initial nucleon's momentum in x,y,z
  const double P1_xyz_sigma = 0.2; // standard devidation from mean of initial nucleon's momentum in x,y,z
  const double phi_max = 2*M_PI;   // Maximum in phi: angle around Z; radians
  const double phi_min = 0.;       // Minimum in phi: angle around Z; radians

// Implicit Assumptions made in math below:
  const double Pb_x = 0.;  // initial electron momentum 'x'; Gev  (if non-zero, need to change code)
  const double Pb_y = 0.;  // initial electron momentum 'y'; Gev (if non-zero, need to change code) 

// Physical Constants
  const double MA = 4.0026032 * 0.93149 - me*2;   // Mass of initial nucleus; Helium-4  Gev
  const double BE_A = 0.028295673;                // Binding energy of initial nucleus; Helium-4; from kaeri; Gev


int main(int argc, char ** argv){
  
// Print error message if not enough / too much arguments inputted. Note: for 3 input variables, argc != 4.
// arguments: 0: unimportant; 1: initial electron z momentum P_b_z (GeV); 2: number of events to record; 3: where to store data (../DATA/gen_file.root);
  if( argc != 4)
    {
      std::cerr << "Wrong number of arguments. Instead try: \n\t"
	   << "./gen_root [Beam Energy] [Total Events] [outfile = ..DATA/gen_file.root]\n";
      
      return -1;
    }
  
   
// Input variables
  double Pb_z = atof(argv[1]);                 // initial electron momentum 'z'; usually between 4 - 12 Gev; P_z = E_b;
  int nEvents = atoi(argv[2]);                 // number of times to loop and find differential cross section
// Define changing variables for later use
  double MA_prime;                             // Mass of final nucleus; Tritium or He-3; Gev
  double Mn;                                   // Define mass of exiting nucleon
  bool isProton;                               // Define boolean of whether we have a proton or neutron
  double BE_A_prime;                           // Binding energy of recoil nucleus
// Define variables of incoming electron
  double Pb_mag = sqrt(Pb_x*Pb_x + Pb_y*Pb_y + Pb_z*Pb_z);  // Magnitude of incoming elecron momentum; Gev
  double Eb     = sqrt(Pb_mag*Pb_mag + me*me);              // initial electron beam energy (about = P_z); Gev

  
// TTree to save values during the loop
  TFile * outfile = new TFile(argv[3], "RECREATE");
  TTree * outtree = new TTree("genT", "Generator Tree");
  double weight, QSq, xB, E_miss, w, q_mag,P1_mag, Phi_Eb, pLead_mag, theta_Eb, pe_mag;
  double P1_vec[3], pe[3], pLead[3],  q_vec[3];
  int lead_type;
  
//Variables we are saving: "name", &variable, "name with elements"/double
  outtree->Branch("weight",&weight,"weight/D");             // Weighed differential cross section; proton and neutron
  outtree->Branch("xB",&xB,"xB/D");                         // Bjorken scaling variable X_b, randomly selected
  outtree->Branch("QSq",&QSq,"QSq/D");                      // Q^2, randomly selected
  outtree->Branch("lead_type",&lead_type,"lead_type/I");    // 2212 if proton, 2112 if neutron
  //outtree->Branch("E_miss",&E_miss,"E_miss/D");           // Excited energy of nucleus, goes in spectral function
  //outtree->Branch("q_vec",q_vec,"q_vec[3]/D");            // Transfered momentum to nucleon <x, y, z>
  //outtree->Branch("w",&w,"w/D");                          // Transfered Energy to nucleon

  outtree->Branch("pe",pe,"pe[3]/D");                       // Momentum of scattered electron <x, y, z>
  //outtree->Branch("theta_Eb",&theta_Eb,"theta_Eb/D");     // angle of scattered electron; Spherical: to Z-axis
  //outtree->Branch("Phi_Eb",&Phi_Eb,"Phi_Eb/D");           // Angle of scattered electron; Spherical: X-Y Plane
 
  outtree->Branch("P1_vec",P1_vec,"P1_vec[3]/D");           // Momentum of initial nucleon <x, y, z>
  outtree->Branch("pLead",pLead,"pLead[3]/D");              // Momentum of ejected nucleus <x, y, z>

  
// Initialize classes
  spec_info Spec_func;                  // Access the spectral function interpolator spec_find from the spec_info class
  Cross_Sections Sig(onshell, dipole);  // Access cross section sigma_en from the Cross_Sections class (given by andrew denniston)
  TRandom3 myRand(0);                   // Get random number

  
// loop to find the differential cross section "nEvents" amount of times.
for (int event = 0 ; event < nEvents ; event++){
  
// Reset tree variable to zero after every loop
  weight = QSq = xB = E_miss = lead_type = q_mag = w = 0.;
  pe_mag = theta_Eb = Phi_Eb = P1_mag = pLead_mag = 0.;
  memset(pLead, 0, sizeof(pLead));
  memset(pe, 0, sizeof(pe));
  memset(P1_vec,0, sizeof(P1_vec));
  memset(q_vec, 0, sizeof(q_vec));

// Define Random variables
  Phi_Eb = phi_min + myRand.Rndm()*(phi_max-phi_min);  // electron scatter angle: radians. Random value between 0 and 2*pi
  QSq = QSq_min + myRand.Rndm()*(QSq_max-QSq_min);     // Gev^2; Random; Typical: between 1 Gev^2 and 4 Gev^2
  xB = xB_min + myRand.Rndm()*(xB_max-xB_min);         // Bjorken scaling variable: unitless, random; Typical: [0.5, 1.5]
  P1_vec[0] = myRand.Gaus(P1_xyz_mean, P1_xyz_sigma);  // initial nucleon momentum 'x'; Random; Typical: between [-0.155,0.155] GeV
  P1_vec[1] = myRand.Gaus(P1_xyz_mean, P1_xyz_sigma);  // initial nucleon momentum 'y'; Random; Typical: between [-0.155,0.155] GeV
  P1_vec[2] = myRand.Gaus(P1_xyz_mean, P1_xyz_sigma);  // initial nucleon momentum 'z'; Random; Typical: between [-0.155,0.155] GeV
  int n_or_p = (int)(myRand.Rndm() + prob_neutron);    // A random number 0 (proton) or 1 (neutron)

  
// Choose the proton or neutral spectral function set (50:50 chance)
  if (n_or_p == 0){
    isProton = true;                        // Boolean input needed in sigma_eN
    Mn = 0.93827231;	                    // Mass of exiting proton; Gev
    lead_type = 2212;                       // As specified by some random institute for a proton
    MA_prime = 3.0160493 * 0.93149 - me;    // Mass of final nucleus; Tritium; Gev
    BE_A_prime = 0.008481821;}              // Binding energy of H-3; from kaeri; Gev
  else{
    isProton = false;                       // Boolean input needed in sigma_eN
    Mn = 0.93957;                           // Mass of exiting neutron; Gev
    lead_type = 2112;                       // As specified by some random institute for a neutron
    MA_prime = 3.0160293 * 0.93149 - 2*me;  // Mass of final nucleus; He-3; Gev
    BE_A_prime = 0.007718058;}              // Binding energy of He-3; from kaeri; Gev

  
// Define variables from electron scatter
  w = QSq / (2 * Mn * xB);	                // Energy transfered to nucleon from scatter event; Gev
  if (me + w > Eb) continue;                    // discontinue if unphysical (more energy transfered than initially had)
  double Eb_prime = Eb - w;	                // Energy of scattered electron; Gev
  pe_mag = sqrt(Eb_prime*Eb_prime - me*me);     // Momentum of scattered Electron; Gev
  q_mag = sqrt(QSq + w*w);                      // Momentum transfered to nucleon from photon; Gev  
  double cos_theta_Eb = (Pb_mag*Pb_mag + pe_mag*pe_mag - q_mag*q_mag)/(2*Pb_mag*pe_mag);  // cos(angle) of scattered electron with Z-axis
  if (fabs(cos_theta_Eb) > 1.) continue;        // discontinue if unphysical
  theta_Eb = acos(cos_theta_Eb);                // angle of scattered electron with z axis
  pe[0] = pe_mag*sin(theta_Eb)*cos(Phi_Eb);     // Scattered electron's X-momentum; Gev
  pe[1] = pe_mag*sin(theta_Eb)*sin(Phi_Eb);     // Scattered electron's Y-momentum; Gev
  pe[2] = pe_mag*cos(theta_Eb);                 // Scattered electron's Z-momentum; Gev
  q_vec[0] = -pe[0];                            // X-momentum transfered to nucleon from photon; Gev
  q_vec[1] = -pe[1];                            // Y-momentum transfered to nucleon from photon; Gev
  q_vec[2] = Pb_z - pe[2];                      // Z-momentum transfered to nucleon from photon Gev

  
// Define nucleus momentum and enery variables
  P1_mag = sqrt(P1_vec[0]*P1_vec[0] + P1_vec[1]*P1_vec[1] + P1_vec[2]*P1_vec[2]); // Momentum of initial nucleon to be ejected; Gev
  pLead[0] = P1_vec[0] + q_vec[0];                                         // final X-momentum of ejected nucleon; Gev   
  pLead[1] = P1_vec[1] + q_vec[1];                                         // final Y-momentum of ejected nucleon; Gev
  pLead[2] = P1_vec[2] + q_vec[2];                                         // final Z-momentum of ejected nucleon; Gev
  pLead_mag = sqrt(pLead[0]*pLead[0] + pLead[1]*pLead[1] + pLead[2]*pLead[2]);  // final Momentum of ejected nucleon; Gev
  double E1_prime = sqrt(pLead_mag*pLead_mag + Mn*Mn);  // final Energy of ejected nucleon; Gev
  double PA_prime_mag = P1_mag;                         // Momentum of final recoil nucleus; Gev
  double PA_prime_x = -P1_vec[0];                       // Momentum of final recoil nucleus 'x'; Gev
  double PA_prime_y = -P1_vec[1];                       // Momentum of final recoil nucleus 'y'; Gev
  double PA_prime_z = -P1_vec[2];                       // Momentum of final recoil nucleus 'z'; Gev
  double PA = 0;                                        // Momentum of initial nucleus; approx as zero; Gev
  double EA = sqrt(PA*PA + MA*MA);                      // Energy of initial nucleus; Gev
  double EA_prime = EA + w - E1_prime;                  // Energy of final recoil nucleus; Gev
  E_miss = Mn - (E1_prime - w);                         // Missing excitation energy for spectral function (assume PA = 0); GeV
  if ((EA_prime*EA_prime - PA_prime_mag*PA_prime_mag) < MA_prime*MA_prime) continue;  // Do not record unphysical event (need enough for rest mass)

  
// Variables that require outside calculation: spectral functin and differential cross section
  TVector3 pe_TVec(pe[0], pe[1], pe[2]); 
  TVector3 pLead_TVec(pLead[0], pLead[1], pLead[2]);    
  double sigma = Sig.sigma_eN(Eb, pe_TVec, pLead_TVec, isProton);  // cross section sigma_eN from Andrew and Jackson: Takes in (double Ebeam, TVector3 ejected electron, TVector3 ejected nucleon, bool isProton) returns in cm^2
  double spec = (pow(0.197345,3))*(0.001)*Spec_func.spec_find((double)P1_mag/0.197345, (double)(E_miss)*1000., lead_type);  // Spectral function for (P_1,E*); Uses spec_find: P_1: inverse femtometers, E*: MeV; Final unit in 1/(geV^4);

// Find differential Cross Section
  double diff_cross = spec * sigma * (w / (2 * Eb * Eb_prime * xB));   

// Creating weighting function: prob of using those exact random variables
  double P1x_prob = ((1./sqrt(2*M_PI*P1_xyz_sigma*P1_xyz_sigma))*exp(-0.5*pow((P1_vec[0] - P1_xyz_mean)/P1_xyz_sigma, 2)));
  double P1y_prob = ((1./sqrt(2*M_PI*P1_xyz_sigma*P1_xyz_sigma))*exp(-0.5*pow((P1_vec[1] - P1_xyz_mean)/P1_xyz_sigma, 2)));
  double P1z_prob = ((1./sqrt(2*M_PI*P1_xyz_sigma*P1_xyz_sigma))*exp(-0.5*pow((P1_vec[2] - P1_xyz_mean)/P1_xyz_sigma, 2)));
  double normalize_range = P1x_prob*P1y_prob*P1z_prob/((xB_max-xB_min)*(QSq_max-QSq_min)*(phi_max-phi_min));

// Divide Differential Cross Section by associated weighting
  weight = diff_cross/normalize_range;  // Weighted Differential Cross Section;

		
// Only record physical values
  if (weight > 0){
    outtree->Fill();}  
  }

    
// Clean up
  outtree->Write();
  outfile->Close();
  return 0;
  
}
