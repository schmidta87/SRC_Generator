#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using namespace std;

const double Ebeam=5.766;
const double deg_margin=1.;

int main(int argc, char ** argv)
{						
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tfomin /path/to/input /path/to/output\n\n";
      return -1;
    }

  TFile * inF = new TFile(argv[1]);
  TTree * inT = (TTree*) inF->Get("T");
  const int nEvents = inT->GetEntries();


  cout << "Hello world.\n";
  return 0;
}
