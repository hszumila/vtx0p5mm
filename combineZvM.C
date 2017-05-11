#include "TH1.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TTree.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TROOT.h"
#include <vector>
#include "unistd.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH2.h"
#include "TCut.h"
#include "Riostream.h"
#include "TNtupleD.h"
#include "TCut.h"
#include "TMultiGraph.h"
#include "TLine.h"

using namespace std;
int main(){
  
  TFile *f1 = new TFile("../output_L1L1/ub/cutPlots.root");
  TFile *f2 = new TFile("../output_L1L2/ub/cutPlots.root");
  TFile *f3 = new TFile("../output_L2L2/ub/cutPlots.root");

  TFile *fout = new TFile("../finalZvMub.root","RECREATE");
  gROOT->SetBatch(true);

  //get plot
  TH2F *h1 = new TH2F("h1","zVertex vs Mass",200,0,0.1,200,-100,100);
  h1 = (TH2F*)f1->Get("h_zvm");
  TH2F *h2 = new TH2F("h2","zVertex vs Mass",200,0,0.1,200,-100,100);
  h2 = (TH2F*)f2->Get("h_zvm");
  TH2F *h3 = new TH2F("h3","zVertex vs Mass",200,0,0.1,200,-100,100);
  h3 = (TH2F*)f3->Get("h_zvm");
  h1->SetName("h1");
  h2->SetName("h2");
  h3->SetName("h3");

  TH2F *h4 = new TH2F("h4","combined Z vs M",200,0,0.1,200,-100,100);
  h4 = (TH2F*)h1->Clone("h4");
  h4->Add(h2);
  h4->Add(h3);

  
  
  fout->Write();
  fout->Close(); 


}
