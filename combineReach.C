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
  
  TFile *f1 = new TFile("../output_L1L1/ub/outputFits.root");
  TFile *f2 = new TFile("../output_L1L2/ub/outputFits.root");
  TFile *f3 = new TFile("../output_L2L2/ub/outputFits.root");

  TFile *fout = new TFile("../finalReachCombinedub.root","RECREATE");
  gROOT->SetBatch(true);

  //get plot
  const int nBinsT = 198;
  const int neps = 400;
  TH2F *h1 = new TH2F("h1","reach 0p5mm L1L1",nBinsT,0,0.1,neps,1E-10,4E-8);
  TH2F *h2 = new TH2F("h2","reach 0p5mm L1L2",nBinsT,0,0.1,neps,1E-10,4E-8);
  TH2F *h3 = new TH2F("h3","reach 0p5mm L2L2",nBinsT,0,0.1,neps,1E-10,4E-8);
  TH2F *h5 = new TH2F("h5","reach combined",nBinsT,0,0.1,neps,1E-10,4E-8);
  h1 = (TH2F*)f1->Get("h_reach");
  h2 = (TH2F*)f2->Get("h_reach");
  h3 = (TH2F*)f3->Get("h_reach");

  h1->SetName("h1");
  h2->SetName("h2");
  h3->SetName("h3");
  h5->SetName("h5");

  h5 = (TH2F*)h1->Clone("h5");
  h5->Add(h2);
  h5->Add(h3);

  TCanvas *cc = new TCanvas("cc","Reach",1000,900);
  cc->SetLogy();
  cc->SetLogx();
  cc->cd();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(58);
  h5->GetZaxis()->SetRangeUser(0.0,0.25);

  h5->Draw("colz");

  h5->GetXaxis()->SetRangeUser(1E-2,1E-1);
  h5->GetXaxis()->SetTitle("Mass [GeV]");
  h5->GetXaxis()->SetTitleOffset(1.3);
  h5->GetYaxis()->SetTitle("#epsilon^{2}");
  cc->Update();
  cc->SaveAs("../reach0p5ub.C");
  cc->Close();
  
  
  fout->Write();
  fout->Close(); 


}
