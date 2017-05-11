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
  
  TFile *f1 = new TFile("../finalReachCombinedub.root");
  TFile *f2 = new TFile("../../vertex1p5/finalReachCombinedub.root");

  TFile *fout = new TFile("../finalReachCombinedAllub.root","RECREATE");
  gROOT->SetBatch(true);

  //get plot
  const int nBinsT = 198;
  const int neps = 400;
  TH2F *h_0p5 = new TH2F("h_0p5","reach 0p5mm all",nBinsT,0,0.1,neps,1E-10,4E-8);
  TH2F *h_1p5 = new TH2F("h_1p5","reach 1p5mm all",nBinsT,0,0.1,neps,1E-10,4E-8);
  TH2F *h_all = new TH2F("h_all","reach combined",nBinsT,0,0.1,neps,1E-10,4E-8);
  h_0p5 = (TH2F*)f1->Get("h5");
  h_1p5 = (TH2F*)f2->Get("h5");
  h_0p5->SetName("h_0p5");
  h_1p5->SetName("h_1p5");
  h_all->SetName("h_all");

  h_all = (TH2F*)h_0p5->Clone("h_all");
  h_all->Add(h_1p5);

  TCanvas *cc = new TCanvas("cc","Reach",1000,900);
  cc->SetLogy();
  cc->SetLogx();
  cc->cd();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(58);
  h_all->GetZaxis()->SetRangeUser(0.0,0.4);

  h_all->Draw("colz");

  h_all->GetXaxis()->SetRangeUser(1E-2,1E-1);
  h_all->GetXaxis()->SetTitle("Mass [GeV]");
  h_all->GetXaxis()->SetTitleOffset(1.3);
  h_all->GetYaxis()->SetTitle("#epsilon^{2}");
  cc->Update();
  cc->SaveAs("../reachAllDataBothSetsub.C");
  cc->Close();
  
  
  fout->Write();
  fout->Close(); 


}
