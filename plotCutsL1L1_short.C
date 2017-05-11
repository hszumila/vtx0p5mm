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

using namespace std;


int main(){

  //input and output files
  TFile *f = new TFile("../data/hps_tkp6_L1L1_15Mar.root");
  TTree *tt=(TTree*)f->Get("ntuple");
  TFile *fout = new TFile("../output_5772/cutPlots_short.root","RECREATE");

  //cuts listed in order
  TCut cut1 = "eleHasL1&&posHasL1&&eleHasL2&&posHasL2&&uncP>0.845&&isPair1&&eleClY*posClY<0";
  TCut cut2 = "max(eleTrkChisq,posTrkChisq)<30";
  TCut cut3 = "eleP<0.75*1.056";
  TCut cut4 = "min(eleMinPositiveIso+0.5*(eleTrkZ0-5*elePY/eleP)*sign(elePY),posMinPositiveIso+0.5*(posTrkZ0-5*posPY/posP)*sign(posPY))>0";
  TCut cut5 = "bscChisq<10";
  TCut cut6 = "bscChisq-uncChisq<5";
  TCut cut7 = "uncP<1.15*1.056";
  TCut cut8 = "max(eleMatchChisq,posMatchChisq)<10";
  TCut cut9 = "abs(eleClT-eleTrkT-43)<4&&abs(posClT-posTrkT-43)<4";
  TCut cut10 = "abs(eleClT-posClT)<2";
  TCut cut11 = "abs(eleP-posP)/(eleP+posP)<0.4";
  TCut cut12 = "posTrkD0<1.5";
  TCut cut13 = "(posMaxHitsShared<5&&eleMaxHitsShared<5)";//(nPos==1||(posMaxHitsShared<4&&eleMaxHitsShared<4))";
  
  gROOT->SetBatch(true);

 

  //plot vz:mass
  TH2F *h_zvm = new TH2F("h_zvm","z vtx vs mass",200,0,0.1,200,-100,100);
  tt->Draw("uncVZ:uncM>>h_zvm",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12 && cut13);
  

 fout->Write();
 fout->Close();
 
}
