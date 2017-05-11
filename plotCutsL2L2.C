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
  TFile *f = new TFile("../data/hps_L2L2_blind_10Feb17.root");
  TTree *tt=(TTree*)f->Get("ntuple");
  TFile *fout = new TFile("../output_L2L2/cutPlots.root","RECREATE");

  //cuts listed in order
  TCut cut1 = "!eleHasL1&&!posHasL1&&eleHasL2&&posHasL2&&uncP>0.845&&isPair1&&abs(posFirstHitY-100*posPY/posP)<1.5&&abs(eleFirstHitY-100*elePY/eleP)<1.5";
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
  TCut cut13 = "(posMaxHitsShared<4&&eleMaxHitsShared<4)";//(nPos==1||(posMaxHitsShared<4&&eleMaxHitsShared<4))";
  
  gROOT->SetBatch(true);

  //plots for standardizing:
  TH1F *h_tdiff1 = new TH1F("h_tdiff1","t1-t2",200,-9,9);
  TH1F *h_tdiff2 = new TH1F("h_tdiff2","t1-t2",200,-9,9);
  TH1F *h_tdiff3 = new TH1F("h_tdiff3","t1-t2",200,-9,9);
  TH1F *h_tdiff4 = new TH1F("h_tdiff4","t1-t2",200,-9,9);
  TH1F *h_tdiff5 = new TH1F("h_tdiff5","t1-t2",200,-9,9);
  TH1F *h_tdiff6 = new TH1F("h_tdiff6","t1-t2",200,-9,9);
  TH1F *h_tdiff7 = new TH1F("h_tdiff7","t1-t2",200,-9,9);
  TH1F *h_tdiff8 = new TH1F("h_tdiff8","t1-t2",200,-9,9);
  TH1F *h_tdiff9 = new TH1F("h_tdiff9","t1-t2",200,-9,9);
  TH1F *h_tdiff11 = new TH1F("h_tdiff11","t1-t2",200,-9,9);
  TH1F *h_tdiff12 = new TH1F("h_tdiff12","t1-t2",200,-9,9);
  TH1F *h_tdiff13 = new TH1F("h_tdiff13","t1-t2",200,-9,9);
  TH1F *h_tdiff14 = new TH1F("h_tdiff14","t1-t2",200,-9,9);

  TH1F *h_z1 = new TH1F("h_z1","zvtx",200,-100,100);
  TH1F *h_z2 = new TH1F("h_z2","zvtx",200,-100,100);
  TH1F *h_z3 = new TH1F("h_z3","zvtx",200,-100,100);
  TH1F *h_z4 = new TH1F("h_z4","zvtx",200,-100,100);
  TH1F *h_z5 = new TH1F("h_z5","zvtx",200,-100,100);
  TH1F *h_z6 = new TH1F("h_z6","zvtx",200,-100,100);
  TH1F *h_z7 = new TH1F("h_z7","zvtx",200,-100,100);
  TH1F *h_z8 = new TH1F("h_z8","zvtx",200,-100,100);
  TH1F *h_z9 = new TH1F("h_z9","zvtx",200,-100,100);
  TH1F *h_z10 = new TH1F("h_z10","zvtx",200,-100,100);
  TH1F *h_z11 = new TH1F("h_z11","zvtx",200,-100,100);
  TH1F *h_z12 = new TH1F("h_z12","zvtx",200,-100,100);
  TH1F *h_z13 = new TH1F("h_z13","zvtx",200,-100,100);
  TH1F *h_z14 = new TH1F("h_z14","zvtx",200,-100,100);

  TH1F *h_m1 = new TH1F("h_m1","mass",200,0,0.1);
  TH1F *h_m2 = new TH1F("h_m2","mass",200,0,0.1);
  TH1F *h_m3 = new TH1F("h_m3","mass",200,0,0.1);
  TH1F *h_m4 = new TH1F("h_m4","mass",200,0,0.1);
  TH1F *h_m5 = new TH1F("h_m5","mass",200,0,0.1);
  TH1F *h_m6 = new TH1F("h_m6","mass",200,0,0.1);
  TH1F *h_m7 = new TH1F("h_m7","mass",200,0,0.1);
  TH1F *h_m8 = new TH1F("h_m8","mass",200,0,0.1);
  TH1F *h_m9 = new TH1F("h_m9","mass",200,0,0.1);
  TH1F *h_m10 = new TH1F("h_m10","mass",200,0,0.1);
  TH1F *h_m11 = new TH1F("h_m11","mass",200,0,0.1); 
  TH1F *h_m12 = new TH1F("h_m12","mass",200,0,0.1); 
  TH1F *h_m13 = new TH1F("h_m13","mass",200,0,0.1); 
  TH1F *h_m14 = new TH1F("h_m14","mass",200,0,0.1); 

  Double_t signal[13];
  Double_t background[13];
  
  //cut on track chi2
  TH1F *h_trkChi2 = new TH1F("h_trkChi2","trk chi2",200,0,50);
  TH2F *h2_trkChi_zCut = new TH2F("h2_trkChi_zCut","trk #chi^{2} vs Zvtx",200,0,50,200,-5,85);
  tt->Draw("uncVZ:eleTrkChisq>>h2_trkChi_zCut", cut1 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 );
  tt->Draw("uncVZ:posTrkChisq>>h2_trkChi_zCut", cut1 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 );
  tt->Draw("eleTrkChisq>>h_trkChi2",cut1);
  tt->Draw("posTrkChisq>>h_trkChi2",cut1);
  tt->Draw("uncVZ>>h_z1",cut1);
  tt->Draw("(eleClT-posClT)>>h_tdiff1",cut1);
  tt->Draw("uncM>>h_m1",cut1);
  signal[0] = h_z1->Integral(h_z1->FindBin(-5),h_z1->FindBin(7));
  background[0] = h_z1->Integral(h_z1->FindBin(-5),h_z1->FindBin(85))-signal[0];	      
  cout<<"cut 1"<<endl;

  //cut on Pmax for single track
  TH1F *h_eleMom = new TH1F("h_eleMom","e- momentum",200,0,1.3);
  tt->Draw("eleP>>h_eleMom",cut1 && cut2);
  tt->Draw("uncVZ>>h_z2",cut1 && cut2);
  tt->Draw("(eleClT-posClT)>>h_tdiff2",cut1 && cut2);
  tt->Draw("uncM>>h_m2",cut1 && cut2);
  signal[1] = h_z2->Integral(h_z2->FindBin(-5),h_z2->FindBin(7));
  background[1] = h_z2->Integral(h_z2->FindBin(-5),h_z2->FindBin(85))-signal[1];
  cout<<"cut 2"<<endl;

  //apply isolation cut
  TH1F *h_iso = new TH1F("h_iso","isolation",200,-20,20);
  tt->Draw("min(eleMinPositiveIso+0.5*(eleTrkZ0-5*elePY/eleP)*sign(elePY),posMinPositiveIso+0.5*(posTrkZ0-5*posPY/posP)*sign(posPY))>>h_iso",cut1 && cut2 && cut3);
  tt->Draw("uncVZ>>h_z3",cut1 && cut2 && cut3);
  tt->Draw("(eleClT-posClT)>>h_tdiff3",cut1 && cut2 && cut3);
  tt->Draw("uncM>>h_m3",cut1 && cut2 && cut3);
  signal[2] = h_z3->Integral(h_z3->FindBin(-5),h_z3->FindBin(7));
  background[2] = h_z3->Integral(h_z3->FindBin(-5),h_z3->FindBin(85))-signal[2];
  cout<<"cut 3"<<endl;

  //apply vertex chi2 cuts
  TH1F *h_bscChi2 = new TH1F("h_bscChi2","bsc chi2",200,0,30);
  TH2F *h2_bscChi_zCut = new TH2F("h2_bscChi_zCut","bsc #chi^{2} vs Zvtx",200,0,30,200,-5,85);
  tt->Draw("uncVZ:bscChisq>>h2_bscChi_zCut",cut1 && cut2 && cut3 && cut4 && cut7 && cut8 && cut9 && cut10 && cut11 );
  tt->Draw("bscChisq>>h_bscChi2",cut1 && cut2 && cut3 && cut4);
  tt->Draw("uncVZ>>h_z4",cut1 && cut2 && cut3 && cut4);
  tt->Draw("(eleClT-posClT)>>h_tdiff4",cut1 && cut2 && cut3 && cut4);
  tt->Draw("uncM>>h_m4",cut1 && cut2 && cut3 && cut4);
  signal[3] = h_z4->Integral(h_z4->FindBin(-5),h_z4->FindBin(7));
  background[3] = h_z4->Integral(h_z4->FindBin(-5),h_z4->FindBin(85))-signal[3];
  cout<<"cut 4"<<endl;

  TH1F *h_bmuChi2 = new TH1F("h_bmuChi2","bsc - unc chi2",200,0,20);
  tt->Draw("(bscChisq-uncChisq)>>h_bmuChi2",cut1 && cut2 && cut3 && cut4 && cut5);
  tt->Draw("uncVZ>>h_z5",cut1 && cut2 && cut3 && cut4 && cut5);
  tt->Draw("(eleClT-posClT)>>h_tdiff5",cut1 && cut2 && cut3 && cut4 && cut5);
  tt->Draw("uncM>>h_m5",cut1 && cut2 && cut3 && cut4 && cut5);
  signal[4] = h_z5->Integral(h_z5->FindBin(-5),h_z5->FindBin(7));
  background[4] = h_z5->Integral(h_z5->FindBin(-5),h_z5->FindBin(85))-signal[4];
  cout<<"cut 5"<<endl;

  //apply Psum max cut
  TH1F *h_pSum = new TH1F("h_pSum","momentum sum",200,0,1.3);
  tt->Draw("uncP>>h_pSum", cut1 && cut2 && cut3 && cut4 && cut5 && cut6);
  tt->Draw("uncVZ>>h_z6",cut1 && cut2 && cut3 && cut4 && cut5 && cut6);
  tt->Draw("(eleClT-posClT)>>h_tdiff6",cut1 && cut2 && cut3 && cut4 && cut5 && cut6);
  tt->Draw("uncM>>h_m6",cut1 && cut2 && cut3 && cut4 && cut5 && cut6);
  signal[5] = h_z6->Integral(h_z6->FindBin(-5),h_z6->FindBin(7));
  background[5] = h_z6->Integral(h_z6->FindBin(-5),h_z6->FindBin(85))-signal[5];
  cout<<"cut 6"<<endl;

  //apply chi2 matching cut
  TH1F *h_MmatchChi2 = new TH1F("h_MmatchChi2","e- matching chi2",100,0,20);
  TH1F *h_PmatchChi2 = new TH1F("h_PmatchChi2","e+ matching chi2",100,0,20);
  TH2F *h2_matchChi_zCut = new TH2F("h2_matchChi_zCut","matching #chi^{2} vs Zvtx",200,0,20,200,-5,85);
  tt->Draw("uncVZ:eleMatchChisq>>h2_matchChi_zCut", cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut9 && cut10 && cut11 );
  tt->Draw("uncVZ:posMatchChisq>>h2_matchChi_zCut", cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut9 && cut10 && cut11 );
  tt->Draw("eleMatchChisq>>h_MmatchChi2", cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7);
  tt->Draw("posMatchChisq>>h_PmatchChi2", cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7);
  tt->Draw("uncVZ>>h_z7",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7);
  tt->Draw("(eleClT-posClT)>>h_tdiff7",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7);
  tt->Draw("uncM>>h_m7",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7);
  signal[6] = h_z7->Integral(h_z7->FindBin(-5),h_z7->FindBin(7));
  background[6] = h_z7->Integral(h_z7->FindBin(-5),h_z7->FindBin(85))-signal[6];
  cout<<"cut 7"<<endl;

  //apply track-ecal timing cut
  TH1F *h_cltkTdiff = new TH1F("h_cltkTdiff","cluster, track time difference - 43 ns",200,-10,10);
  tt->Draw("(eleClT-eleTrkT-43)>>h_cltkTdiff",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8);
  tt->Draw("(posClT-posTrkT-43)>>h_cltkTdiff",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8);
  tt->Draw("uncVZ>>h_z8",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8);
  tt->Draw("(eleClT-posClT)>>h_tdiff8",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8);
  tt->Draw("uncM>>h_m8",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8);
  signal[7] = h_z8->Integral(h_z8->FindBin(-5),h_z8->FindBin(7));
  background[7] = h_z8->Integral(h_z8->FindBin(-5),h_z8->FindBin(85))-signal[7];
  cout<<"cut 8"<<endl;
    
  //apply cluster time difference cut
  TH1F *h_clTdiff  = new TH1F("h_clTdiff","2 cluster time difference",200,-8,8);
  tt->Draw("(eleClT-posClT)>>h_clTdiff",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9);
  tt->Draw("(eleClT-posClT)>>h_tdiff9",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9);
  tt->Draw("uncVZ>>h_z9",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9);
  tt->Draw("uncM>>h_m9",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9);
  signal[8] = h_z9->Integral(h_z9->FindBin(-5),h_z9->FindBin(7));
  background[8] = h_z9->Integral(h_z9->FindBin(-5),h_z9->FindBin(85))-signal[8];
  cout<<"cut 9"<<endl;

  //apply momentum asymmetry cut
  TH1F *h_momAsy = new TH1F("h_momAsy","momentum asymmetry",200,0,1);
  TH2F *h_pSumDiff = new TH2F("h_pSumDiff","momentum sum vs difference",200,0,1.3,200,0,0.65);
  tt->Draw("abs(eleP-posP)/(eleP+posP)>>h_momAsy",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10);
  tt->Draw("abs(eleP-posP):uncP>>h_pSumDiff",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10);
  tt->Draw("uncVZ>>h_z10",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10);
  tt->Draw("uncM>>h_m10",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10);
  signal[9] = h_z10->Integral(h_z10->FindBin(-5),h_z10->FindBin(7));
  background[9] = h_z10->Integral(h_z10->FindBin(-5),h_z10->FindBin(85))-signal[9];
  cout<<"cut 10"<<endl;

   //plot coplanarity
  TH1F *h_coplanarity = new TH1F("h_coplanarity","coplanarity",200,0,300);
  tt->Draw("(atan(eleClY/(eleClX-42.5))-atan(posClY/(posClX-42.5)))*180/3.14+180>>h_coplanarity",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 );
  tt->Draw("uncVZ>>h_z11",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11);
  tt->Draw("uncM>>h_m11",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11);
  tt->Draw("(eleClT-posClT)>>h_tdiff11",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 &&cut11);
  signal[10] = h_z11->Integral(h_z11->FindBin(-5),h_z11->FindBin(7));
  background[10] = h_z11->Integral(h_z11->FindBin(-5),h_z11->FindBin(85))-signal[10];
  cout<<"cut 11"<<endl;

  TH1F *h_coplanaritySVT = new TH1F("h_coplanaritySVT","coplanarity",200,0,300);
  tt->Draw("(atan(eleTrkEcalY/(eleTrkEcalX-42.5))-atan(posTrkEcalY/(posTrkEcalX-42.5)))*180/3.14+180>>h_coplanaritySVT",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11);

  //cut e+ track d0
  TH1F *h_d0 = new TH1F("h_d0","e+ track d0",200,-4,4);
  tt->Draw("posTrkD0>>h_d0",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 &&cut11 &&cut12);
  tt->Draw("uncVZ>>h_z12",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&&cut12);
  tt->Draw("uncM>>h_m12",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 &&cut12);
  tt->Draw("(eleClT-posClT)>>h_tdiff12",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut11 &&cut12);
  cout<<"cut 12"<<endl;


  signal[11] = h_z12->Integral(h_z12->FindBin(-5),h_z12->FindBin(7));
  background[11] = h_z12->Integral(h_z12->FindBin(-5),h_z12->FindBin(85))-signal[11];

  //cut max shared hits
  TH1F *h_PmaxShared = new TH1F("h_PmaxShared","e+ max shared hits",10,0,7);
  TH1F *h_MmaxShared = new TH1F("h_MmaxShared","e- max shared hits",10,0,7);
  tt->Draw("posMaxHitsShared>>h_PmaxShared",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 &&cut11 && cut12 &&cut13);
  tt->Draw("eleMaxHitsShared>>h_MmaxShared",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 &&cut11 && cut12&&cut13);
  tt->Draw("uncVZ>>h_z13",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12&&cut13);
  tt->Draw("uncM>>h_m13",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12 &&cut13);
  tt->Draw("(eleClT-posClT)>>h_tdiff13",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut11 &&cut12 &&cut13);
  signal[12] = h_z13->Integral(h_z13->FindBin(-5),h_z13->FindBin(7));
  background[12] = h_z13->Integral(h_z13->FindBin(-5),h_z13->FindBin(85))-signal[12];
  cout<<"cut 13"<<endl;

  /////////////////
 
  
  //plot vy:vx
  TH2F *h_yvx = new TH2F("h_yvx","y vtx vs x vtx",300,-3,3,300,-3,3);
  tt->Draw("uncVY:uncVX>>h_yvx",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12 && cut13);
  
  //plot vy:mass
  TH2F *h_yvm = new TH2F("h_yvm","y vtx vs mass",300,0,0.1,300,-3,3);
  tt->Draw("uncVY:uncM>>h_yvm",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12 && cut13);
  
  //plot vx:mass
  TH2F *h_xvm = new TH2F("h_xvm","x vtx vs mass",300,0,0.1,300,-3,3);
  tt->Draw("uncVX:uncM>>h_xvm",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12 && cut13);

  //plot vz:mass
  TH2F *h_zvm = new TH2F("h_zvm","z vtx vs mass",400,0,0.1,400,-100,100);
  tt->Draw("uncVZ:uncM>>h_zvm",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12 && cut13);
  
  //Cut correlation plots
  TH2F *h2_bscChi_bscuncChi = new TH2F("h2_bscChi_bscuncChi","#chi_{bsc}^{2} vs #chi^2 difference (bsc-unc)",200,0,20,200,0,30);
  tt->Draw("bscChisq:(bscChisq-uncChisq)>>h2_bscChi_bscuncChi",cut1);
  TH2F *h2_bscChi_Psum = new TH2F("h2_bscChi_Psum","#chi_{bsc}^{2} vs P_{sum}",200,0,1.3,200,0,30);
  tt->Draw("bscChisq:uncP>>h2_bscChi_Psum",cut1);
  TH2F *h2_matchChi_bscChi = new TH2F("h2_matchChi_bscChi","matching #chi^{2} vs #chi^{2}_{bsc}",200,0,30,200,0,30);
  tt->Draw("eleMatchChisq:bscChisq>>h2_matchChi_bscChi",cut1);
  TH2F *h2_matchChi_tecalTcut = new TH2F("h2_matchChi_tecalTcut","matching #chi^{2} vs track-ecal time difference",200,-10,10,200,0,30);
  tt->Draw("eleMatchChisq:(eleClT-eleTrkT-43)>>h2_matchChi_tecalTcut",cut1);
  TH2F *h2_matchChi_ecalT = new TH2F("h2_matchChi_ecalT","matching #chi^{2} vs 2 cluster tdiff",200,-10,10,200,0,30);
  tt->Draw("eleMatchChisq:(eleClT-posClT)>>h2_matchChi_ecalT",cut1);
  TH2F *h2_momAsy_copl = new TH2F("h2_momAsy_copl","Momentum asymmetry vs coplanarity",200,0,1,200,0,220);
  tt->Draw("abs(eleP-posP)/(eleP+posP):((atan(eleTrkEcalY/(eleTrkEcalX-42.5))-atan(posTrkEcalY/(posTrkEcalX-42.5)))*180/3.14+180)>>h2_momAsy_copl",cut1);
  TH2F *h2_tecalT_copl = new TH2F("h2_tecalT_copl","Track-ecal tdiff vs coplanarity",200,0,220,200,-10,10);
  tt->Draw("(eleClT-eleTrkT-43):(atan(eleTrkEcalY/(eleTrkEcalX-42.5))-atan(posTrkEcalY/(posTrkEcalX-42.5)))*180/3.14+180>>h2_tecalT_copl",cut1);
  TH2F *h2_tecalT_ecalT = new TH2F("h2_tecalT_ecalT","Track-ecal tdiff vs 2 cluster tdiff",200,-10,10,200,-10,10);
  tt->Draw("(eleClT-eleTrkT-43):(eleClT-posClT)>>h2_tecalT_ecalT",cut1);

  //plot e+ track d0
  TH1F *h_epd0 = new TH1F("h_epd0",";e+ d0",100,-5,5);
  tt->Draw("posTrkD0>>h_epd0",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11);

  cout<<"finished all plots"<<endl;

  ///////////////////
  //plot the t1-t2 on the same plot
  TCanvas *ctime = new TCanvas("ctime","ctime",800,600);
  ctime->cd();
  h_tdiff1->SetLineColor(kRed);
  h_tdiff1->SetLineWidth(2);
  h_tdiff2->SetLineColor(kMagenta);
  h_tdiff2->SetLineWidth(2);
  h_tdiff3->SetLineColor(kBlue);
  h_tdiff3->SetLineWidth(2);
  h_tdiff4->SetLineColor(kGreen);
  h_tdiff4->SetLineWidth(2);
  h_tdiff5->SetLineColor(kCyan);
  h_tdiff5->SetLineWidth(2);
  h_tdiff6->SetLineColor(kRed+2);
  h_tdiff6->SetLineWidth(2);
  h_tdiff7->SetLineColor(kMagenta+2);
  h_tdiff7->SetLineWidth(2);
  h_tdiff8->SetLineColor(kBlue+2);
  h_tdiff8->SetLineWidth(2);
  h_tdiff9->SetLineColor(kGreen+2);
  h_tdiff9->SetLineWidth(2);
  h_tdiff11->SetLineColor(kYellow);
  h_tdiff11->SetLineWidth(2);
  h_tdiff12->SetLineColor(kYellow+2);
  h_tdiff12->SetLineWidth(2);
  h_tdiff13->SetLineColor(kCyan+4);
  h_tdiff13->SetLineWidth(2);
  TLegend *ltime = new TLegend(0.65,0.55,0.9,0.85);
  h_tdiff1->Draw();
  h_tdiff2->Draw("same");
  h_tdiff3->Draw("same");
  h_tdiff4->Draw("same");
  h_tdiff5->Draw("same");
  h_tdiff6->Draw("same");
  h_tdiff7->Draw("same");
  h_tdiff8->Draw("same");
  h_tdiff9->Draw("same");
  h_tdiff11->Draw("same");
  h_tdiff12->Draw("same");
  h_tdiff13->Draw("same");
  ltime->AddEntry("h_tdiff1","L2L2","l");
  ltime->AddEntry("h_tdiff2","+Trk chi2 cut","l");
  ltime->AddEntry("h_tdiff3","+Pmax single track","l");
  ltime->AddEntry("h_tdiff4","+Isolation cut","l");
  ltime->AddEntry("h_tdiff5","+Bsc vtx chi2 cut","l");
  ltime->AddEntry("h_tdiff6","+BscUnc chi2 diff cut","l");
  ltime->AddEntry("h_tdiff7","+Psum <1.2 cut","l");
  ltime->AddEntry("h_tdiff8","+match chi2 cut","l");
  ltime->AddEntry("h_tdiff9","+track-ecal timing cut","l");
  ltime->AddEntry("h_tdiff11","+momentum asymmetry cut","l");
  ltime->AddEntry("h_tdiff12","+e+ D0 cut","l");
  ltime->AddEntry("h_tdiff13","+shared hits<4","l");
  ltime->Draw();
  ctime->Update();
  ctime->SaveAs("../output_L2L2/tdiff_cuts.C");
  ctime->Close();

  //plot zvtx on the same plot
  TCanvas *czvtx = new TCanvas("czvtx","czvtx",800,600);
  czvtx->cd();
  h_z1->SetLineColor(kRed);
  h_z1->SetLineWidth(2);
  h_z2->SetLineColor(kMagenta);
  h_z2->SetLineWidth(2);
  h_z3->SetLineColor(kBlue);
  h_z3->SetLineWidth(2);
  h_z4->SetLineColor(kGreen);
  h_z4->SetLineWidth(2);
  h_z5->SetLineColor(kCyan);
  h_z5->SetLineWidth(2);
  h_z6->SetLineColor(kRed+2);
  h_z6->SetLineWidth(2);
  h_z7->SetLineColor(kMagenta+2);
  h_z7->SetLineWidth(2);
  h_z8->SetLineColor(kBlue+2);
  h_z8->SetLineWidth(2);
  h_z9->SetLineColor(kGreen+2);
  h_z9->SetLineWidth(2);
  h_z10->SetLineColor(kCyan+2);
  h_z10->SetLineWidth(2);
  h_z11->SetLineColor(kYellow);
  h_z11->SetLineWidth(2);
  h_z12->SetLineColor(kYellow+2);
  h_z12->SetLineWidth(2);
  h_z13->SetLineColor(kCyan+4);
  h_z13->SetLineWidth(2);
  TLegend *lzvtx = new TLegend(0.65,0.55,0.9,0.85);
  h_z1->Draw();
  h_z2->Draw("same");
  h_z3->Draw("same");
  h_z4->Draw("same");
  h_z5->Draw("same");
  h_z6->Draw("same");
  h_z7->Draw("same");
  h_z8->Draw("same");
  h_z9->Draw("same");
  h_z10->Draw("same");
  h_z11->Draw("same");
  h_z12->Draw("same");
  h_z13->Draw("same");
  lzvtx->AddEntry("h_z1","L2L2","l");
  lzvtx->AddEntry("h_z2","+Trk chi2 cut","l");
  lzvtx->AddEntry("h_z3","+Pmax single track","l");
  lzvtx->AddEntry("h_z4","+Isolation cut","l");
  lzvtx->AddEntry("h_z5","+Bsc vtx chi2 cut","l");
  lzvtx->AddEntry("h_z6","+BscUnc chi2 diff cut","l");
  lzvtx->AddEntry("h_z7","+Psum <1.2 cut","l");
  lzvtx->AddEntry("h_z8","+match chi2 cut","l");
  lzvtx->AddEntry("h_z9","+track-ecal timing cut","l");
  lzvtx->AddEntry("h_z10","+cluster time diff cut","l");
  lzvtx->AddEntry("h_z11","+momentum asymmetry","l");
  lzvtx->AddEntry("h_z12","+e+ track d0<1.5","l");
  lzvtx->AddEntry("h_z13","+shared hits<4","l");
  lzvtx->Draw();
  czvtx->Update();
  czvtx->SaveAs("../output_L2L2/zvtx_cuts.C");
  czvtx->Close();
  
  //plot the mass on the same plot
  TCanvas *cmass = new TCanvas("cmass","cmass",800,600);
  cmass->cd();
  h_m1->SetLineColor(kRed);
  h_m1->SetLineWidth(2);
  h_m2->SetLineColor(kMagenta);
  h_m2->SetLineWidth(2);
  h_m3->SetLineColor(kBlue);
  h_m3->SetLineWidth(2);
  h_m4->SetLineColor(kGreen);
  h_m4->SetLineWidth(2);
  h_m5->SetLineColor(kCyan);
  h_m5->SetLineWidth(2);
  h_m6->SetLineColor(kRed+2);
  h_m6->SetLineWidth(2);
  h_m7->SetLineColor(kMagenta+2);
  h_m7->SetLineWidth(2);
  h_m8->SetLineColor(kBlue+2);
  h_m8->SetLineWidth(2);
  h_m9->SetLineColor(kGreen+2);
  h_m9->SetLineWidth(2);
  h_m10->SetLineColor(kCyan+2);
  h_m10->SetLineWidth(2);
  h_m11->SetLineColor(kYellow);
  h_m11->SetLineWidth(2);
  h_m12->SetLineColor(kYellow+2);
  h_m12->SetLineWidth(2);
  h_m13->SetLineColor(kCyan+4);
  h_m13->SetLineWidth(2);
  
  TLegend *lmass = new TLegend(0.65,0.55,0.9,0.85);
  h_m1->Draw();
  h_m2->Draw("same");
  h_m3->Draw("same");
  h_m4->Draw("same");
  h_m5->Draw("same");
  h_m6->Draw("same");
  h_m7->Draw("same");
  h_m8->Draw("same");
  h_m9->Draw("same");
  h_m10->Draw("same");
  h_m11->Draw("same");
  h_m12->Draw("same");
  h_m13->Draw("same");
  lmass->AddEntry("h_m1","L2L2","l");
  lmass->AddEntry("h_m2","+Trk chi2 cut","l");
  lmass->AddEntry("h_m3","+Pmax single track","l");
  lmass->AddEntry("h_m4","+Isolation cut","l");
  lmass->AddEntry("h_m5","+Bsc vtx chi2 cut","l");
  lmass->AddEntry("h_m6","+BscUnc chi2 diff cut","l");
  lmass->AddEntry("h_m7","+Psum <1.2 cut","l");
  lmass->AddEntry("h_m8","+match chi2 cut","l");
  lmass->AddEntry("h_m9","+track-ecal timing cut","l");
  lmass->AddEntry("h_m10","+cluster time diff cut","l");
  lmass->AddEntry("h_m11","+momentum asymmetry","l");
  lmass->AddEntry("h_m12","+ e+ track d0<1.5","l");
  lmass->AddEntry("h_m13","+shared hits<4","l");
  lmass->Draw();
  cmass->Update();
  cmass->SaveAs("../output_L2L2/mass_cuts.C");
  cmass->Close();

  /////////////////
 
  TCanvas *vertex = new TCanvas("vertex","vertex",1200,800);
 vertex->Divide(2,2);
 vertex->cd(1);
 h_yvx->Draw("colz");
 vertex->cd(2);
 h_yvm->Draw("colz");
 vertex->cd(3);
 h_xvm->Draw("colz");
 vertex->cd(4);
 h_zvm->Draw("colz");
 vertex->Update();
 vertex->SaveAs("../output_L2L2/vertex_cuts.C");
 vertex->Close();


 ///////////////////
 //Plot each cut after all other cuts applied:
 TH1F *h_trkChi2_post = new TH1F("h_trkChi2_post","trk chi2",200,0,40);
 tt->Draw("eleTrkChisq>>h_trkChi2_post",cut1 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);
 tt->Draw("posTrkChisq>>h_trkChi2_post",cut1 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);
 TH1F *h_trkChi2_z = new TH1F("h_trkChi2_z","zvtx",200,-100,100);
 tt->Draw("uncVZ>>h_trkChi2_z",cut1 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);

 TH1F *h_eleMom_post = new TH1F("h_eleMom_post","e- momentum",200,0,1.3);
 tt->Draw("eleP>>h_eleMom_post",cut1 && cut2 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);
 TH1F *h_eleMom_z = new TH1F("h_eleMom_z","zvtx",200,-100,100);
 tt->Draw("uncVZ>>h_eleMom_z",cut1 && cut2 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);

 TH1F *h_bscChi2_post = new TH1F("h_bscChi2_post","bsc chi2",200,0,30);
 tt->Draw("bscChisq>>h_bscChi2_post",cut1 && cut2 && cut3 && cut4 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);
 TH1F *h_bscChi2_z = new TH1F("h_bscChi2_z","zvtx",200,-100,100);
 tt->Draw("uncVZ>>h_bscChi2_z",cut1 && cut2 && cut3 && cut4 && cut6 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);

 TH1F *h_bmuChi2_post = new TH1F("h_bmuChi2_post","bsc - unc chi2",200,0,20);
 tt->Draw("(bscChisq-uncChisq)>>h_bmuChi2_post",cut1 && cut2 && cut3 && cut4 && cut5 && cut7 && cut8 && cut9 && cut10 && cut11&& cut12 && cut13);
 TH1F *h_bmuChi2_z = new TH1F("h_bmuChi2_z","zvtx",200,-100,100);
 tt->Draw("uncVZ>>h_bmuChi2_z",cut1 && cut2 && cut3 && cut4 && cut5 && cut7 && cut8 && cut9 && cut10 && cut11 && cut12 && cut13 );


 TH1F *h_momAsy_post = new TH1F("h_momAsy_post","momentum asymmetry",200,0,1);
 TH2F *h_pSumDiff_post = new TH2F("h_pSumDiff_post","momentum sum vs difference",200,0,1.3,200,0,0.65);
 tt->Draw("abs(eleP-posP)/(eleP+posP)>>h_momAsy_post",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10&& cut12 && cut13);
 tt->Draw("abs(eleP-posP):uncP>>h_pSumDiff_post",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut12 && cut13 );
 TH1F *h_momAsy_z = new TH1F("h_momAsy_z","zvtx",200,-100,100);
 tt->Draw("uncVZ>>h_momAsy_z",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut9 && cut10 && cut12 && cut13 );

 TH1F *h_zvtxPeak = new TH1F("h_zvtxPeak",";z vertex +/-2ns cluster time diff",100,-100,100);
 tt->Draw("uncVZ>>h_zvtxPeak",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut11 && cut12 && cut13 &&"abs(eleClT-posClT)<1");
 TH1F *h_zvtxAcc = new TH1F("h_zvtxAcc",";z vertex >+/-2ns cluster time diff",100,-100,100);
 tt->Draw("uncVZ>>h_zvtxAcc",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut11 && cut12 && cut13 &&"abs(eleClT-posClT)>3 && abs(eleClT-posClT)<5");
 TH2F *h_zVmPeak = new TH2F("h_zVmPeak","Within +/- 2ns cluster time diff;mass [GeV];z vertex [mm]",200,0,0.1,200,-100,100);
 tt->Draw("uncVZ:uncM>>h_zVmPeak",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut11 && cut12 && cut13 &&"abs(eleClT-posClT)<2");
 TH2F *h_zVmAcc = new TH2F("h_zVmAcc","> +/- 2ns cluster time diff; mass [GeV];z vertex [mm]",200,0,0.1,200,-100,100);
 tt->Draw("uncVZ:uncM>>h_zVmAcc",cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 && cut11 && cut12 && cut13 && "abs(eleClT-posClT)>3 && abs(eleClT-posClT)<9");

 //Plot the z vertex for the missing cuts
 TCanvas *ccz = new TCanvas("ccz","z vertex for cuts",800,800);
 ccz->cd();
 h_trkChi2_z->SetLineColor(kRed);
 h_trkChi2_z->SetLineWidth(2);
 h_eleMom_z->SetLineColor(kMagenta);
 h_eleMom_z->SetLineWidth(2);
 h_bscChi2_z->SetLineColor(kBlue);
 h_bscChi2_z->SetLineWidth(2);
 h_bmuChi2_z->SetLineColor(kGreen);
 h_bmuChi2_z->SetLineWidth(2);
 h_momAsy_z->SetLineColor(kCyan);
 h_momAsy_z->SetLineWidth(2);
 TLegend *llz = new TLegend(0.65,0.55,0.9,0.85);
 h_trkChi2_z->Draw();
 h_eleMom_z->Draw("same");
 h_bscChi2_z->Draw("same");
 h_bmuChi2_z->Draw("same");
 h_momAsy_z->Draw("same");
 llz->AddEntry("h_trkChi2_z","all cuts except #chi^{2}_{trk}","l");
 llz->AddEntry("h_eleMom_z","all cuts except P_{max} for e- track","l");
 llz->AddEntry("h_bscChi2_z","all cuts except #chi^{2}_{bsc}","l");
 llz->AddEntry("h_bmuChi2_z","all cuts except #chi^{2}_{bsc}-#chi^{2}_{unc}","l");
 llz->AddEntry("h_momAsy_z","all cuts except P_{asymmetry}","l");
 llz->Draw();
 ccz->Update();
 ccz->SaveAs("../output_L2L2/zVertexCutEffects.C");
 ccz->Close();

 

 //Plot the change in number of events from cuts in the peak versus the background
 std::fstream myfile("../output_L2L2/cutEffects.txt", std::ios_base::out);
  
 Double_t ratio[11];
 Double_t cutN[11];
 Double_t totalEff[11];
 TCanvas *ccr = new TCanvas("ccr","Cut efficiencies",800,800);
 ccr->cd();
 for (int ii=1; ii<13; ii++){
   double num = ((signal[ii-1]-signal[ii])/signal[ii-1]);
   double denom = ((background[ii-1]-background[ii])/background[ii-1]);
   double totalN = signal[ii]+background[ii];
   double totalNchange = ((signal[ii-1]+background[ii-1])-totalN)/(signal[ii-1]+background[ii-1]);
   cout <<"ii:\t"<<ii<<"\t"<<totalNchange*100<<"\t"<<num*100<<"\t"<<denom*100<<endl;
   myfile <<"ii:\t"<<ii<<"\t"<<totalNchange*100<<"\t"<<num*100<<"\t"<<denom*100<<endl;
   
   ratio[ii-1] = num/denom;
   cutN[ii-1] = ii;
 }
 TGraph *eff = new TGraph(11,cutN,ratio);
 eff->SetMarkerStyle(22);
 eff->SetMarkerColor(kBlue);
 eff->Draw("ap");
 ccr->Update();
 ccr->SaveAs("../output_L2L2/cutEfficiencies.C");
 ccr->Close();


 //plot the z vtx for the stuff in and out of time
 TCanvas *cczs = new TCanvas("cczs","z vertex for cuts",800,800);
 cczs->cd();
 h_zvtxPeak->SetLineColor(kRed);
 h_zvtxPeak->SetLineWidth(2);
 h_zvtxAcc->SetLineColor(kBlue);
 h_zvtxAcc->SetLineWidth(2);
 TLegend *llzs = new TLegend(0.65,0.55,0.9,0.85);
 h_zvtxPeak->Draw();
 h_zvtxAcc->Draw("same");
 llzs->AddEntry("h_zvtxPeak","+/-2ns peak","l");
 llzs->AddEntry("h_zvtxAcc","> +/-2ns peak","l");
 llzs->Draw();
 cczs->Update();
 cczs->SaveAs("../output_L2L2/zVertexAccidentals.C");
 cczs->Close();

 
 
 fout->Write();
 fout->Close();
 
}
