

void plotKinkCutsL1L2(){

 
  

/*This file reads the slic reconstructed output files and makes an ntuple that
can be used for FEE calibration. 
*/

/////////////////////
  
  TFile *f = new TFile("../data/hps_L1L2_blind_10Feb17.root");
  TTree *tt=(TTree*)f->Get("ntuple");

  TCanvas *canvas = new TCanvas("canvas","Kink Plots", 700, 700);
 canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  
  TFile *fout = new TFile("../output_L1L2/kinkPlotsL1L2_noKinkCut.root","RECREATE");
  std::string pdf_file_name = "../output_L1L2/kinkPlotsL1L2_noKinkCut.pdf";
  gROOT->SetBatch(true);


  //Setup plots:
  ////////////////////
  /////TrackStuff/////
  ////////////////////
  cout<<"make plots"<<endl;
 
  //plot kink
  TH1F *hm_phiKink2 = new TH1F("hm_phiKink2","e- phi Kink L2;e- phi kink L2",100,-0.005,0.005);
  TH1F *hp_phiKink2 = new TH1F("hp_phiKink2","e+ phi Kink L2;e+ phi kink L2",100,-0.005,0.005);
  TH1F *hm_phiKink1 = new TH1F("hm_phiKink1","e- phi Kink L1;e- phi kink L1",100,-0.0005,0.0005);
  TH1F *hp_phiKink1 = new TH1F("hp_phiKink1","e+ phi Kink L1;e+ phi kink L1",100,-0.0005,0.0005);

  TH1F *hm_LKink2 = new TH1F("hm_LKink2","e- lambda Kink L2;e- lambda kink L2",100,-0.01,0.01);
  TH1F *hp_LKink2 = new TH1F("hp_LKink2","e+ lambda Kink L2;e+ lambda kink L2",100,-0.01,0.01);
  TH1F *hm_LKink1 = new TH1F("hm_LKink1","e- lambda Kink L1;e- lambda kink L1",100,-0.005,0.005);
  TH1F *hp_LKink1 = new TH1F("hp_LKink1","e+ lambda Kink L1;e+ lambda kink L1",100,-0.005,0.005);


  TH2F *hm_LPkink1 = new TH2F("hm_LPkink1","  e- lambda vs phi kink L1;lambda kink L1;phi kink L1",100,-0.005,0.005,100,-0.0005,0.0005);
  TH2F *hm_LPkink2 = new TH2F("hm_LPkink2","  e- lambda vs phi kink L2;lambda kink L2;phi kink L2",100,-0.01,0.01,100,-0.005,0.005);
  TH2F *hm_Lkink12 = new TH2F("hm_Lkink12","  e- lambda kink;lambda kink L1;lambda kink L2",100,-0.005,0.005,100,-0.01,0.01);
  TH2F *hm_Pkink12 = new TH2F("hm_Pkink12","  e- phi kink;phi kink L1;phi kink L2",100,-0.0005,0.0005,100,-0.005,0.005);

  TH2F *hp_LPkink1 = new TH2F("hp_LPkink1","  e+ lambda vs phi kink L1;lambda kink L1;phi kink L1",100,-0.005,0.005,100,-0.0005,0.0005);
  TH2F *hp_LPkink2 = new TH2F("hp_LPkink2","  e+ lambda vs phi kink L2;lambda kink L2;phi kink L2",100,-0.01,0.01,100,-0.005,0.005);
  TH2F *hp_Lkink12 = new TH2F("hp_Lkink12","  e+ lambda kink;lambda kink L1;lambda kink L2",100,-0.005,0.005,100,-0.01,0.01);
  TH2F *hp_Pkink12 = new TH2F("hp_Pkink12","  e+ phi kink;phi kink L1;phi kink L2",100,-0.0005,0.0005,100,-0.005,0.005);

 //plot kink
  TH1F *hm_phiKink3 = new TH1F("hm_phiKink3","e- phi Kink L3;e- phi kink L3",100,-0.005,0.005);
  TH1F *hp_phiKink3 = new TH1F("hp_phiKink3","e+ phi Kink L3;e+ phi kink L3",100,-0.005,0.005);

  TH1F *hm_LKink3 = new TH1F("hm_LKink3","e- lambda Kink L3;e- lambda kink L3",100,-0.01,0.01);
  TH1F *hp_LKink3 = new TH1F("hp_LKink3","e+ lambda Kink L3;e+ lambda kink L3",100,-0.01,0.01);

  TH2F *hm_LPkink3 = new TH2F("hm_LPkink3","  e- lambda vs phi kink L3;lambda kink L2;phi kink L3",100,-0.01,0.01,100,-0.005,0.005);

  TH2F *hp_LPkink3 = new TH2F("hp_LPkink3","  e+ lambda vs phi kink L3;lambda kink L2;phi kink L3",100,-0.01,0.01,100,-0.005,0.005);
  


  
   /////////////////////////////////////////////////////////////////////////////////////////////
   TCut evCut= "((eleHasL1&&!posHasL1&&abs(posFirstHitY-100*posPY/posP)<1.5)||(!eleHasL1&&posHasL1&&abs(eleFirstHitY-100*elePY/eleP)<1.5))&&posHasL2&&eleHasL2&&uncP>0.845&&isPair1&&max(eleTrkChisq,posTrkChisq)<30&&eleP<0.75*1.056&&min(eleMinPositiveIso+0.5*(eleTrkZ0-5*elePY/eleP)*sign(elePY),posMinPositiveIso+0.5*(posTrkZ0-5*posPY/posP)*sign(posPY))>0&&bscChisq<10&&bscChisq-uncChisq<4&&uncP<1.15*1.056&&max(eleMatchChisq,posMatchChisq)<10&&abs(eleClT-eleTrkT-43)<4&&abs(posClT-posTrkT-43)<4&&abs(eleClT-posClT)<2&&posTrkD0<1.5&&abs(eleP-posP)/(eleP+posP)<0.4&&eleMaxHitsShared<5&&posMaxHitsShared<5";
   
   int ii=0;
   TCut signal = "uncVZ>-5";
   
     tt->Draw("elePhiKink2>>hm_phiKink2",signal && evCut);
     tt->Draw("posPhiKink2>>hp_phiKink2",signal && evCut);
     tt->Draw("elePhiKink1>>hm_phiKink1",signal && evCut);
     tt->Draw("posPhiKink1>>hp_phiKink1",signal && evCut);
     tt->Draw("eleLambdaKink2>>hm_LKink2",signal && evCut);
     tt->Draw("posLambdaKink2>>hp_LKink2",signal && evCut);
     tt->Draw("eleLambdaKink1>>hm_LKink1",signal && evCut);
     tt->Draw("posLambdaKink1>>hp_LKink1",signal && evCut);
     tt->Draw("elePhiKink1:eleLambdaKink1>>hm_LPkink1",signal && evCut);
     tt->Draw("elePhiKink2:eleLambdaKink2>>hm_LPkink2",signal && evCut);
     tt->Draw("eleLambdaKink2:eleLambdaKink1>>hm_Lkink12",signal && evCut);
     tt->Draw("elePhiKink2:elePhiKink1>>hm_Pkink12",signal && evCut);
     tt->Draw("posPhiKink1:posLambdaKink1>>hp_LPkink1",signal && evCut);
     tt->Draw("posPhiKink2:posLambdaKink2>>hp_LPkink2",signal && evCut);
     tt->Draw("posLambdaKink2:posLambdaKink1>>hp_Lkink12",signal && evCut);
     tt->Draw("posPhiKink2:posPhiKink1>>hp_Pkink12",signal && evCut);    
     tt->Draw("elePhiKink3>>hm_phiKink3",signal && evCut);
     tt->Draw("posPhiKink3>>hp_phiKink3", signal && evCut);
     tt->Draw("eleLambdaKink3>>hm_LKink3",signal && evCut);
     tt->Draw("posLambdaKink3>>hp_LKink3",signal && evCut);
     tt->Draw("elePhiKink3:eleLambdaKink3>>hm_LPkink3",signal && evCut);
     tt->Draw("posPhiKink3:posLambdaKink3>>hp_LPkink3",signal && evCut);

    


  canvas->Update();

  hm_phiKink2->SetLineColor(kBlue);
  hm_phiKink2->SetLineWidth(2);
  hm_phiKink2->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hp_phiKink2->SetLineColor(kBlue);
  hp_phiKink2->SetLineWidth(2);
  hp_phiKink2->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hm_phiKink1->SetLineColor(kBlue);
  hm_phiKink2->SetLineWidth(2);
  hm_phiKink1->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hp_phiKink1->SetLineColor(kBlue);
  hp_phiKink1->SetLineWidth(2);
  hp_phiKink1->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hm_LKink2->SetLineColor(kBlue);
  hm_LKink2->SetLineWidth(2);
  hm_LKink2->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hp_LKink2->SetLineColor(kBlue);
  hp_LKink2->SetLineWidth(2);
  hp_LKink2->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hm_LKink1->SetLineColor(kBlue);
  hm_LKink1->SetLineWidth(2);
  hm_LKink1->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hp_LKink1->SetLineColor(kBlue);
  hp_LKink1->SetLineWidth(2);
  hp_LKink1->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hm_LPkink1->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());


  hm_LPkink2->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());


  hm_Lkink12->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());


  hm_Pkink12->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());


  hp_LPkink1->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());


  hp_LPkink2->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());


  hp_Lkink12->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());


  hp_Pkink12->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  
   hm_phiKink3->SetLineColor(kBlue);
   hm_phiKink3->SetLineWidth(2);
  hm_phiKink3->Draw();
    canvas->Print( (pdf_file_name + "(").c_str());


   hp_phiKink3->SetLineColor(kBlue);
   hp_phiKink3->SetLineWidth(2);
  hp_phiKink3->Draw();
    canvas->Print( (pdf_file_name + "(").c_str());


   hm_LKink3->SetLineColor(kBlue);
   hm_LKink3->SetLineWidth(2);
  hm_LKink3->Draw();
    canvas->Print( (pdf_file_name + "(").c_str());


   hp_LKink3->SetLineColor(kBlue);
   hp_LKink3->SetLineWidth(2);
  hp_LKink3->Draw();
    canvas->Print( (pdf_file_name + "(").c_str());

  hm_LPkink3->Draw("colz");
    canvas->Print( (pdf_file_name + "(").c_str());

  hp_LPkink3->Draw("colz");
    canvas->Print( (pdf_file_name + ")").c_str());


 fout->Write();
 fout->Close();






 
}//end main







