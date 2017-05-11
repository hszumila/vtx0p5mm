double fitZcutUB(double mass){
  
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
  return -10.856+6413.01*mass -180028*pow(mass,2)+1439930*pow(mass,3);
 
}

void studyCutsL1L2p(){

  TFile *f = new TFile("../data/ub_0p5cut.root");
  TTree *tt=(TTree*)f->Get("ntuple");
  
  TCanvas *canvas = new TCanvas("canvas"," cut Plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  
  TFile *fout = new TFile("../output_L1L2/ub/wiggleCutsP.root","RECREATE");
  std::string pdf_file_name = "../output_L1L2/ub/wiggleCutsP.pdf";
  gROOT->SetBatch(true);


  const int nPts = 100;
  double mass[nPts];
  double zCut[nPts];
  for (int i=0; i<nPts; i++){
    mass[i] = (double) i/1000;
    zCut[i] = fitZcutUB(mass[i]);
    cout<<"mass\t"<<mass[i]<<"\tzCut\t"<<zCut[i]<<endl;
  }
  

  TGraph *gx = new TGraph(nPts,mass,zCut);
  gx->SetLineColor(kRed);
  gx->SetLineWidth(2);
  

  TCut evCut="!eleHasL1&&posHasL1&&abs(eleTrkExtrpYL1)<1.5&&eleHasL2&&posHasL2&&uncP>0.845&&isPair1&&eleClY*posClY<0&&max(eleTrkChisq,posTrkChisq)<30&&eleP<0.75*1.056&&min(eleMinPositiveIso+0.5*(eleTrkZ0-5*elePY/eleP)*sign(elePY),posMinPositiveIso+0.5*(posTrkZ0-5*posPY/posP)*sign(posPY))>0&&bscChisq<10&&bscChisq-uncChisq<5&&uncP<1.15*1.056&&max(eleMatchChisq,posMatchChisq)<10&&abs(eleClT-eleTrkT-43)<4&&abs(posClT-posTrkT-43)<4&&abs(eleClT-posClT)<2&&abs(eleP-posP)/(eleP+posP)<0.4&&posTrkD0<1.5&&(posMaxHitsShared<5&&eleMaxHitsShared<5)";

  //plot the initial distribution
  TH2F *h_initial = new TH2F("h_initial","No new cuts from 10%;corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_initial",evCut);
  cout<<"Finished initial plot"<<endl;

  //plot the layer 2 iso cut
  TH2F *h_iso = new TH2F("h_iso","Applied L2 iso cut;corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  TCut isoCut = "min(eleMinPositiveIsoL2+0.33*(eleTrkZ0-5*elePY/eleP)*sign(elePY),posMinPositiveIsoL2+0.33*(posTrkZ0-5*posPY/posP)*sign(posPY))>0";
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_iso",evCut && isoCut);
  cout<<"Finished L2 iso cut"<<endl;
  
  //show what remains if using diagonal trk chi2 cut
  TH2F *h_trkchi = new TH2F("h_trkchi","Applied trk chi2 diagonal cut and L2 iso cut;corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  TCut trkChiCut = "eleTrkChisq+posTrkChisq<30";
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_trkchi",evCut && trkChiCut && isoCut);
  cout<<"Finished trkChisq cut"<<endl;

  //show what remains if using kink cuts
  TH2F *h_kinkS = new TH2F("h_kinkS","Applied kink cuts square and L2 iso cut;corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  TCut kinkSCut = "abs(elePhiKink2)<0.002&&abs(posPhiKink2)<0.002&&abs(elePhiKink1)<0.0001&&abs(posPhiKink1)<0.0001&&abs(eleLambdaKink2)<0.004&&abs(posLambdaKink2)<0.004&&abs(elePhiKink3)<0.002&&abs(posPhiKink3)<0.002&&abs(posLambdaKink3)<0.004&&abs(eleLambdaKink3)<0.004";
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_kinkS",evCut && kinkSCut && isoCut);
  cout<<"Finished square kink cut"<<endl;
  
  TH2F *h_kinkE = new TH2F("h_kinkE","Applied kink cuts elliptical and L2 iso cut;corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  TCut kinkECut = "(pow((elePhiKink2*-0.98-eleLambdaKink2*0.2)/0.002,2)+pow((eleLambdaKink2*-0.98+elePhiKink2*0.2)/0.0045,2))<1&&(pow((posPhiKink2*-0.98-posLambdaKink2*0.2)/0.002,2)+pow((posLambdaKink2*-0.98+posPhiKink2*0.2)/0.0045,2))<1&&(pow((elePhiKink3*0.98-eleLambdaKink3*0.2)/0.002,2)+pow((eleLambdaKink3*0.98+elePhiKink3*0.2)/0.0045,2))<1&&(pow((posPhiKink3*0.98-posLambdaKink3*0.2)/0.002,2)+pow((posLambdaKink3*0.98+posPhiKink3*0.2)/0.0045,2))<1&&abs(elePhiKink1)<0.0001&&abs(posPhiKink1)<0.0001&&abs(eleLambdaKink1)<0.0013&&abs(posLambdaKink1)<0.0013";
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_kinkE",evCut && kinkECut && isoCut);  
  cout<<"Finished elliptical kink cut"<<endl;
  
  //show what remains if using diagonal matching cut
  TH2F *h_match = new TH2F("h_match","Applied matching diagonal cut and L2 iso cut;corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  TCut matchCut = "eleMatchChisq+posMatchChisq<10";
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_match",evCut && matchCut && isoCut);
  cout<<"Finished matching cut"<<endl;
  
  //show sum of all cuts
  TH2F *h_all = new TH2F("h_all","All new cuts (square kink);corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  TH2F *h_allE = new TH2F("h_allE","All new cuts (elliptical kink);corrected mass [GeV];unconstrained z vertex [mm]",200,0,0.1,200,-70,70);
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_all",evCut && trkChiCut && kinkSCut && matchCut && isoCut);
  tt->Draw("uncVZ:uncM - 0.15e-3 * (elePX/eleP - posPX/posP) * uncVZ/uncM>>h_allE",evCut && trkChiCut && kinkECut && matchCut && isoCut);
  cout<<"Finished all plots, now printing"<<endl;

  /////////////////////////////////
  canvas->Update();
  h_initial->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_iso->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_trkchi->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_kinkS->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + "(").c_str());
 
  h_kinkE->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + "(").c_str());
 
  h_match->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + "(").c_str());
 
  h_all->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + "(").c_str());
 
  h_allE->Draw("colz");
  gx->Draw("lsame");
  canvas->Print( (pdf_file_name + ")").c_str());

  fout->Write();
  fout->Close();

 
}
