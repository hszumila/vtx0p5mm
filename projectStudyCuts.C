void projectStudyCuts(){

  //read the input file
  TFile *f = new TFile("../output_L1L1/ub/wiggleCuts.root");

  //read the input plots
  TH2F *h2_initial = new TH2F("h2_initial","zVertex vs Mass",200,0.0,0.1,200,-60,60);
  TH2F *h2_iso = new TH2F("h2_iso","zVertex vs Mass",200,0.0,0.1,200,-60,60);
  TH2F *h2_trkchi = new TH2F("h2_trkchi","zVertex vs Mass",200,0.0,0.1,200,-60,60);
  TH2F *h2_kinkS = new TH2F("h2_kinkS","zVertex vs Mass",200,0.0,0.1,200,-60,60);
  TH2F *h2_kinkE = new TH2F("h2_kinkE","zVertex vs Mass",200,0.0,0.1,200,-60,60);
  TH2F *h2_match = new TH2F("h2_match","zVertex vs Mass",200,0.0,0.1,200,-60,60);
  TH2F *h2_all = new TH2F("h2_all","zVertex vs Mass",200,0.0,0.1,200,-60,60);
  TH2F *h2_allE = new TH2F("h2_allE","zVertex vs Mass",200,0.0,0.1,200,-60,60);

  h2_initial = (TH2F*)f->Get("h_initial");
  h2_iso = (TH2F*)f->Get("h_iso");
  h2_trkchi = (TH2F*)f->Get("h_trkchi");
  h2_kinkS = (TH2F*)f->Get("h_kinkS");
  h2_kinkE = (TH2F*)f->Get("h_kinkE");
  h2_match = (TH2F*)f->Get("h_match");
  h2_all = (TH2F*)f->Get("h_all");
  h2_allE = (TH2F*)f->Get("h_allE");
    
    //make the output plots
    TH1D *h1Initial = new TH1D("h1Initial",";z vertex [mm]",200,-60,60);
    TH1D *h1Iso = new TH1D("h1Iso",";z vertex [mm]",200,-60,60);
    TH1D *h1Trkchi = new TH1D("h1Trkchi",";z vertex [mm]",200,-60,60);
    TH1D *h1KinkS = new TH1D("h1KinkS",";z vertex [mm]",200,-60,60);
    TH1D *h1KinkE = new TH1D("h1KinkE",";z vertex [mm]",200,-60,60);
    TH1D *h1Match = new TH1D("h1Match",";z vertex [mm]",200,-60,60);
    TH1D *h1All = new TH1D("h1All",";z vertex [mm]",200,-60,60);
    TH1D *h1AllE = new TH1D("h1AllE",";z vertex [mm]",200,-60,60);

    //Project the 2d plots
    h1Initial = (TH1D*)h2_initial->ProjectionY("h1Initial",0,200);
    h1Iso = (TH1D*)h2_iso->ProjectionY("h1Iso",0,200);
    h1Trkchi = (TH1D*)h2_trkchi->ProjectionY("h1Trkchi",0,200);
    h1KinkS = (TH1D*)h2_kinkS->ProjectionY("h1KinkS",0,200);
    h1KinkE = (TH1D*)h2_kinkE->ProjectionY("h1KinkE",0,200);
    h1Match = (TH1D*)h2_match->ProjectionY("h1Match",0,200);
    h1All = (TH1D*)h2_all->ProjectionY("h1All",0,200);
    h1AllE = (TH1D*)h2_allE->ProjectionY("h1AllE",0,200);

    //Superimpose and plot
    TCanvas *cc = new TCanvas("cc","cc",800,600);
    cc->cd();
    h1Initial->SetLineColor(kRed);
    h1Initial->SetLineWidth(2);
    h1Iso->SetLineColor(kMagenta);
    h1Iso->SetLineWidth(2);
    h1Trkchi->SetLineColor(kBlue);
    h1Trkchi->SetLineWidth(2);
    h1KinkS->SetLineColor(kGreen);
    h1KinkS->SetLineWidth(2);
    h1KinkE->SetLineColor(kCyan);
    h1KinkE->SetLineWidth(2);
    h1Match->SetLineColor(kCyan+4);
    h1Match->SetLineWidth(2);
    h1All->SetLineColor(kRed+2);
    h1All->SetLineWidth(2);
    h1AllE->SetLineColor(kMagenta+2);
    h1AllE->SetLineWidth(2);
    TLegend *leg = new TLegend(0.65,0.55,0.9,0.85);
    h1Initial->Draw();
    h1Iso->Draw("same");
    h1Trkchi->Draw("same");
    h1KinkS->Draw("same");
    h1KinkE->Draw("same");
    h1Match->Draw("same");
    h1All->Draw("same");
    h1AllE->Draw("same");
    leg->AddEntry("h1Initial","Initial cuts form 10%","l");
    leg->AddEntry("h1Iso","+L2 iso cut","l");
    leg->AddEntry("h1Trkchi","diagonal trk chi2 cut","l");
    leg->AddEntry("h1KinkS","square kink cuts","l");
    leg->AddEntry("h1KinkE","elliptical kink cuts","l");
    leg->AddEntry("h1Match","diagonal matching cuts","l");
    leg->AddEntry("h1All","all cuts (square kink)","l");
    leg->AddEntry("h1AllE","all cuts (elliptical kink)","l");
    leg->Draw();
    cc->Update();
    cc->SaveAs("../output_L1L1/ub/proectedCuts.C");
    cc->Close();

}
