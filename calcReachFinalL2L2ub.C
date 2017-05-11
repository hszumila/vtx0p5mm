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

/////////////////////////////////////////////////////////
//Stuff I need to run this code://///////////////////////
////////////////////////////////////////////////////////

//Mass resolution:
double mresA = 0.1565;
double mresB = -0.00177;


const double epsBin = 0.838;
const double bWidth = 1.4*2;

//mass in GeV
double fitZcut(double mass){
  //return (15.86 + 1857*mass - 62743.2*mass*mass + 536084*mass*mass*mass-5.0);
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
  return 0.0;
 
}


//fitting function
Double_t fitvtx(Double_t *x,Double_t *par){
  return par[0]*exp(((x[0]-par[1])<par[3])*(-0.5*pow((x[0]-par[1]),2)/pow(par[2],2))+((x[0]-par[1])>=par[3])*(-0.5*pow(par[3],2)/pow(par[2],2)-(x[0]-par[1]-par[3])/par[4]));
}


const double frad = 0.1;
const double alpha = 1./137;
const double zMax = 100; //mm to first layer


//E0 [GeV], mA [MeV]
//returns [mm]
double gct(double E0, double eps2, double mA){
  return 8*(E0/10)*(pow(0.0001,2)/eps2)*pow(100/mA,2);
}

// Crystal ball peak function
Double_t CrystalBall(double z, double mass) {
  double par0 = 1.0;
  double par1 = 141.5;
  double par2 = -35.922+4.80304*mass-0.045296*pow(mass,2);
  double par3 = 7.89029+0.0669*mass;
  double par4 = -0.70241+0.0458932*mass-0.000525039*pow(mass,2);
    
  Double_t sigma = par3;
  if (par3<0){sigma = -par3;}
  Double_t t = (z-par2)/sigma;
  if (par0 < 0) t = -t;
  Double_t absAlpha = fabs(par0);
  if (t >= -absAlpha) {
    return par4*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par1/absAlpha,par1)*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par1/absAlpha - absAlpha; 

    return par4*(a/TMath::Power(b - t, par1));
  }
}



double vertEff(double gct, double zz, double mA){
  //     return 1./gct*exp((-5-zz)/gct)*(-0.03242+0.00514*mA)*exp(-pow( ( zz-(12.4881+1.6875*mA))/(2*(4.7982+0.11158*mA)) ,2) );

  if (mA<22){mA=22;}
  else if (mA>60){mA=60;}
  return 1./gct*exp((-5-zz)/gct)*CrystalBall(zz,mA);
  
}


//mA in MeV, gct in mm, zCut in mm
double integrateZ(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEff(gct,iz,mA)+4*vertEff(gct,(2*iz+step)/2,mA)+vertEff(gct,iz+step,mA));
  }
  return integral;
}



double calcSignal(double E0, double eps2, double mA, double zCut, double Nbin, double mrange){
  double gg = gct(E0,eps2,mA);
  if (Nbin !=0){
    return frad*Nbin*(3*3.14*eps2/(2*alpha))*(mA/1000/mrange)*epsBin*integrateZ(gg,mA,zCut);
  }
  else {return 0;}
}

/*
 * This is the code which reads in the 2D histogram, slices it, and plots the reach.
 */

//Read in the 2d histo
int main(){
  
  TFile *f = new TFile("../output_L2L2/ub/cutPlots.root");
  TFile *f1 = new TFile("../finalZvMub.root");
  TFile *fout = new TFile("../output_L2L2/ub/outputFits.root","RECREATE");
  gROOT->SetBatch(true);

  //get z vs mass plot
  TH2F *h2 = new TH2F("h2","zVertex vs Mass",400,0,0.1,400,-100,100);
  h2 = (TH2F*)f->Get("h_zvm");


  TH2F *hc = new TH2F("hc","combined z vs m",200,0,0.1,200,-100,100);
  hc = (TH2F*)f1->Get("h4");
  hc->SetName("hc");

  const int binning = 8;
  const int nBinsT = h2->GetXaxis()->GetNbins()-2;
  TCanvas *mslice[nBinsT];  
  Double_t centermass[nBinsT];
  double massrange[nBinsT];
  float reslimited[nBinsT];
  TH1 *v1[nBinsT];
  TF1 *c1[nBinsT];
  Double_t zmean[nBinsT];
  Double_t zsigma[nBinsT];
  Double_t zconst[nBinsT];
  Double_t zCut[nBinsT];
  Double_t zCutScaled[nBinsT];
  Double_t nEventsBin[nBinsT];
  const int neps = 400;
  double yield[neps][nBinsT];
  TH2D *h_reach = new TH2D("h_reach","Reach plot",nBinsT,0,0.1,neps,1E-10,4E-8);
  Double_t par3[nBinsT];
  Double_t par4[nBinsT];
  Double_t zCutMax[nBinsT];
  const double quantile = 0.5;
  

  TCanvas *canvas = new TCanvas("canvas"," vertex slice plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLogy();
  std::string pdf_file_name = "../output_L2L2/ub/vertexSliceFits.pdf";

  
  ////////////////////////////////////////
  //slice the 2d histo by mass resolution
  ////////////////////////////////////////
  canvas->Update();
  for (int ii=0; ii<nBinsT; ii++){
    
    ///////////////////////////////////////////////////////////////////////////////
  
    centermass[ii] = h2->GetXaxis()->GetBinCenter(ii);
    double sigM = mresA*centermass[ii]+mresB;
    reslimited[ii] = bWidth*sigM;
    int lowBin = h2->GetXaxis()->FindBin(centermass[ii]-reslimited[ii]/2);
    int highBin = h2->GetXaxis()->FindBin(centermass[ii]+reslimited[ii]/2);
    double lowedge = h2->GetXaxis()->GetBinLowEdge(lowBin);
    double highedge = h2->GetXaxis()->GetBinUpEdge(highBin);
    massrange[ii] = highedge-lowedge;
     
    ////////////////////////////////////////////////////////////////////////////////
    v1[ii] = (TH1*)h2->ProjectionY(Form("slice_%d",ii),lowBin,highBin);
    v1[ii]->GetXaxis()->SetRangeUser(-50,50);
    v1[ii]->SetTitle(Form("Z Vtx, Mass[%.4f, %.4f],%.4f",lowedge, highedge, centermass[ii]));
    
    v1[ii]->Draw();    
    c1[ii] = new TF1(Form("c%d",ii),"gaus",lowedge,highedge);
    c1[ii]->SetLineColor(kViolet);
    v1[ii]->Fit(c1[ii],"QS");

    zconst[ii] = c1[ii]->GetParameter(0);
    zmean[ii] = c1[ii]->GetParameter(1);
    zsigma[ii] = c1[ii]->GetParameter(2);

    //fit each slice for zCut
    TF1 *vFit = new TF1("vFit",fitvtx,zmean[ii]-2*zsigma[ii],zmean[ii]+10*zsigma[ii],5);
    vFit->SetLineColor(kRed);
    vFit->SetParameters(zconst[ii],zmean[ii],zsigma[ii],3*zsigma[ii],5);
    vFit->SetParLimits(4,0,10);
    //LSQM
    v1[ii]->Fit(vFit,"+QMR","",zmean[ii]-2*zsigma[ii],zmean[ii]+10*zsigma[ii]);//+QR
    c1[ii]->Draw("same");
    vFit->Draw("lsame");
    if(v1[ii]->GetEntries()>100){
      canvas->Print( (pdf_file_name + "(").c_str());
    }
    par3[ii] = vFit->GetParameter(3);
    par4[ii] = vFit->GetParameter(4);
    
    //nEventsBin[ii] = v1[ii]->Integral(v1[ii]->FindBin(vFit->GetParameter(1)),v1[ii]->FindLastBinAbove(0));
    nEventsBin[ii] = hc->Integral(lowBin,highBin,0,200);
    zCut[ii] = fitZcut(centermass[ii]);
    //zCut[ii] = vFit->GetX(0.5*massrange[ii]/reslimited[ii]/vFit->GetParameter(4),zmean[ii],v1[ii]->GetBinCenter(v1[ii]->FindLastBinAbove(0)));
    zCutMax[ii] = zCut[ii]-TMath::Log(1.-quantile);//*5.2;//vFit->GetParameter(4);
    zCutScaled[ii] = zCut[ii]-TMath::Log(0.1);

    gStyle->SetOptFit(111);
       
    //calculate the signal yield at each mass
    for (int ll=0; ll<neps; ll++){
      double ieps = h_reach->GetYaxis()->GetBinCenter(ll);
      double imass = h_reach->GetXaxis()->GetBinCenter(ii);
      //   yield[ll][ii] = calcSignal(1.05, ll*pow(10,-10), centermass[ii]*1000, zCut[ii], nEventsBin[ii], massrange[ii]);
      if (imass>0.012){
	yield[ll][ii] = calcSignal(1.05, ieps, imass*1000, zCutScaled[ii], nEventsBin[ii], massrange[ii]);
      //cout<<"mass\t"<<imass<<"\teps2\t"<<ieps<<"\tyield\t"<<yield[ll][ii]<<endl;
      }
      else {yield[ll][ii]=0.0;}
      h_reach->SetBinContent(ii,ll,yield[ll][ii]);
    }//end loop neps   
  }//end loop mass
  canvas->Print( (pdf_file_name + ")").c_str());
	   
  
  //plot zVm with the zcut shown
  const int nPts = 160;
  const int nPts2 = 40;
  double cmassX[nPts];
  double cmassX2[nPts2];
  double zcutX[nPts];
  double zcutmaxX[nPts];
  double zcutSX[nPts];
  double zcutX2[nPts2];
  int jk=0;
  for (int ijk=0; ijk<nPts; ijk++){
    cmassX[ijk] = centermass[ijk+80];
    zcutX[ijk] = zCut[ijk+80];
    zcutmaxX[ijk] = zCutMax[ijk+80];
    zcutSX[ijk] = zCutScaled[ijk+80];
    if (ijk%4==0){
      cmassX2[jk] = centermass[ijk+80];
      zcutX2[jk] = zCut[ijk+80];
      jk++;
    }

  }

  //plot the zcut values
  TCanvas *ccz = new TCanvas("ccz","zcut",1000,900);
  ccz->cd();
  gStyle->SetOptStat(0);
  TGraph *zzcut = new TGraph();
  double prevCut = 0;
  int ic = 0;
  for (int jj=0; jj<nBinsT; jj++){
    if (zCut[jj]!=prevCut){
      zzcut->SetPoint(ic,centermass[jj],zCut[jj]);
      prevCut=zCut[jj];
      ic++;
    }
  }
  zzcut->Draw("alp");
  ccz->Update();
  ccz->SaveAs("../output_L2L2/ub/zcut.C");
  ccz->Close();

  
  TCanvas *kk1 = new TCanvas("kk1","z slice mass",800,800);
  kk1->cd();
  h2->Draw("colz");
  TGraph *gx = new TGraph(nPts,cmassX,zcutX);//nBinsT,centermass,zCut);
  TGraph *gxMax = new TGraph(nPts,cmassX,zcutmaxX);
  TGraph *gxS = new TGraph(nPts,cmassX,zcutSX);
  gx->SetLineColor(kRed);
  gx->SetLineWidth(2);
  gx->Draw("lsame");
  gxMax->SetLineColor(kRed);
  gxMax->SetLineStyle(2);
  gxMax->SetLineWidth(2);
  gxMax->Draw("lsame");
  gxS->SetLineColor(kMagenta+3);
  gxS->SetLineWidth(2);
  gxS->Draw("lsame");
  kk1->Update();
  kk1->SaveAs("../output_L2L2/ub/zVm_zCut.C");
  kk1->Close();

  TCanvas *kk2 = new TCanvas("kk2","par3 for mass",800,800);
  kk2->cd();
  TGraph *gx2 = new TGraph(nBinsT,centermass,par3);
  gx2->Draw("");
  kk2->Update();
  kk2->SaveAs("../output_L2L2/ub/par3_zCut.C");
  kk2->Close();

  TCanvas *kk3 = new TCanvas("kk3","par4 for mass",800,800);
  kk3->cd();
  TGraph *gx3 = new TGraph(nBinsT,centermass,par4);
  gx3->Draw("");
  kk3->Update();
  kk3->SaveAs("../output_L2L2/ub/par4_zCut.C");
  kk3->Close();
  
 
  /////////////////////////////////////
  //Plot reach////////////////////////			   
  TCanvas *cc = new TCanvas("cc","Reach",1000,900);
  cc->SetLogy();
  cc->SetLogx();
  cc->cd();
  gStyle->SetOptStat(0);
  h_reach->GetZaxis()->SetRangeUser(0.0,0.25);
  h_reach->Draw("colz");
  h_reach->GetXaxis()->SetRangeUser(0.001,0.1);
  h_reach->GetXaxis()->SetTitle("Mass [GeV]");
  h_reach->GetXaxis()->SetTitleOffset(1.3);
  h_reach->GetYaxis()->SetTitle("#epsilon^{2}");
  cc->Update();
  cc->SaveAs("../output_L2L2/ub/reach.C");
  cc->Close();  

  fout->Write();
  fout->Close(); 

}//end main



