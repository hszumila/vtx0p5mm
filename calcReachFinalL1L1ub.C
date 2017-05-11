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
//Sho: 0.0382m+0.00072
//Me: 0.034m+0.00082
double mresA = 0.0243;//0.03;
double mresB = 0.0008;


const double epsBin = 0.838;
const double bWidth = 1.4*2;

//mass in GeV
double fitZcut(double mass){
  
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
 
  return 22.4+896.463*mass-29940.4*pow(mass,2)+213639*pow(mass,3);
 
}

double fitZcutUB(double mass){
  
  if (mass<0.018){
    mass = 0.018;
    return -27.2+4462.8*mass-84970.9*pow(mass,2);//38.3+270.7*mass-20481.1*pow(mass,2)+172067*pow(mass,3);
  }
  else if (mass>0.018 && mass < 0.024){
    return -27.2+4462.8*mass-84970.9*pow(mass,2);
  }
  else if (mass > 0.024 && mass < 0.07){
    return 38.3+270.7*mass-20481.1*pow(mass,2)+172067*pow(mass,3);
  }
  else if (mass>0.07){
    mass = 0.07;
    return 38.3+270.7*mass-20481.1*pow(mass,2)+172067*pow(mass,3);
  } 
}


//fitting function
Double_t fitvtx(Double_t *x,Double_t *par){
  return par[0]*exp(((x[0]-par[1])<par[3])*(-0.5*pow((x[0]-par[1]),2)/pow(par[2],2))+((x[0]-par[1])>=par[3])*(-0.5*pow(par[3],2)/pow(par[2],2)-(x[0]-par[1]-par[3])/par[4]));
}


const double frad = 0.1;
const double alpha = 1./137;
const double zMax = 100; //mm to first layer


double radFrac(double mA){
  return 0.1;////4*(mA/1000)+0.02;
  //0.1

}

//E0 [GeV], mA [MeV]
//returns [mm]
double gct(double E0, double eps2, double mA){
  return 8*(E0/10)*(pow(0.0001,2)/eps2)*pow(100/mA,2);
}


double vertEff(double gct, double zz, double mA){
    return 1./gct*exp((-5-zz)/gct)*exp((-1.08E-5+8.1E-7*mA-1.392E-8*mA*mA)*pow(zz,3)
				       +(0.001578-0.000113*mA+1.567E-6*mA*mA)*pow(zz,2)
				      +(-0.066555+0.002836*mA-3.1037E-5*mA*mA)*zz
				       +(-0.316+0.0135*mA-0.000154*mA*mA));
  
}

double vertEffMRS(double gct, double zz, double mA){
  mA = mA/1000;
  /*
  return 1./gct*exp((-5-zz)/gct)*exp(
				     (-5.684E-6+0.0002262*mA-0.003204*mA*mA+0.01531*mA*mA*mA)*pow(zz,4)
				     +(0.0001985-0.006858*mA+0.08334*mA*mA-0.2806*mA*mA*mA)*pow(zz,3)
				     +(-0.004715+0.1889*mA-2.627*mA*mA+8.643*mA*mA*mA)*pow(zz,2)
				     +(-0.02017+0.3605*mA)*zz
				     +(-0.0737+1.918*mA)
				     );
  */
  /*
  return 1./gct*exp((-5-zz)/gct)*exp(
				     (-7.047E-5+0.001008*mA)*pow(zz,3)
				     +(0.002708-0.0446*mA)*pow(zz,2)
				     +(-0.03599+0.5486*mA)*zz
				     +(-0.2564+4.003*mA)
  				     );
  */

   return 1./gct*exp((-5-zz)/gct)*exp(
				     (-0.0002078+0.008753*mA-0.1396*mA*mA+0.8077*mA*mA*mA)*pow(zz,3)
				     +(-0.001201+0.1404*mA-2.614*mA*mA+10.65*mA*mA*mA)*pow(zz,2)
				     +(-0.03537+0.5395*mA)*zz
				     +(-0.2359+3.606*mA)
  				     );
}


//mA in MeV, gct in mm, zCut in mm
double integrateZ(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffMRS(gct,iz,mA)+4*vertEffMRS(gct,(2*iz+step)/2,mA)+vertEffMRS(gct,iz+step,mA));
  }
  return integral;
}



double calcSignal(double E0, double eps2, double mA, double zCut, double Nbin, double mrange){
  double gg = gct(E0,eps2,mA);
  if (Nbin !=0){
    return radFrac(mA)*Nbin*(3*3.14*eps2/(2*alpha))*(mA/1000/mrange)*epsBin*integrateZ(gg,mA,zCut);
  }
  else {return 0;}
}

/*
 * This is the code which reads in the 2D histogram, slices it, and plots the reach.
 */

//Read in the 2d histo
int main(){
  
  TFile *f = new TFile("../output_L1L1/ub/cutPlots.root");
  TFile *f1 = new TFile("../finalZvMub.root");
  TFile *fout = new TFile("../output_L1L1/ub/outputFits.root","RECREATE");
  gROOT->SetBatch(true);

  //get z vs mass plot
  TH2F *h2 = new TH2F("h2","zVertex vs Mass",400,0.0,0.1,400,-100,100);
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
  Double_t par0[nBinsT];
  Double_t par1[nBinsT];
  Double_t par3[nBinsT];
  Double_t par4[nBinsT];
  Double_t par2[nBinsT];
  Double_t zCutMax[nBinsT];
  const double quantile = 0.9;
  

  TCanvas *canvas = new TCanvas("canvas"," vertex slice plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLogy();
  std::string pdf_file_name = "../output_L1L1/ub/vertexSliceFits.pdf";

  
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
    par0[ii] = vFit->GetParameter(0);
    par1[ii] = vFit->GetParameter(1);
    par2[ii] = vFit->GetParameter(2);
    par3[ii] = vFit->GetParameter(3);
    par4[ii] = vFit->GetParameter(4);

    //nEventsBin[ii] = v1[ii]->Integral();//
    nEventsBin[ii] = hc->Integral(lowBin,highBin,0,200);

    // nEventsBin[ii] = v1[ii]->Integral(v1[ii]->FindBin(vFit->GetParameter(1)),v1[ii]->FindLastBinAbove(0));//hc->Integral(lowBin,highBin,95,100);
    //zCut[ii] = fitZcut(centermass[ii]);
    //zCutMax[ii] = zCut[ii]-TMath::Log(1.-quantile);///*vFit->GetParameter(4);
    //zCutScaled[ii] = zCut[ii]-vFit->GetParameter(4)*TMath::Log(0.1);
  

    //Use the fit:
    zCut[ii] = fitZcutUB(centermass[ii]);//vFit->GetX(0.5*massrange[ii]/reslimited[ii]/vFit->GetParameter(4),zmean[ii],v1[ii]->GetBinCenter(v1[ii]->FindLastBinAbove(0)));
    //zCut[ii] = vFit->GetX(0.5*massrange[ii]/reslimited[ii]/vFit->GetParameter(4),zmean[ii],v1[ii]->GetBinCenter(v1[ii]->FindLastBinAbove(0)));
    zCutScaled[ii] = fitZcut(centermass[ii])-TMath::Log(0.1);
    zCutMax[ii] = zCut[ii]-TMath::Log(1.-quantile);
    
    gStyle->SetOptFit(111);
       
    //calculate the signal yield at each mass
    for (int ll=0; ll<neps; ll++){
      double ieps = h_reach->GetYaxis()->GetBinCenter(ll);
      double imass = h_reach->GetXaxis()->GetBinCenter(ii);      
      //   yield[ll][ii] = calcSignal(1.05, ll*pow(10,-10), centermass[ii]*1000, zCut[ii], nEventsBin[ii], massrange[ii]);
      yield[ll][ii] = calcSignal(1.05, ieps, imass*1000, zCut[ii], nEventsBin[ii], massrange[ii]);
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
  ccz->SaveAs("../output_L1L1/ub/zcut.C");
  ccz->Close();

 //plot the z resolution values
  TCanvas *czr = new TCanvas("czr","z resolution",1000,900);
  czr->cd();
  gStyle->SetOptStat(0);
  TGraph *zzres = new TGraph(nBinsT,centermass,par2);
  zzres->Draw("alp");
  czr->Update();
  czr->SaveAs("../output_L1L1/ub/zres.C");
  czr->Close();

  
  
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
  //gxMax->Draw("lsame");
  gxS->SetLineColor(kMagenta);
  gxS->SetLineWidth(2);
  gxS->Draw("lsame");
  kk1->Update();
  kk1->SaveAs("../output_L1L1/ub/zVm_zCutMRS.C");
  kk1->Close();

  TCanvas *kk2 = new TCanvas("kk2","par3 for mass",800,800);
  kk2->cd();
  TGraph *gx2 = new TGraph(nBinsT,centermass,par3);
  gx2->Draw("");
  kk2->Update();
  kk2->SaveAs("../output_L1L1/ub/par3_zCut.C");
  kk2->Close();

  TCanvas *kk3 = new TCanvas("kk3","par4 for mass",800,800);
  kk3->cd();
  TGraph *gx3 = new TGraph(nBinsT,centermass,par4);
  gx3->Draw("");
  kk3->Update();
  kk3->SaveAs("../output_L1L1/ub/par4_zCut.C");
  kk3->Close();


  TCanvas *kk4 = new TCanvas("kk4","par0 for mass",800,800);
  kk4->cd();
  TGraph *gx4 = new TGraph(nBinsT,centermass,par0);
  gx4->Draw("");
  kk4->Update();
  kk4->SaveAs("../output_L1L1/ub/par0_zCut.C");
  kk4->Close();

  TCanvas *kk5 = new TCanvas("kk5","par1 for mass",800,800);
  kk5->cd();
  TGraph *gx5 = new TGraph(nBinsT,centermass,par1);
  gx5->Draw("");
  kk5->Update();
  kk5->SaveAs("../output_L1L1/ub/par1_zCut.C");
  kk5->Close();

  TCanvas *kk6 = new TCanvas("kk6","par2 for mass",800,800);
  kk6->cd();
  TGraph *gx6 = new TGraph(nBinsT,centermass,par2);
  gx6->Draw("");
  kk6->Update();
  kk6->SaveAs("../output_L1L1/ub/par2_zCut.C");
  kk6->Close();
  
 
  /////////////////////////////////////
  //Plot reach////////////////////////			   
  TCanvas *cc = new TCanvas("cc","Reach",1000,900);
  cc->SetLogy();
  //cc->SetLogx();
  cc->cd();
  gStyle->SetOptStat(0);
  //h_reach->GetZaxis()->SetRangeUser(0.0,0.25);//25);
  h_reach->GetXaxis()->SetRangeUser(0.01,0.07);
  double maxSig = h_reach->GetMaximum();
  const int nContours = 10;
  double contours[nContours];
  for (int kk=0; kk<nContours; kk++){
    double mult = kk*0.1;
    contours[kk] = mult*maxSig;
  }
  h_reach->SetContour(nContours,contours);
  h_reach->Draw("contz");
  h_reach->GetXaxis()->SetRangeUser(0.001,0.1);
  h_reach->GetXaxis()->SetTitle("Mass [GeV]");
  h_reach->GetXaxis()->SetTitleOffset(1.3);
  h_reach->GetYaxis()->SetTitle("#epsilon^{2}");
  cc->Update();
  cc->SaveAs("../output_L1L1/ub/reach.C");
  cc->Close();  

  fout->Write();
  fout->Close(); 

}//end main



