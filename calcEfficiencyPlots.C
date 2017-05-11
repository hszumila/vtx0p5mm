

const double zMax = 100;

//for the L1L1 dataset
double vertEffL1L1(double gct, double zz, double mA){
   return 1./gct*exp((-5-zz)/gct)*exp((-1.08E-5+8.1E-7*mA-1.392E-8*mA*mA)*pow(zz,3)
				     +(0.001578-0.000113*mA+1.567E-6*mA*mA)*pow(zz,2)
				      +(-0.066555+0.002836*mA-3.1037E-5*mA*mA)*zz
				      +(-0.316+0.0135*mA-0.000154*mA*mA));
    }

double fitZcutL1L1(double mass){
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
  return 22.4+896.463*mass-29940.4*pow(mass,2)+213639*pow(mass,3);
}
/////////////////////////////////////////////////////
// Crystal ball peak function
Double_t CrystalBall(double z, double mass) {
 
  
  double par0 = 1.0;
  double par1 = 141.5;
  double par2 = -35.1007+3.26202*mass-0.0272688*pow(mass,2);
  double par3 = 10.0466+0.0661*mass;
  double par4 = -0.708591+0.0545465*mass-0.000664232*pow(mass,2);
    
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





double vertEffL1L2(double gct, double zz, double mA){
 
  //if (mA<20){mA=20;}
  //else if (mA>50){mA=50;}
  return 1./gct*exp((-5-zz)/gct)*CrystalBall(zz,mA);
  
}













//////////////////////////////////////////////////////
//for the L1L2 dataset

double vertEffL1L2G(double gct, double zz, double mA){
   return 1./gct*exp((-5-zz)/gct)*(-0.5526+0.044312*mA-0.0005447*pow(mA,2))*exp(-pow( ( zz-(-4.5487+1.3203*mA))/(2*(6.9132+0.092927*mA)) ,2) );
    }


double fitZcutL1L2(double mass){
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
  return -23.16+5906*mass-153311*pow(mass,2)+1080900*pow(mass,3);
}


//for the L2L2 dataset
double vertEffL2L2(double gct, double zz, double mA){
   return 1./gct*exp((-5-zz)/gct)*(-0.03242+0.00514*mA)*exp(-pow( ( zz-(12.4881+1.6875*mA))/(2*(4.7982+0.11158*mA)) ,2) );
}


double fitZcutL2L2(double mass){
  return -17.76 +8563*mass-293498*mass*mass+2678060*mass*mass*mass;
}


//here's the integral part
//mA in MeV, gct in mm, zCut in mm
double integrateZ0(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffL1L1(gct,iz,mA)+4*vertEffL1L1(gct,(2*iz+step)/2,mA)+vertEffL1L1(gct,iz+step,mA));
  }
  return integral;
}

double integrateZ1(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  double integralG = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffL1L2(gct,iz,mA)+4*vertEffL1L2(gct,(2*iz+step)/2,mA)+vertEffL1L2(gct,iz+step,mA));
  }
  for (double iz=zCut; iz<=zMax; iz+=step){
    integralG += (step/6)*(vertEffL1L2G(gct,iz,mA)+4*vertEffL1L2G(gct,(2*iz+step)/2,mA)+vertEffL1L2G(gct,iz+step,mA));
  }
  double diff = integralG - integral;
  cout<<"gct:\t"<<gct<<"\tmA:\t"<<mA<<"\tzCut:\t"<<zCut<<"\tintegral\t"<<integral<<"\tgaussianint\t"<<integralG<<"\tdifference\t"<<diff<<endl;
  return integral;
}
double integrateZ2(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    //integral += (step/6)*(vertEffL2L2(gct,iz,mA)+4*vertEffL2L2(gct,(2*iz+step)/2,mA)+vertEffL2L2(gct,iz+step,mA));
    integral += (step/6)*(vertEffL1L2G(gct,iz,mA)+4*vertEffL1L2G(gct,(2*iz+step)/2,mA)+vertEffL1L2G(gct,iz+step,mA));
  }
  return integral;
}
//E0 [GeV], mA [MeV]
//returns [mm]
double gct(double E0, double eps2, double mA){
  return 8*(E0/10)*(pow(0.0001,2)/eps2)*pow(100/mA,2);
}


void calcEfficiencyPlots(){
  gROOT->SetBatch(true);

  //define by layer (L1L1=0,L1L2=2,L2L2=3)
  const int ntypes = 3;
  
  //loop over mass
  const int nmass = 90;
  double mass[nmass];
  for (int ii=0;ii<nmass;ii++){
    mass[ii] = ii;
  }
  
 

  //loop over eps2
  const int neps2 = 13;
  double eps2[neps2] = {7E-10, 8E-10, 9E-10, 1E-9, 2E-9, 3E-9, 4E-9, 5E-9, 6E-9, 7E-9, 8E-9, 9E-9, 1E-8};
  int eps2cheat[neps2] = {7,8,9,1,2,3,4,5,6,7,8,9,1};
  int eps2cheat2[neps2] = {10,10,10,9,9,9,9,9,9,9,9,9,8};
  double gctf[neps2][nmass];
  const int nzpos = 115;
  int zpos[nzpos];
  for (int jj=0; jj<nzpos; jj++){
    zpos[jj] = jj-5;//mm
  }
  
  //store the integral values
  double integral[ntypes][nzpos][neps2][nmass];  

  
  for(int kk=0; kk<nzpos; kk++){
    for (int ii=0; ii<nmass; ii++){
      for (int jj=0; jj<neps2;jj++){
	gctf[jj][ii] = gct(1.056,eps2[jj],mass[ii]);
	double zz = zpos[kk];
	integral[0][kk][jj][ii] = integrateZ0(gctf[jj][ii], mass[ii], zz);
	integral[1][kk][jj][ii] = integrateZ1(gctf[jj][ii], mass[ii], zz);
	integral[2][kk][jj][ii] = integrateZ2(gctf[jj][ii], mass[ii], zz);
      }
    }
  }


  
  //make canvas
  TCanvas *canvas = new TCanvas("canvas"," Integral Plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLeftMargin(0.1);
  canvas->SetRightMargin(0.2);
  std::string pdf_file_name = "../integralResults2.pdf";
  canvas->Update();
  gStyle->SetOptStat(0);

  ///////////////////////////////////////////////
  //Point for drawing the zCut
 
  const int nPts = 90;//55;//35;
  const int nPts2 = 90;//20;
  double cmassX[nPts];
  double zcutX[ntypes][nPts];
  for (int ijk=0; ijk<nPts; ijk++){
    cmassX[ijk] = ijk;
    zcutX[0][ijk] = fitZcutL1L1(cmassX[ijk]/1000);
    zcutX[1][ijk] = fitZcutL1L2(cmassX[ijk]/1000);
    cout<<"mass:\t"<<cmassX[ijk]<<"\t zCut:\t"<<zcutX[0][ijk]<<endl;
    
    //zcutX[2][ijk] = fitZcutL2L2(cmassX[ijk]/1000);
  }
  for (int ijk=0; ijk<nPts2; ijk++){
    cmassX[ijk] = ijk;//+15+0.5*ijk;
    zcutX[2][ijk] = fitZcutL2L2(cmassX[ijk]/1000);


  }
  


  //plots to make:  // double integral[ntypes][nzpos][neps2][nmass];  
  TH2F *h2_mVz[ntypes];
  TH2F *h2_epsVz[ntypes];
  const char *l1l1 = "L1L1, #epsilon^{2}=5E-9; vertex z position [mm]; mass [MeV];full integral where zCut = z vertex position";
  const char *l1l2 = "L1L2, #epsilon^{2}=5E-9; vertex z position [mm]; mass [Mev];full (CB) integral where zCut = z vertex position";
  const char *l2l2 = "L1L2, #epsilon^{2}=5E-9; vertex z position [mm]; mass [MeV];full (Gaussian) integral where zCut = z vertex position";
  const char *l1l1e = "L1L1, mass=35MeV; vertex z position [mm]; #epsilon^{2};full integral where zCut = z vertex position";
  const char *l1l2e = "L1L2, mass=35MeV; vertex z position [mm]; #epsilon^{2};full (CB) integral where zCut = z vertex position";
  const char *l2l2e = "L1L2, mass=35MeV; vertex z position [mm]; #epsilon^{2};full (Gaussian approx) integral where zCut = z vertex position";
 
  for (int kk=0; kk<ntypes; kk++){
    const char *layer;
    const char *layere;
    if (kk==0){layer = l1l1;layere = l1l1e;}
    else if (kk==1) {layer=l1l2;layere = l1l2e;}
    else {layer = l2l2;layere = l2l2e;}

    h2_mVz[kk] = new TH2F(Form("h2_mVz[%d]",kk),layer,117, -5, 100,100,0,90);
    h2_epsVz[kk] = new TH2F(Form("h2_epsVz[%d]",kk),layere,117, -5, 100,200,1E-10,4E-8);
    h2_mVz[kk]->GetZaxis()->SetRangeUser(1E-5,1E0); 
    h2_epsVz[kk]->GetZaxis()->SetRangeUser(1E-5,1E0); 
  }

 	
  //mass vs z, weighted by efficiency
  for (int ik=0; ik<ntypes; ik++){
    for (int iz=0; iz<nzpos; iz++){
      int zbin = h2_mVz[ik]->GetXaxis()->FindBin(zpos[iz]);

      for (int ii=0; ii<nmass;ii++){
	int mbin = h2_mVz[ik]->GetYaxis()->FindBin(mass[ii]);
	h2_mVz[ik]->SetBinContent(zbin,mbin,integral[ik][iz][7][ii]);
      }
      for (int jj=0; jj<neps2;jj++){
	int ebin = h2_epsVz[ik]->GetYaxis()->FindBin(eps2[jj]);
	h2_epsVz[ik]->SetBinContent(zbin,ebin,integral[ik][iz][jj][6]);
      }
      
    }
  }

  //gct (fix epsilon) vs z, weighted by efficiency

  //gct (fix mass) vs z, weighted by efficiency
  canvas->SetLogz();
  TGraph *gx[ntypes];

  for (int ik=0; ik<ntypes; ik++){

    h2_mVz[ik]->Draw("colz");
    h2_mVz[ik]->GetZaxis()->SetTitleOffset(1.5);
    if (ik!=2){
      gx[ik] = new TGraph(nPts,zcutX[ik],cmassX);
      gx[ik]->SetLineColor(kRed);
      gx[ik]->SetLineWidth(3);
      gx[ik]->Draw("lsame");
    }

    else{
      gx[ik] = new TGraph(nPts2,zcutX[ik],cmassX);
      gx[ik]->SetLineColor(kRed);
      gx[ik]->SetLineWidth(3);
      //gx[ik]->Draw("lsame");
    }
    canvas->Print( (pdf_file_name + "(").c_str());

  }
  canvas->SetLogy();

  for (int ik=0; ik<ntypes; ik++){

    h2_epsVz[ik]->Draw("colz");
    h2_epsVz[ik]->GetZaxis()->SetTitleOffset(1.5);
    canvas->Print( (pdf_file_name + "(").c_str());

  }
    
  canvas->Print( (pdf_file_name + ")").c_str());


}
