

const double zMax = 100;

//for the L1L1 dataset
double vertEffL1L1(double gct, double zz, double mA){
   return 1./gct*exp((-5-zz)/gct)*exp((-1.08E-5+8.1E-7*mA-1.392E-8*mA*mA)*pow(zz,3)
				     +(0.001578-0.000113*mA+1.567E-6*mA*mA)*pow(zz,2)
				      +(-0.066555+0.002836*mA-3.1037E-5*mA*mA)*zz
				      +(-0.316+0.0135*mA-0.000154*mA*mA));
    }

double fitZcutL1L1(double mass){
  //return (15.86 + 1857*mass - 62743.2*mass*mass + 536084*mass*mass*mass-5.0);
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
  return 22.4+896.463*mass-29940.4*pow(mass,2)+213639*pow(mass,3);
}

//for the L1L2 dataset
double vertEffL1L2(double gct, double zz, double mA){
   return 1./gct*exp((-5-zz)/gct)*(-0.5526+0.044312*mA-0.0005447*pow(mA,2))*exp(-pow( ( zz-(-4.5487+1.3203*mA))/(2*(6.9132+0.092927*mA)) ,2) );
    }

double fitZcutL1L2(double mass){
  return -7.88 +5171*mass-149040*mass*mass+1162420*mass*mass*mass;
}

//for the L2L2 dataset
double vertEffL2L2(double gct, double zz, double mA){
   return 1./gct*exp((-5-zz)/gct)*(0.03242+0.00514*mA)*exp(-pow( ( zz-(12.4881+1.6875*mA))/(2*(4.7982+0.11158*mA)) ,2) );
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
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffL1L2(gct,iz,mA)+4*vertEffL1L2(gct,(2*iz+step)/2,mA)+vertEffL1L2(gct,iz+step,mA));
  }
  return integral;
}
double integrateZ2(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffL2L2(gct,iz,mA)+4*vertEffL2L2(gct,(2*iz+step)/2,mA)+vertEffL2L2(gct,iz+step,mA));
  }
  return integral;
}
//E0 [GeV], mA [MeV]
//returns [mm]
double gct(double E0, double eps2, double mA){
  return 8*(E0/10)*(pow(0.0001,2)/eps2)*pow(100/mA,2);
}


void calcIntegral(){
  gROOT->SetBatch(true);

  //define by layer (L1L1=0,L1L2=2,L2L2=3)
  const int ntypes = 3;
  
  //loop over mass
  const int nmass = 11;
  double mass[nmass] = {20,22,24,26,28,30,35,40,50,60,70};

  //loop over eps2
  const int neps2 = 5;
  double eps2[neps2] = {1E-9,3E-9,5E-9,7E-9,9E-9};
  int eps2cheat[neps2] = {1,3,5,7,9};
  double gctf[neps2][nmass];

  //store the integral values
  double integral[ntypes][neps2][nmass];

  for (int ii=0; ii<nmass; ii++){
    for (int jj=0; jj<neps2;jj++){
      gctf[jj][ii] = gct(1.056,eps2[jj],mass[ii]);
      integral[0][jj][ii] = integrateZ0(gctf[jj][ii], mass[ii], fitZcutL1L1(mass[ii]/1000));
      integral[1][jj][ii] = integrateZ1(gctf[jj][ii], mass[ii], fitZcutL1L2(mass[ii]/1000));
      integral[2][jj][ii] = integrateZ2(gctf[jj][ii], mass[ii], fitZcutL2L2(mass[ii]/1000));
      //cout<<"gct\t"<<gctf<<"\tintegral\t"<<integral[0][jj][ii]<<"\t"<<integral[1][jj][ii]<<"\t"<<integral[2][jj][ii]<<endl;
    }
  }

  //make canvas
  TCanvas *canvas = new TCanvas("canvas"," Integral Plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLogy();
  std::string pdf_file_name = "../integralResults.pdf";
  canvas->Update();

  //plot the integrated efficiency values for each dataset as function of mass for fixed gct
  TMultiGraph *mg[neps2];
  //for (int kk=0;kk<neps2; kk++){
  int kk=0;
    mg[kk] = new TMultiGraph();
    mg[kk]->SetTitle(Form("#epsilon^{2} = %dE-9;A' Mass [MeV];Integral [zCut, 100mm]",eps2cheat[kk]));
    TGraphErrors *g0 = new TGraphErrors(nmass,mass,integral[0][0]);//kk]);
    TGraphErrors *g1 = new TGraphErrors(nmass,mass,integral[0][2]);//kk]);
    TGraphErrors *g2 = new TGraphErrors(nmass,mass,integral[0][4]);//kk]);
    g0->SetLineColor(kRed);
    g0->SetLineWidth(2);
    g1->SetLineColor(kBlue);
    g1->SetLineWidth(2);
    g2->SetLineColor(kCyan+2);
    g2->SetLineWidth(2);
    mg[kk]->Add(g0);
    mg[kk]->Add(g1);
    mg[kk]->Add(g2);
    mg[kk]->Draw("a");
    mg[kk]->GetYaxis()->SetRangeUser(0.0001,1.0);

    canvas->Print( (pdf_file_name + "(").c_str());
    //}
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle(";A'mass [MeV];decay length [mm]");
  TGraphErrors *gg[neps2];
  TLegend *lv = new TLegend(0.65,0.55,0.9,0.85);
  TLegendEntry *le[neps2];
  for (int ss = 0; ss<neps2; ss++){
    gg[ss] = new TGraphErrors(nmass,mass,gctf[ss]);
    gg[ss]->SetLineColor(kCyan+ss);
    gg[ss]->SetLineWidth(2);
    mg1->Add(gg[ss]);
    le[ss] = lv->AddEntry(Form("gg[%d]",ss),Form("#epsilon^{2} = %dE-9",eps2cheat[ss]),"0");
    le[ss]->SetTextColor(kCyan+ss);
  }
  mg1->Draw("a");
  lv->Draw(); 
  canvas->Print( (pdf_file_name + ")").c_str());


}
