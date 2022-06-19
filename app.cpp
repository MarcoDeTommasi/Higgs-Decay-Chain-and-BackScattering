
#include "CompositeParticle.h"
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom3.h"

using std::cout;
using std::endl;
int main(){
  //Particle("name",Id,Mass,Charge,Momentum)
  //CompositeParticle("Name",Id,Mass,Charge)
  CompositeParticle h1("Higgs1",25,125,0);
  Particle phot("Fotone",22,0,0);
  Particle phot1("Fotone",22,0,0);
  h1.add(&phot);
  h1.add(&phot1);
  CompositeParticle h2("Higgs2",25,125,0),z1("Z1",32,91,0),z2("Z2",32,31,0);
  Particle ele("elettrone",11,0.00051,-1);
  Particle pos("positrone",-11,0.00051,+1);
  Particle mu("muone",13,0.139,-1);
  Particle mum("antimuone",-13,0.139,+1);
  h2.add(&z1);
  h2.add(&z2);
  z1.add(&ele);
  z1.add(&pos);
  z2.add(&mu);
  z2.add(&mum);
  
  //generazione del grafico richiesto per il decadimento dei fotoni (h1)

  double Pz0=50;
  double resol=0.1;
  int nDecay=10000;

  int nbins = 100;
  double xlo = -0.5*resol;
  double xhi = M_PI*(1+0.5*resol);
  double binwidth = (xhi-xlo) / nbins;
  cout << "# bins: " << nbins << "\t bin width: " << binwidth << std::endl;
  TH1F hx1("hx1", "Distribuzione degli angoli alpha di emissione dei fotoni", nbins, xlo, xhi );

  
  TRandom3* gen = new TRandom3();
  gen->SetSeed(0);

  for(int i=0;i<nDecay;i++){
    double temp= gen->Gaus(Pz0,Pz0*resol);
    Vector3D pz(0,0,temp);
    h1.setMomentum(pz);
    double theta=acos(phot.Momentum().ScalarProduct(phot1.Momentum())/(phot.Momentum().mod()*phot1.Momentum().mod()));
    hx1.Fill(theta);
  }

  xlo=-0.5*resol;
  xhi=125+0.8*resol;
  
  TH1F elet("elet", "Distrbuzione del Modulo degli impulsi P per elettroni", nbins, xlo, xhi );
  TH1F posit("posit", "Distrbuzione del Modulo degli impulsi P per positorni", nbins, xlo, xhi );
  TH1F muon("muon", "Distrbuzione del Modulo degli impulsi P per muoni", nbins, xlo, xhi );
  TH1F muonm("muonm", "Distrbuzione del Modulo degli impulsi P per antimuoni", nbins, xlo, xhi );


  //Generazione Istogrammi per la catena di decadimenti H-> Z1,Z2  Z1->e+,e-
  // Z2->mu+,mu-
  
  resol=0.2;
  Pz0=100;
  
  for(int i=0;i<nDecay;i++){
    double temp= gen->Gaus(Pz0,Pz0*resol);
    Vector3D pz(0,0,temp);
    h2.setMomentum(pz);
    elet.Fill(ele.Momentum().mod());
    posit.Fill(pos.Momentum().mod());
    muon.Fill(mu.Momentum().mod());
    muonm.Fill(mum.Momentum().mod());
  }
  
  //Realizzo gli istogrammi e li salvo come pdf. 
  TCanvas canv("canv", "canvas for plotting", 1280, 1024);

  //fotoni
  hx1.GetXaxis()->SetTitle("Distribution of alpha[rad]");
  hx1.GetYaxis()->SetTitle("Conteggi");
  hx1.Draw();
  hx1.SetLineColor(kBlue);
  canv.SaveAs("angoli_emissione.pdf");

  //prodotti di decadimento Z1 e Z2
  canv.Clear();
  canv.Divide(2,2);
  canv.cd(1);
  elet.Draw();
  elet.GetXaxis()->SetTitle("Modulo di P");
  elet.GetYaxis()->SetTitle("Conteggi");
  canv.cd(2);
  posit.Draw();
  posit.GetXaxis()->SetTitle("Modulo di P");
  posit.GetYaxis()->SetTitle("Conteggi");
  canv.cd(3);
  muon.Draw();
  muon.GetXaxis()->SetTitle("Modulo di P");
  muon.GetYaxis()->SetTitle("Conteggi");
  canv.cd(4);
  muonm.Draw();
  muonm.GetXaxis()->SetTitle("Modulo di P");
  muonm.GetYaxis()->SetTitle("Conteggi");
  canv.SaveAs("magnitude.pdf");


  delete gen;
 
}
