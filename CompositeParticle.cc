#include "CompositeParticle.h"
#include <iostream>
#include <string>
#include "TRandom3.h"
#include <math.h>


using std::cout;
using std::endl;
using std::ostream;
using std::string;
CompositeParticle::CompositeParticle(const string& name,const int& id, const double& mass, const int& charge):Particle(name,id,mass,charge,0.){}

//Lorentz Boost 3D
Vector3D CompositeParticle::Boost(const Vector3D& betacdm, const Vector3D& momento, const double& Energia){
   //ho la formula generale del boost dal lab al sistema in moto, per ottenere la matrice inversa cambio segno al boost
  Vector3D beta=betacdm*-1;
  double gamma=1/sqrt(1-beta.mod()*beta.mod()); //gamma cdm

  //equzioni generali del boost
  double px=-gamma*beta.x()*Energia + (1+(gamma-1)*(beta.x()*beta.x())/(beta.mod()*beta.mod()))*momento.x() + (gamma-1)*(beta.x()*beta.y())/(beta.mod()*beta.mod())*momento.y() + (gamma-1)*(beta.x()*beta.z())/(beta.mod()*beta.mod())*momento.z();
  double py=-gamma*beta.y()*Energia + (1+(gamma-1)*(beta.y()*beta.y())/(beta.mod()*beta.mod()))*momento.y() + (gamma-1)*(beta.x()*beta.y())/(beta.mod()*beta.mod())*momento.x() + (gamma-1)*(beta.y()*beta.z())/(beta.mod()*beta.mod())*momento.z();
  double pz=-gamma*beta.z()*Energia + (1+(gamma-1)*(beta.z()*beta.z())/(beta.mod()*beta.mod()))*momento.z() + (gamma-1)*(beta.z()*beta.y())/(beta.mod()*beta.mod())*momento.y() + (gamma-1)*(beta.x()*beta.z())/(beta.mod()*beta.mod())*momento.x();
  //ritorno le componenti boostate
  return Vector3D(px,py,pz);
}

void CompositeParticle::setMomentum(const Vector3D& momento){
  Particle::setMomentum(momento);
  if(particles_.size()==2){
    double EnergyFather=sqrt(mass()*mass()+Momentum().ScalarProduct(Momentum()));
    Vector3D betacdm=Momentum()/EnergyFather;
    TRandom3*  gen = new TRandom3();
    gen->SetSeed(0);
    double phi= gen->Uniform(0,2*M_PI);
    double theta=gen->Uniform(0,M_PI);
    double moduloP=sqrt(pow(mass(),4)+pow(particles_[0]->mass()*particles_[0]->mass()-particles_[1]->mass()*particles_[1]->mass(),2) -2*pow(mass(),2)*(particles_[0]->mass()*particles_[0]->mass() + particles_[1]->mass()*particles_[1]->mass()))/(2*mass());
    Vector3D p1(sin(theta)*cos(phi)*moduloP , sin(theta)*sin(phi)*moduloP , cos(theta)*moduloP);
    Vector3D p2=-1*p1;
    particles_[0]->setMomentum(Boost(betacdm,p1,sqrt(pow(particles_[0]->mass(),2)+pow(moduloP,2))));
    particles_[1]->setMomentum(Boost(betacdm,p2,sqrt(pow(particles_[1]->mass(),2)+pow(moduloP,2))));

    //calcolo massa invariante
    
    /* double Minv=sqrt(pow((particles_[0]->Energy()+particles_[1]->Energy()),2)-(pow(particles_[0]->Momentum().mod(),2)+pow(particles_[1]->Momentum().mod(),2)+2*particles_[0]->Momentum().ScalarProduct(particles_[1]->Momentum())));
       cout<<"Massa invariante di "<<particles_[1]->name()<<" e "<<particles_[0]->name()<<" di valore: "<<Minv<<endl;*/
        
    delete gen;
  }else{
    cout <<"Impossibile eseguire il decadimento a 2 corpi, numero di particelle figlie introdotte diverso da 2."<< endl;
  }
}


void CompositeParticle::add(Particle* part){
  if(particles_.size()>2){
    cout<<"Numero di Particelle figlie inserite maggiore di 2"<<endl;
    particles_.push_back(part);
  }else{
    particles_.push_back(part);
  }
}

void CompositeParticle::print()const{
  if(particles_.size()==2){
    cout<<"Particella Padre: "<<name()
	<<" con Id:  "<<id()
	<<", massa:  "<<mass()
	<<" GeV impulso: "<<Momentum()
	<<" GeV e carica " <<charge()
	<<endl;
    cout<<"Particelle figlie prodotte:"<<endl;
    for(std::vector<Particle*>::const_iterator it =particles_.begin(); it !=particles_.end(); ++it) {
      cout<<"Particella 1: "<< (*it)->name()
	  <<", con massa: "<< (*it)->mass()
	  <<" GeV, Impulso: "<<(*it)->Momentum()
	  <<" GeV ed Energia: "<<(*it)->Energy()
	  <<endl;
    }
  }else{
    cout<<"Impossibile eseguire il decadimento a 2 corpi"<<endl;
  }
}
