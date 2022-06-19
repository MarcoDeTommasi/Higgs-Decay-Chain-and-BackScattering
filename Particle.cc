#include "Particle.h"
#include <math.h>

using std::cout;
using std::endl;
using std::string;
using std::ostream;

Particle::Particle(const string& name, const int& id, const double& mass,const int& charge ,const Vector3D& momento){
  name_=name;
  id_=id;
  mass_=mass;
  charge_=charge;
  momento_=momento;
}
Particle::Particle(const string& name, const int& id, const double& mass, const int& charge, const double& px,const double& py, const double& pz){
  name_=name;
  id_=id;
  mass_=mass;
  charge_=charge;
  momento_.setCoord(px,py,pz);
}
void Particle::setMomentum(const Vector3D& momento){
  momento_=momento;
}
void Particle::setName(const string& name){
  name_=name;
}
void Particle::setMass(const double& mass){
  mass_=mass;
}
void Particle::setId(const int& id){
  id_=id;
}
void Particle::setCharge(const int& charge){
  charge_=charge;
}
double Particle::Energy()const{
  return sqrt(mass_*mass_+momento_.ScalarProduct(momento_));
}

ostream& operator<<(ostream& os, const Particle& p){
  os<<"Particle: "<<p.name()
    <<" con ID: "<<p.id()
    <<" , carica: "<<p.charge()
    <<" ,  massa: "<<p.mass()
    <<" GeV,  momento: "<<p.Momentum()
    <<" GeV ed energia: "<<p.Energy()
    <<" GeV";
  return os;
}
