#ifndef Particle_h
#define Particle_h

#include "Vector3D.h"
#include <string>
#include <iostream>

using std::string;
using std::ostream;

class Particle{
 public:
  Particle(const string& name,const int& id ,const double& mass, const int& charge,const Vector3D& momento);
  Particle(const string& name, const int& id, const double& mass, const int& charge, const double& px=0.,const double& py=0.,const double& pz=0.);
  
  virtual void setMomentum(const Vector3D& momento);
  
  void setName(const string& name);
  void setMass(const double& mass);
  void setId(const int& id);
  void setCharge(const int& charge);
  

  double mass()const {return mass_;}
  string name()const {return name_;}
  int id()const {return id_;}
  int charge()const {return charge_;}
  Vector3D Momentum()const {return momento_;}
  double Energy()const;

  friend ostream& operator<<(ostream& os, const Particle& p);
  
 private:
  double mass_;
  int charge_,id_;
  Vector3D momento_ ;
  string name_;

  
};
#endif
