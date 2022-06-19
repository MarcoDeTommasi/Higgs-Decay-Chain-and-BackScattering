#ifndef CompositeParticle_h
#define CompositeParticle_h

#include "Particle.h"
#include <vector>

using std::vector;

class CompositeParticle: public Particle{
 public:
  CompositeParticle(const string& name, const int& id, const double& mass, const int& charge);

  virtual void setMomentum(const Vector3D& momento);
  Vector3D Boost(const Vector3D& beta, const Vector3D& momento, const double& Energia);
  
  void add(Particle* part);
  void print() const;
  
 private:
  vector<Particle*> particles_;

};
#endif
