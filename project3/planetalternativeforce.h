#ifndef PLANETALTERNATIVEFORCE_H
#define PLANETALTERNATIVEFORCE_H
#include "planet.h"

class PlanetAlternativeForce : public Planet
{
private:
    double beta;
public:
    PlanetAlternativeForce();
    PlanetAlternativeForce(double mass_,
                      double xPosition_,
                      double yPosition_,
                      double xVelocity_,
                      double yVelocity_,
                      string filename_,
                      string planetName_, double beta_);
    virtual vector<double> getAcceleration( vector<Planet*> planets_, int numberOfPlanets_) override;

};

#endif // PLANETALTERNATIVEFORCE_H
