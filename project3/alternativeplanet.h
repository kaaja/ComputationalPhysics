#ifndef ALTERNATIVEPLANET_H
#define ALTERNATIVEPLANET_H
#include "planet.h"

class AlternativePlanet : public Planet
{
private:
    double beta;
public:
    AlternativePlanet();
    AlternativePlanet(double mass_,
                      double xPosition_,
                      double yPosition_,
                      double xVelocity_,
                      double yVelocity_,
                      string filename_,
                      string planetName_, double beta_);
    virtual vector<double> getAcceleration( vector<Planet*> planets_, int numberOfPlanets_) override;

};

#endif // ALTERNATIVEPLANET_H
