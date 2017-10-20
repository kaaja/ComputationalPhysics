#ifndef PLANETMERCURY_H
#define PLANETMERCURY_H
#include "planet.h"

class PlanetMercury : public Planet
{
public:
    PlanetMercury();
    PlanetMercury(double mass_,
                  double xPosition_,
                  double yPosition_,
                  double xVelocity_,
                  double yVelocity_,
                  string filename_,
                  string planetName_);
    virtual vector<double> getAcceleration( vector<Planet*> planets_, int numberOfPlanets_) override;
    virtual void writeTofile(double timeUsed_, double centerOfMassX_, double centerOfMassY_) override;


};

#endif // PLANETMERCURY_H
