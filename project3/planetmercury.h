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
    void getAcceleration( vector<Planet> planets_, double *accelerationX_, double *accelerationY_, int numberOfPlanets_);


};

#endif // PLANETMERCURY_H
