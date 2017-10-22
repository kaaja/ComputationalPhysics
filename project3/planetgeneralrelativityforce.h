#ifndef PLANETGENERALRELATIVITYFORCE_H
#define PLANETGENERALRELATIVITYFORCE_H
#include "planet.h"

class PlanetGeneralRelativityForce : public Planet
{
public:
    PlanetGeneralRelativityForce();
    PlanetGeneralRelativityForce(double mass_,
                  double xPosition_,
                  double yPosition_,
                  double xVelocity_,
                  double yVelocity_,
                  string filename_,
                  string planetName_);
    virtual vector<double> getAcceleration( vector<Planet*> planets_, int numberOfPlanets_) override;
    virtual void writeTofile(double timeUsed_, double centerOfMassX_, double centerOfMassY_) override;
};
#endif // PLANETGENERALRELATIVITYFORCE_H
