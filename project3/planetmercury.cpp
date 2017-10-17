#include "planetmercury.h"


PlanetMercury::PlanetMercury() : Planet(){}

PlanetMercury::PlanetMercury(double mass_,
                             double xPosition_,
                             double yPosition_,
                             double xVelocity_,
                             double yVelocity_,
                             string filename_,
                             string planetName_)
    : Planet(mass_,
             xPosition_,
             yPosition_,
             xVelocity_,
             yVelocity_,
             filename_,
             planetName_){}

void PlanetMercury::getAcceleration(vector<Planet> planets_, double *accelerationX_, double *accelerationY_, int numberOfPlanets_)
{
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    double rPlanetDistance;
    double c = 173.0/365.25;
    double l, relativisticCorrection;
    *accelerationX_ = 0.;
    *accelerationY_ = 0.;
    int start  = 0;

    for (int planetNumber = start; planetNumber < numberOfPlanets_; planetNumber++)
    {
        if (planets_[planetNumber].getMass() != mass){
            rPlanetDistance = getRadialDistance(planets_[planetNumber]);
            l = sqrt((xPosition*yVelocity)*(xPosition*yVelocity) + (- yPosition*xVelocity)*(- yPosition*xVelocity));
            relativisticCorrection = (3*l*l/(pow(rPlanetDistance,3)*c*c) );
            *accelerationX_ += -FourPi2*planets_[planetNumber].getMass()/planets_[0].getMass()*(xPosition - planets_[planetNumber].getXPosition())/pow(rPlanetDistance,3)*relativisticCorrection;
            *accelerationY_ += -FourPi2*planets_[planetNumber].getMass()/planets_[0].getMass()*(yPosition - planets_[planetNumber].getYPosition())/pow(rPlanetDistance,3)*relativisticCorrection;
        }
    }
}
