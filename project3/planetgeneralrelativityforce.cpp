#include "planetgeneralrelativityforce.h"


PlanetGeneralRelativityForce::PlanetGeneralRelativityForce() : Planet(){}

PlanetGeneralRelativityForce::PlanetGeneralRelativityForce(double mass_,
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

vector<double> PlanetGeneralRelativityForce::getAcceleration(vector<Planet*> planets_, int numberOfPlanets_)
{
    vector<double> accelerations;
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    double rPlanetDistance;
    double c = 173.0*365.;
    double l, relativisticCorrection;
    double accelerationX_ = 0.;
    double accelerationY_ = 0.;
    int start  = 0;

    for (int planetNumber = start; planetNumber < numberOfPlanets_; planetNumber++)
    {
        if (planets_[planetNumber]->getMass() != mass){
            rPlanetDistance = getRadialDistance(*planets_[planetNumber]); // check the star!
            l = sqrt(pow((xPosition*yVelocity - yPosition*xVelocity), 2));
            relativisticCorrection = (3*l*l/(pow(rPlanetDistance,2)*c*c) );
            accelerationX_ += -FourPi2*planets_[planetNumber]->getMass()/planets_[0]->getMass()*(xPosition - planets_[planetNumber]->getXPosition())/pow(rPlanetDistance,3)*(1.+ relativisticCorrection);
            accelerationY_ += -FourPi2*planets_[planetNumber]->getMass()/planets_[0]->getMass()*(yPosition - planets_[planetNumber]->getYPosition())/pow(rPlanetDistance,3)*(1.+ relativisticCorrection);
        }
    }
    accelerations.push_back(accelerationX_);
    accelerations.push_back(accelerationY_);
    return accelerations;
}

void PlanetGeneralRelativityForce:: writeTofile(double timeUsed_, double centerOfMassX_, double centerOfMassY_)
{
    if(planetName=="Mercury")
    {
        if( (time>99.75 && step < 0.001) || (step >= 0.001))
        {
            Planet::writeTofile(timeUsed_, centerOfMassX_, centerOfMassY_);
        }
    }

}
