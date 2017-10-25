#include "planetalternativeforce.h"

PlanetAlternativeForce::PlanetAlternativeForce() : Planet(){}

PlanetAlternativeForce::PlanetAlternativeForce(double mass_,
                                     double xPosition_,
                                     double yPosition_,
                                     double xVelocity_,
                                     double yVelocity_,
                                     string filename_,
                                     string planetName_,
                                     double beta_)
    : Planet(mass_,
             xPosition_,
             yPosition_,
             xVelocity_,
             yVelocity_,
             filename_,
             planetName_)
{
    beta = beta_;
}
vector<double> PlanetAlternativeForce::getAcceleration(vector<Planet*> planets_, int numberOfPlanets_)
{
    vector<double> accelerations;
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    double rPlanetDistance;
    double accelerationX_ = 0.;
    double accelerationY_ = 0.;
    int start  = 0;

    for (int planetNumber = start; planetNumber < numberOfPlanets_; planetNumber++)
    {
        if (planets_[planetNumber]->getMass() != mass){
            rPlanetDistance = getRadialDistance(*planets_[planetNumber]);
            accelerationX_ += -FourPi2*planets_[planetNumber]->getMass()/planets_[0]->getMass()*(xPosition - planets_[planetNumber]->getXPosition())/pow(rPlanetDistance,beta);
            accelerationY_ += -FourPi2*planets_[planetNumber]->getMass()/planets_[0]->getMass()*(yPosition - planets_[planetNumber]->getYPosition())/pow(rPlanetDistance,beta);
        }
    }
    accelerations.push_back(accelerationX_);
    accelerations.push_back(accelerationY_);
    return accelerations;
}
