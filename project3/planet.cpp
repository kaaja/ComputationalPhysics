#include "planet.h"

Planet:: Planet () { mass = xPosition = yPosition = xVelocity = yVelocity = 0.0;}

Planet:: Planet (double mass_, double xPosition_, double yPosition_, double xVelocity_, double yVelocity_)
{
    mass = mass_;
    xVelocity = xVelocity_;
    yVelocity = yVelocity_;
    xPosition = xPosition_;
    yPosition = yPosition_;
    radialDistance = sqrt(xPosition*xPosition + yPosition*yPosition);

}

double Planet::getMass() const
{
    return mass;
}

double Planet:: getXPosition() {return xPosition;}
double Planet:: getYPosition() {return yPosition;}
double Planet:: getRadialDistance(Planet otherPlanet_)
{
    double xPositionOtherPlanet = otherPlanet_.getXPosition();
    double yPositionOtherPlanet = otherPlanet_.getYPosition();
    double xDistance = xPosition - xPositionOtherPlanet;
    double yDistance = yPosition - yPositionOtherPlanet;
    double r = sqrt(xDistance*xDistance + yDistance*yDistance);
    return r;
}
double Planet:: getXVelocity() {return xVelocity;}
double Planet:: getYVelocity() {return yVelocity;}
double Planet:: getKineticEnergy(double mass_, double vx_, double vy_)
{
    double velocity2 = vx_*vx_ + vy_*vy_;
    return 0.5*mass_*velocity2;
}
double Planet:: getPotentialEnergy(double r_, double mass_)
{
    return 6.67*pow(10,-11)*mass/r_;
}
double Planet:: getAngularMomentum(double r_, double mass_, double vx_, double vy_)
{
    double velocity = sqrt(vx_*vx_ + vy_*vy_);
    return r_*mass_*velocity;
}
void Planet:: getAcceleration(vector<Planet> planets_, double *accelerationX_, double *accelerationY_, int numberOfPlanets_)
{
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    double r;
    double forceX, forceY;
    forceX = -FourPi2*xPosition*mass/getRadialDistance(planets_[0]);
    forceY = -FourPi2*yPosition*mass/getRadialDistance(planets_[0]);
    for (int planetNumber = 1; planetNumber < numberOfPlanets_; planetNumber++)
    {
        if (planetNumber != 1){
            r = getRadialDistance(planets_[planetNumber]);
            forceX += FourPi2*planets_[planetNumber].getMass()/planets_[0].getMass()*(xPosition - planets_[planetNumber].getXPosition())/pow(r,3);
            forceY += FourPi2*planets_[planetNumber].getMass()/planets_[0].getMass()*(yPosition - planets_[planetNumber].getYPosition())/pow(r,3);
        }
    }
    *accelerationX_ = forceX/mass;
    *accelerationY_ = forceY/mass;
}

void Planet:: getAlternativeForce(double mass_, double x_, double y_, double r_, double *forceX_, double *forceY_, double beta_)
{
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    *forceX_ = -FourPi2*x_ *mass/(pow(r_,beta_));
    *forceY_ = -FourPi2*y_*mass/(pow(r_,beta_));
}

void Planet:: setXposition(double x_){ xPosition = x_;}
void Planet:: setYposition(double y_){yPosition = y_;}
void Planet:: setDistance(double r_){radialDistance = r_;}
void Planet:: setXVelociy(double vx_){xVelocity = vx_;}
void Planet:: setYVelociy(double vy_){yVelocity = vy_;}
