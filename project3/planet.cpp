#include "planet.h"

Planet:: Planet () { mass = xPosition = yPosition = xVelocity = yVelocity = 0.0;}

Planet:: Planet (double mass_, double xPosition_, double yPosition_, double xVelocity_, double yVelocity_)
{
    mass = mass_;
    xVelocity = xVelocity_;
    yVelocity = yVelocity_;
    xPosition = xPosition_;
    yPosition = yPosition_;

}

double Planet::getMass() const
{
    return mass;
}

double Planet:: getInitialXPosition() const {return xPosition;}
double Planet:: getInitialYPosition() const {return yPosition;}
double Planet:: getInitialXVelocity() const {return xVelocity;}
double Planet:: getInitialYVelocity() const {return yVelocity;}
double Planet:: getKineticEnergy(double vx_, double vy_)
{
    double velocity2 = vx_*vx_ + vy_*vy_;
    return 0.5*mass*velocity2;
}
double Planet:: getPotentialEnergy(double r_)
{
    return 6.67*pow(10,-11)*mass/r_;
}
