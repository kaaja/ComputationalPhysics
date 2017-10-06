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
