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
void Planet:: getForce(Planet otherPlanet, double *forceX_, double *forceY_)
{
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    *forceX_ = -FourPi2*xPosition *mass/radialDistance;
    *forceY_ = -FourPi2*yPosition*mass/radialDistance;

}
void Planet:: getAcceleration(double mass_, double *accelerationX, double *accelerationY, double forceX_, double forceY_)
{
    *accelerationX = forceX_/mass_;
    *accelerationY = forceY_/mass_;
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
