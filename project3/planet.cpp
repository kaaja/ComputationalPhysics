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
void Planet:: getForce(double mass_, double x_, double y_, double r_, double *forceX_, double *forceY_)
{
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    *forceX_ = -FourPi2*x_ *mass/(r_*r_*r_);
    *forceY_ = -FourPi2*y_*mass/(r_*r_*r_);

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
