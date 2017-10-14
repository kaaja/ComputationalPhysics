#ifndef PLANET_H
#define PLANET_H

#include<iostream>
#include<new>
#include<vector>
#include<cmath>

using std::vector;
using namespace std;

class Planet
{
private:
    double mass = 0.0;
    //vector<double> position(2);
    //vector<double> velocity(2);
    double xPosition = 0.0;
    double yPosition = 0.0;
    double xVelocity = 0.0;
    double yVelocity = 0.0;
    double xForce, yForce, radialDistance;

public:
    Planet();
    Planet(double mass_, double xPosition_, double yPosition_, double xVelocity_, double yVelocity_);

    ~Planet () {} // end initializer

    double getMass() const;
    //vector<double> getInitialPosition() const;
    //vector<double> getInitialVelocity() const;

    double getXPosition();
    double getYPosition();
    double getRadialDistance(Planet OtherPlanet );
    double getXVelocity();
    double getYVelocity();
    double getKineticEnergy(double mass_, double vx_, double vy_);
    double getPotentialEnergy(double r_, double mass_);
    double getAngularMomentum(double r_, double mass_, double vx_, double vy_);
    void getAcceleration( vector<Planet> planets_, double *accelerationX_, double *accelerationY_, int numberOfPlanets_);
    void getAlternativeForce(double mass_, double x_, double y_, double r_, double *forceX_, double *forceY_, double beta_);
    void setXposition(double x_);
    void setYposition(double y_);
    void setDistance(double r_);
    void setXVelociy(double vx_);
    void setYVelociy(double vy_);


};

#endif // PLANET_H
