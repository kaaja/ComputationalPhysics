#ifndef PLANET_H
#define PLANET_H

#include<iostream>
#include<new>
#include<vector>
#include<cmath>
#include <fstream>
#include <iomanip>
#include <string>

using std::vector;
using namespace std;

class Planet
{
protected:
    double mass = 0.0;
    //vector<double> position(2);
    //vector<double> velocity(2);
    double xPosition = 0.0;
    double yPosition = 0.0;
    double xVelocity = 0.0;
    double yVelocity = 0.0;
    double time, step;
    double xForce, yForce, radialDistance;
    string filename, planetName;

public:
    Planet();
    Planet(double mass_, double xPosition_, double yPosition_, double xVelocity_, double yVelocity_, string filename_, string planetName_);

    ~Planet () {} // end initializer

    string getPlanetName() const;
    double getMass() const;
    //vector<double> getInitialPosition() const;
    //vector<double> getInitialVelocity() const;

    double getXPosition();
    double getYPosition();
    double getRPosition();
    double getRadialDistance(Planet &OtherPlanet );
    double getXVelocity();
    double getYVelocity();
    double getKineticEnergy();
    double getPotentialEnergy();
    double getAngularMomentum();
    virtual void getAcceleration( vector<Planet*> planets_, double *accelerationX_, double *accelerationY_, int numberOfPlanets_);
    void getAlternativeForce(double mass_, double x_, double y_, double r_, double *forceX_, double *forceY_, double beta_);
    void setStep(double step_);
    void setXposition(double x_);
    void setYposition(double y_);
    void setDistance(double r_);
    void setXVelociy(double vx_);
    void setYVelociy(double vy_);
    void setTime(double time_);
    virtual void writeTofile(double timeUsed_, double centerOfMassX_, double centerOfMassY_);
};

#endif // PLANET_H
