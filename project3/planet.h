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
public:
    Planet();
    Planet(double mass_, double xPosition_, double yPosition_, double xVelocity_, double yVelocity_);

    ~Planet () {} // end initializer

    double getMass() const;
    //vector<double> getInitialPosition() const;
    //vector<double> getInitialVelocity() const;

    double getInitialXPosition() const;
    double getInitialYPosition() const;
    double getInitialXVelocity() const;
    double getInitialYVelocity() const;
};

#endif // PLANET_H
