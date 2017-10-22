#include "solver.h"

Solver:: Solver(string solverType_) {solverType=solverType_;}

void Solver:: solve(double step_, vector<Planet *> planets_, int numberOfplanets_, int iterationStart_)
{
    step = step_;
    if (solverType == "ForwardEuler") forwardEuler(planets_, numberOfplanets_, iterationStart_);
    else if (solverType == "VelocityVerlet") velocityVerlet(planets_, numberOfplanets_, iterationStart_);
}


void Solver:: forwardEuler(vector<Planet *> planets_, int numberOfplanets_, int iterationStart)
{
    // Position
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        x = planets_[planetNumber]->getXPosition();
        y = planets_[planetNumber]->getYPosition();
        vx = planets_[planetNumber]->getXVelocity();
        vy = planets_[planetNumber]->getYVelocity();

        x +=  step*vx;
        y +=  step*vy;
        planets_[planetNumber]->setXposition(x);
        planets_[planetNumber]->setYposition(y);
    }
    //  Accelerations vector all planets_
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        acc= planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_);
        if(iterationStart==1 && planetNumber==iterationStart){
            accXVec.push_back(0.0);
            accYVec.push_back(0.0);
        }
        accXVec.push_back(acc[0]);
        accYVec.push_back(acc[1]);
        acc.clear();
    }
    // Velocity
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        vx = planets_[planetNumber]->getXVelocity();
        vy = planets_[planetNumber]->getYVelocity();

        vx +=  step*accXVec[planetNumber];
        vy +=  step*accYVec[planetNumber];
        planets_[planetNumber]->setXVelociy(vx);
        planets_[planetNumber]->setYVelociy(vy);

        // Write results
        planets_[planetNumber]->writeTofile(NAN, 0.0, 0.0);
    }
    // Reset acceleration vectors
    accXVec.clear();
    accYVec.clear();

}

void Solver:: velocityVerlet(vector<Planet *> planets_, int numberOfplanets_, int iterationStart)
{
    // Vector accelerations all planets
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        acc= planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_);
        if(iterationStart==1 && planetNumber==iterationStart){
            accXVec.push_back(0.0);
            accYVec.push_back(0.0);
        }
        accXVec.push_back(acc[0]);
        accYVec.push_back(acc[1]);
        acc.clear();
    }
    // Position
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        x = planets_[planetNumber]->getXPosition();
        y = planets_[planetNumber]->getYPosition();
        vx = planets_[planetNumber]->getXVelocity();
        vy = planets_[planetNumber]->getYVelocity();

        x +=  step*vx + step*step/2.0* accXVec[planetNumber];
        y +=  step*vy + step*step/2.0* accYVec[planetNumber];
        planets_[planetNumber]->setXposition(x);
        planets_[planetNumber]->setYposition(y);
    }
    // Vector old accelerations
    accXVecOld.clear();   accYVecOld.clear();
    accXVecOld = accXVec; accYVecOld = accYVec;
    accXVec.clear();      accYVec.clear();
    // Acceleration update
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        acc= planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_);
        if(iterationStart==1 && planetNumber==iterationStart){
            accXVec.push_back(0.0);
            accYVec.push_back(0.0);
        }
        accXVec.push_back(acc[0]);
        accYVec.push_back(acc[1]);
        acc.clear();
    }
    // Velocity
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        vx = planets_[planetNumber]->getXVelocity();
        vy = planets_[planetNumber]->getYVelocity();

        vx +=  step/2.0*(accXVecOld[planetNumber] + accXVec[planetNumber]);
        vy +=  step/2.0*(accYVecOld[planetNumber] + accYVec[planetNumber]);
        planets_[planetNumber]->setXVelociy(vx);
        planets_[planetNumber]->setYVelociy(vy);
        // write results
        planets_[planetNumber]->writeTofile(NAN, 0.0, 0.0);
    }
}

/*
void Solver:: velocityVerlet(vector<Planet *> planets_, int numberOfplanets_, int iterationStart)
{
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        x = planets_[planetNumber]->getXPosition();
        y = planets_[planetNumber]->getYPosition();
        vx = planets_[planetNumber]->getXVelocity();
        vy = planets_[planetNumber]->getYVelocity();

        accelerationXOld = planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_)[0];
        accelerationYOld = planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_)[1];

        x +=  step*vx + step*step/2.0* accelerationXOld;
        y +=  step*vy + step*step/2.0* accelerationYOld;

        planets_[planetNumber]->setXposition(x);
        planets_[planetNumber]->setYposition(y);

        accelerationX = planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_)[0];
        accelerationY = planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_)[1];

        vx +=  step/2.0*(accelerationXOld + accelerationX);
        vy +=  step/2.0*(accelerationYOld + accelerationY);

        planets_[planetNumber]->setXVelociy(vx);
        planets_[planetNumber]->setYVelociy(vy);

        planets_[planetNumber]->writeTofile(NAN, 0.0, 0.0);
    }
}

void Solver:: forwardEuler(vector<Planet *> planets_, int numberOfplanets_, int iterationStart)
{
    for (int planetNumber = iterationStart; planetNumber < numberOfplanets_; planetNumber++)
    {
        x = planets_[planetNumber]->getXPosition();
        y = planets_[planetNumber]->getYPosition();
        vx = planets_[planetNumber]->getXVelocity();
        vy = planets_[planetNumber]->getYVelocity();

        x += step*vx;
        y += step*vy;

        planets_[planetNumber]->setXposition(x);
        planets_[planetNumber]->setYposition(y);

        accelerationX = planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_)[0];
        accelerationY = planets_[planetNumber]->getAcceleration(planets_, numberOfplanets_)[1];

        vx += step*accelerationX;
        vy += step*accelerationY;

        planets_[planetNumber]->setXVelociy(vx);
        planets_[planetNumber]->setYVelociy(vy);

        planets_[planetNumber]->writeTofile(NAN, 0.0, 0.0);
    }
    accelerationX = 0.0;
    accelerationY = 0.0;
}*/
