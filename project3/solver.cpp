#include "solver.h"


//ofstream ofile;

Solver:: Solver() { N = finalTime = 0; filename = "";}
Solver:: Solver(int N_, double finalTime_)
{
    N = N_, finalTime = finalTime_;
    step = finalTime/double(N);
    time = 0.0;
    centerOfMassSystem = false;

}

void Solver:: addPlanet(Planet &planet_)
{
    planet = &planet_;
    planets.push_back (planet);
    numberOfPlanets += 1;
    planet->setStep(step);
}


void Solver:: forwardEuler()
{
    start = clock();

    int iterationStart = 1; // If fixed sun

    if (centerOfMassSystem) // Moving sun
    {
        iterationStart = 0;
    }
    // Write initial values
    for(int planetNumber=iterationStart;planetNumber<numberOfPlanets;planetNumber++){
        planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
    }
    while (time < finalTime){
        //  Accelerations vector all planets
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            acc= planets[planetNumber]->getAcceleration(planets, numberOfPlanets);
            if(iterationStart==1 && planetNumber==iterationStart){
                accXVec.push_back(0.0);
                accYVec.push_back(0.0);
            }
            accXVec.push_back(acc[0]);
            accYVec.push_back(acc[1]);
            acc.clear();
        }
        // Position
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            x = planets[planetNumber]->getXPosition();
            y = planets[planetNumber]->getYPosition();
            vx = planets[planetNumber]->getXVelocity();
            vy = planets[planetNumber]->getYVelocity();

            x +=  step*vx;
            y +=  step*vy;
            planets[planetNumber]->setXposition(x);
            planets[planetNumber]->setYposition(y);
        }
        // Velocity
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            vx = planets[planetNumber]->getXVelocity();
            vy = planets[planetNumber]->getYVelocity();

            vx +=  step*accXVec[planetNumber];
            vx +=  step*accYVec[planetNumber];
            planets[planetNumber]->setXVelociy(vx);
            planets[planetNumber]->setYVelociy(vy);
        }
        // Reset acceleration vectors
        accXVec.clear();
        accYVec.clear();
        // Write results
        for(int planetNumber=iterationStart;planetNumber<numberOfPlanets;planetNumber++){
            planets[planetNumber]->setTime(time+step);
            planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
        }
        time += step;
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    // Write time used only
    for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        planets[planetNumber]->writeTofile(timeUsed, getCenterOfMassX(), getCenterOfMassY());
}

void Solver:: velocityVerlet()
{
    start = clock();

    int iterationStart = 1; // If fixed sun

    if (centerOfMassSystem) // Moving sun
    {
        iterationStart = 0;
    }
    // Write initial values
    for(int planetNumber=iterationStart;planetNumber<numberOfPlanets;planetNumber++){
        planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
    }
    while (time < finalTime){
        // Vector accelerations all planets
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            acc= planets[planetNumber]->getAcceleration(planets, numberOfPlanets);
            if(iterationStart==1 && planetNumber==iterationStart){
                accXVec.push_back(0.0);
                accYVec.push_back(0.0);
            }
            accXVec.push_back(acc[0]);
            accYVec.push_back(acc[1]);
            acc.clear();
        }
        // Position
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            x = planets[planetNumber]->getXPosition();
            y = planets[planetNumber]->getYPosition();
            vx = planets[planetNumber]->getXVelocity();
            vy = planets[planetNumber]->getYVelocity();

            x +=  step*vx + step*step/2.0* accXVec[planetNumber];
            y +=  step*vy + step*step/2.0* accYVec[planetNumber];
            planets[planetNumber]->setXposition(x);
            planets[planetNumber]->setYposition(y);
        }
        // Vector old accelerations
        accXVecOld.clear();
        accYVecOld.clear();
        accXVecOld = accXVec;
        accYVecOld = accYVec;
        accXVec.clear();
        accYVec.clear();
        // Acceleration update
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            acc= planets[planetNumber]->getAcceleration(planets, numberOfPlanets);
            if(iterationStart==1 && planetNumber==iterationStart){
                accXVec.push_back(0.0);
                accYVec.push_back(0.0);
            }
            accXVec.push_back(acc[0]);
            accYVec.push_back(acc[1]);
            acc.clear();
        }
        // Velocity
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            vx = planets[planetNumber]->getXVelocity();
            vy = planets[planetNumber]->getYVelocity();

            vx +=  step/2.0*(accXVecOld[planetNumber] + accXVec[planetNumber]);
            vy +=  step/2.0*(accYVecOld[planetNumber] + accYVec[planetNumber]);
            planets[planetNumber]->setXVelociy(vx);
            planets[planetNumber]->setYVelociy(vy);
        }
        // Write results
        for(int planetNumber=iterationStart;planetNumber<numberOfPlanets;planetNumber++){
            planets[planetNumber]->setTime(time+step);
            planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
        }
        time += step;
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    // Write time used only
    for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        planets[planetNumber]->writeTofile(timeUsed, getCenterOfMassX(), getCenterOfMassY());
}

void Solver:: velocityVerletOld()
{
    start = clock();

    int iterationStart = 1;

    if (centerOfMassSystem)
    {
        iterationStart = 0;
    }
    while (time < finalTime){
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            acc = planets[planetNumber]->getAcceleration(planets, numberOfPlanets);
            if(iterationStart==1 && planetNumber==iterationStart){
                accelerationsX.push_back(0.0);
                accelerationsY.push_back(0.0);
            }
            accelerationsX.push_back(acc[0]);
            accelerationsY.push_back(acc[1]);
            acc.clear();
        }
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            x = planets[planetNumber]->getXPosition();
            y = planets[planetNumber]->getYPosition();
            vx = planets[planetNumber]->getXVelocity();
            vy = planets[planetNumber]->getYVelocity();

            if (time==0.0){
                planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
            }

            accelerationXOld = accelerationsX[planetNumber];
            accelerationYOld = accelerationsY[planetNumber];

            x +=  step*vx + step*step/2.0* accelerationXOld;
            y +=  step*vy + step*step/2.0* accelerationYOld;

            planets[planetNumber]->setXposition(x);
            planets[planetNumber]->setYposition(y);

            accelerationX = planets[planetNumber]->getAcceleration(planets, numberOfPlanets)[0];
            accelerationY = planets[planetNumber]->getAcceleration(planets, numberOfPlanets)[1];

            vx +=  step/2.0*(accelerationXOld + accelerationX);
            vy +=  step/2.0*(accelerationYOld + accelerationY);

            planets[planetNumber]->setXVelociy(vx);
            planets[planetNumber]->setYVelociy(vy);

            planets[planetNumber]->setTime(time+step);
            planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
        }
        accelerationsX.clear();
        accelerationsY.clear();
        time += step;
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        planets[planetNumber]->writeTofile(timeUsed, getCenterOfMassX(), getCenterOfMassY());
}

void Solver:: forwardEulerOld()
{
    int iterationStart = 1;

    if (centerOfMassSystem)
    {
        iterationStart = 0;
    }
    start = clock();
    while (time < finalTime){
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            if (time==0.0 )
                planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());

            x = planets[planetNumber]->getXPosition();
            y = planets[planetNumber]->getYPosition();
            vx = planets[planetNumber]->getXVelocity();
            vy = planets[planetNumber]->getYVelocity();


            x += step*vx;
            y += step*vy;

            planets[planetNumber]->setXposition(x);
            planets[planetNumber]->setYposition(y);

            accelerationX = planets[planetNumber]->getAcceleration(planets, numberOfPlanets)[0];
            accelerationY = planets[planetNumber]->getAcceleration(planets, numberOfPlanets)[1];

            vx += step*accelerationX;
            vy += step*accelerationY;

            planets[planetNumber]->setXVelociy(vx);
            planets[planetNumber]->setYVelociy(vy);

            planets[planetNumber]->setTime(time+step);
            planets[planetNumber]->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
        }
        accelerationX = 0.0;
        accelerationY = 0.0;
        time += step;
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        planets[planetNumber]->writeTofile(timeUsed, getCenterOfMassX(), getCenterOfMassY());
}

double Solver:: getCenterOfMassX()
{
    double centerOfMassX = 0.;
    double totalMass = 0.;
    for (int planetNumber = 0; planetNumber < numberOfPlanets; planetNumber++)
    {
        totalMass += planets[planetNumber]->getMass();
        centerOfMassX += planets[planetNumber]->getXPosition()*planets[planetNumber]->getMass();
    }
    return centerOfMassX/totalMass;
}

double Solver::getCenterOfMassY()
{
    double centerOfMassY = 0.;
    double totalMass = 0.;
    for ( int planetNumber = 0; planetNumber < numberOfPlanets; planetNumber++)
    {
        totalMass += planets[planetNumber]->getMass();
        centerOfMassY += planets[planetNumber]->getYPosition()*planets[planetNumber]->getMass();
    }
    return centerOfMassY/totalMass;
}

void Solver::changeToCenterOfMassSystem()
{
    double centerOfMassX = getCenterOfMassX();
    double centerOfMassY = getCenterOfMassY();
    for (int planetNumber = 0; planetNumber < numberOfPlanets; planetNumber++)
    {
        double xPosition = planets[planetNumber]->getXPosition();
        double yPosition = planets[planetNumber]->getYPosition();
        planets[planetNumber]->setXposition(xPosition - centerOfMassX);
        planets[planetNumber]->setYposition(yPosition - centerOfMassY);
    }
    setSunVelocity();
    centerOfMassSystem = true;
}

void Solver::setSunVelocity()
{
    double momentumOfPlanetsX = 0.;
    double momentumOfPlanetsY = 0.;
    for (int planetNumber = 1; planetNumber < numberOfPlanets; planetNumber++)
    {
        momentumOfPlanetsX += planets[planetNumber]->getMass()*planets[planetNumber]->getXVelocity();
        momentumOfPlanetsY += planets[planetNumber]->getMass()*planets[planetNumber]->getYVelocity();
    }
    planets[0]->setXVelociy(-momentumOfPlanetsX/planets[0]->getMass());
    planets[0]->setYVelociy(-momentumOfPlanetsY/planets[0]->getMass());
}
