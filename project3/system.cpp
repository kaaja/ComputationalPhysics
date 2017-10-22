#include "system.h"

System:: System()
{
    centerOfMassSystem = false;
}

void System:: addPlanet(Planet &planet_)
{
    planet = &planet_;
    planets.push_back (planet);
    numberOfPlanets += 1;

}


void System:: simulate(Solver &solver_, int N_, double finalTime_)
{
    N = N_, finalTime = finalTime_;
    step = finalTime/double(N);
    time = 0.0;
    int iterationStart = 1;

    for(int i =0;i<numberOfPlanets;i++)
    {
        planet->setStep(step);
        planet->writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
    }
    if (centerOfMassSystem) iterationStart = 0;

    start = clock();
    while (time < finalTime){
        for(int planetNumber=iterationStart;planetNumber<numberOfPlanets;planetNumber++) planets[planetNumber]->setTime(time+step);
        solver_.solve(step, planets, numberOfPlanets, iterationStart);
        time += step;
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        planets[planetNumber]->writeTofile(timeUsed, getCenterOfMassX(), getCenterOfMassY());
}

double System:: getCenterOfMassX()
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

double System::getCenterOfMassY()
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

void System::changeToCenterOfMassSystem()
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

void System::setSunVelocity()
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
