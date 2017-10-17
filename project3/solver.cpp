#include "solver.h"


//ofstream ofile;

Solver:: Solver() { N = finalTime = 0; filename = "";}
Solver:: Solver(int N_, double finalTime_, string filename_)
{
    N = N_, finalTime = finalTime_;
    step = finalTime/double(N);
    time = 0.0;
    centerOfMassSystem = "False";
    /*filename = filename_+string(".csv");
    ofile.open(filename);
    ofile << "time,x,y,vx/pi,vy/pi,potentialEnergy,kineticEnergy,angularMomentum,timeUsed,logTimeUsed,r" << endl;
    */

}

void Solver:: addPlanet(Planet planet_)
{
    planet = planet_;
    planets.push_back (planet);
    numberOfPlanets += 1;
}

void Solver:: forwardEuler()
{
    x = planet.getXPosition();
    y = planet.getYPosition();
    vx = planet.getXVelocity();
    vy = planet.getYVelocity();
    mass = planet.getMass();
    r = sqrt(x*x + y*y); //r = planet.getDistance()//
    pi = acos(-1.0);
    FourPi2 = 4.*pi*pi;
    //potentialEnergy = planet.getPotentialEnergy(r, mass);
    //kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
    //angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
    //writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);

    start = clock();
    while (time < finalTime){
        x += step*vx;
        y += step*vy;
        vx -= step*FourPi2*x/(r*r*r);
        vy -= step*FourPi2*y/(r*r*r);
        r = sqrt(x*x + y*y);
        //potentialEnergy = planet.getPotentialEnergy(r, mass);
        //kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
        //angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
        time += step;
        //writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    //writeTofile(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, timeUsed);
    //writeTofile(0, 0, 0, 0, 0, 0, 0, 0, timeUsed, 0);

    //ofile.close();
}

void Solver:: velocityVerlet()
{
    pi = acos(-1.0);
    FourPi2 = 4.*pi*pi;

    start = clock();

    int iterationStart = 1;

    /*if (centerOfMassSystem == "True")
        iterationStart = 0;*/

    while (time < finalTime){
        for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        {
            x = planets[planetNumber].getXPosition();
            y = planets[planetNumber].getYPosition();
            vx = planets[planetNumber].getXVelocity();
            vy = planets[planetNumber].getYVelocity();
            mass = planets[planetNumber].getMass();

            if (time==0.0){
                planets[planetNumber].writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
            }
            planets[planetNumber].getAcceleration(planets, &accelerationX, &accelerationY, numberOfPlanets);

            accelerationXOld = accelerationX;
            accelerationYOld = accelerationY;

            x +=  step*vx + step*step/2.0* accelerationXOld;
            y +=  step*vy + step*step/2.0* accelerationYOld;

            planets[planetNumber].setXposition(x);
            planets[planetNumber].setYposition(y);

            planets[planetNumber].getAcceleration(planets, &accelerationX, &accelerationY, numberOfPlanets);

            vx +=  step/2.0*(accelerationXOld + accelerationX);
            vy +=  step/2.0*(accelerationYOld + accelerationY);

            planets[planetNumber].setXVelociy(vx);
            planets[planetNumber].setYVelociy(vy);

            planets[planetNumber].setTime(time+step);
            planets[planetNumber].writeTofile(NAN, getCenterOfMassX(), getCenterOfMassY());
        }
        time += step;
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    for (int planetNumber = iterationStart; planetNumber < numberOfPlanets; planetNumber++)
        planets[planetNumber].writeTofile(timeUsed, getCenterOfMassX(), getCenterOfMassY());
}

void Solver:: alternativeForceVelocityVerlet(double beta_)
{

    pi = acos(-1.0);
    FourPi2 = 4.*pi*pi;

    x = planet.getXPosition();
    y = planet.getYPosition();
    vx = planet.getXVelocity();
    vy = planet.getYVelocity();
    mass = planet.getMass();
    r = sqrt(x*x + y*y);
    planet.getAlternativeForce(mass, x, y, r, &forceX, &forceY, beta_);
    //planet.getAcceleration(mass, &accelerationX, &accelerationY, forceX, forceY);

    //potentialEnergy = planet.getPotentialEnergy(r, mass);
    //kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
    //angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
    //writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);

    start = clock();
    while (time < finalTime){
        accelerationXOld = accelerationX;
        accelerationYOld = accelerationY;
        x +=  step*vx + step*step/2.0* accelerationXOld;
        y +=  step*vy + step*step/2.0* accelerationYOld;
        r = sqrt(x*x + y*y);
        planet.getAlternativeForce(mass, x, y, r, &forceX, &forceY, beta_);
        //planet.getAcceleration2(mass, &accelerationX, &accelerationY, forceX, forceY);
        vx +=  step/2.0*(accelerationXOld + accelerationX);
        vy +=  step/2.0*(accelerationYOld + accelerationY);
        //potentialEnergy = planet.getPotentialEnergy(r, mass);
        //kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
        //angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
        time += step;
        //writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    //writeTofile(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, timeUsed);
    //writeTofile(0, 0, 0, 0, 0, 0, 0, 0, timeUsed, 0);

    //ofile.close();
}

double Solver:: getCenterOfMassX()
{
    double centerOfMassX;
    double totalMass;
    for (int planetNumber = 0; planetNumber < numberOfPlanets; planetNumber++)
    {
        totalMass += planets[planetNumber].getMass();
        centerOfMassX += planets[planetNumber].getXPosition()*planets[planetNumber].getMass();
    }
    return centerOfMassX/totalMass;
}

double Solver::getCenterOfMassY()
{
    double centerOfMassY;
    double totalMass;
    for ( int planetNumber = 0; planetNumber < numberOfPlanets; planetNumber++)
    {
        totalMass += planets[planetNumber].getMass();
        centerOfMassY += planets[planetNumber].getYPosition()*planets[planetNumber].getMass();
    }
    return centerOfMassY/totalMass;
}

void Solver::changeToCenterOfMassSystem()
{
    double centerOfMassX = getCenterOfMassX();
    double centerOfMassY = getCenterOfMassY();
    for (int planetNumber = 0; planetNumber < numberOfPlanets; planetNumber++)
    {
        double xPosition = planets[planetNumber].getXPosition();
        double yPosition = planets[planetNumber].getYPosition();
        planets[planetNumber].setXposition(xPosition - centerOfMassX);
        planets[planetNumber].setYposition(yPosition - centerOfMassY);
        planets[planetNumber].changeToCenterOfMassSystemInAcceleration();
    }
    setSunVelocity();
    centerOfMassSystem = "True";

}

void Solver::setSunVelocity()
{
    double momentumOfPlanetsX, momentumOfPlanetsY;
    for (int planetNumber = 1; planetNumber < numberOfPlanets; planetNumber++)
    {
        momentumOfPlanetsX += planets[planetNumber].getMass()*planets[planetNumber].getXVelocity();
        momentumOfPlanetsY += planets[planetNumber].getMass()*planets[planetNumber].getYVelocity();
    }
    planets[0].setXVelociy(-momentumOfPlanetsX/planets[0].getMass());
    planets[0].setYVelociy(-momentumOfPlanetsY/planets[0].getMass());
}
