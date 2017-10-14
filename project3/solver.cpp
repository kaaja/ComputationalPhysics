#include "solver.h"


ofstream ofile;

Solver:: Solver() { N = finalTime = 0; filename = "";}
Solver:: Solver(int N_, double finalTime_, string filename_)
{
    N = N_, finalTime = finalTime_;
    step = finalTime/double(N);
    time = 0.0;
    filename = filename_+string(".csv");
    ofile.open(filename);
    ofile << "time,x,y,vx/pi,vy/pi,potentialEnergy,kineticEnergy,angularMomentum,timeUsed,logTimeUsed,r" << endl;

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
    potentialEnergy = planet.getPotentialEnergy(r, mass);
    kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
    angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
    writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);

    start = clock();
    while (time < finalTime){
        x += step*vx;
        y += step*vy;
        vx -= step*FourPi2*x/(r*r*r);
        vy -= step*FourPi2*y/(r*r*r);
        r = sqrt(x*x + y*y);
        potentialEnergy = planet.getPotentialEnergy(r, mass);
        kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
        angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
        time += step;
        writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    //writeTofile(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, timeUsed);
    writeTofile(0, 0, 0, 0, 0, 0, 0, 0, timeUsed, 0);

    ofile.close();
}

void Solver:: velocityVerlet()
{
    pi = acos(-1.0);
    FourPi2 = 4.*pi*pi;

    start = clock();
    while (time < finalTime){

        x = planets[1].getXPosition();
        y = planets[1].getYPosition();
        vx = planets[1].getXVelocity();
        vy = planets[1].getYVelocity();
        mass = planets[1].getMass();
        r = planets[1].getRadialDistance(planets[0]);
        //r = sqrt(x*x + y*y);

        potentialEnergy = planets[1].getPotentialEnergy(r, mass);
        kineticEnergy   = planets[1].getKineticEnergy(mass, vx, vy);
        angularMomentum = planets[1].getAngularMomentum(r, mass, vx, vy);

        if (time==0.0)
            writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);
        planets[1].getForce(planets, &forceX, &forceY, numberOfPlanets);
        planets[1].getAcceleration(mass, &accelerationX, &accelerationY, forceX, forceY);

        accelerationXOld = accelerationX;
        accelerationYOld = accelerationY;

        x +=  step*vx + step*step/2.0* accelerationXOld;
        y +=  step*vy + step*step/2.0* accelerationYOld;

        planets[1].setXposition(x);
        planets[1].setYposition(y);

        r = planets[1].getRadialDistance(planets[0]);

        planets[1].getForce(planets, &forceX, &forceY, numberOfPlanets);
        planets[1].getAcceleration(mass, &accelerationX, &accelerationY, forceX, forceY);

        vx +=  step/2.0*(accelerationXOld + accelerationX);
        vy +=  step/2.0*(accelerationYOld + accelerationY);

        planets[1].setXVelociy(vx);
        planets[1].setYVelociy(vy);

        potentialEnergy = planets[1].getPotentialEnergy(r, mass);
        kineticEnergy   = planets[1].getKineticEnergy(mass, vx, vy);
        angularMomentum = planets[1].getAngularMomentum(r, mass, vx, vy);
        time += step;
        writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);


    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    //writeTofile(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, timeUsed);
    writeTofile(0, 0, 0, 0, 0, 0, 0, 0, timeUsed, 0);

    ofile.close();
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
    planet.getAcceleration(mass, &accelerationX, &accelerationY, forceX, forceY);

    potentialEnergy = planet.getPotentialEnergy(r, mass);
    kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
    angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
    writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);

    start = clock();
    while (time < finalTime){
        accelerationXOld = accelerationX;
        accelerationYOld = accelerationY;
        x +=  step*vx + step*step/2.0* accelerationXOld;
        y +=  step*vy + step*step/2.0* accelerationYOld;
        r = sqrt(x*x + y*y);
        planet.getAlternativeForce(mass, x, y, r, &forceX, &forceY, beta_);
        planet.getAcceleration(mass, &accelerationX, &accelerationY, forceX, forceY);
        vx +=  step/2.0*(accelerationXOld + accelerationX);
        vy +=  step/2.0*(accelerationYOld + accelerationY);
        potentialEnergy = planet.getPotentialEnergy(r, mass);
        kineticEnergy   = planet.getKineticEnergy(mass, vx, vy);
        angularMomentum = planet.getAngularMomentum(r, mass, vx, vy);
        time += step;
        writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy, angularMomentum, NAN, r);
    }
    finish = clock();
    timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));
    //writeTofile(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, timeUsed);
    writeTofile(0, 0, 0, 0, 0, 0, 0, 0, timeUsed, 0);

    ofile.close();
}
void Solver:: writeTofile(double time_, double x_, double y_, double vx_, double vy_, double potentialEnergy_, double kineticEnergy_ , double AngularMomentum_, double timeUsed_, double r_)
{
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(16) << time_ << ", ";
    ofile << setw(15) << setprecision(16) << x_ << ", ";
    ofile << setw(15) << setprecision(16) << y_ << ", ";
    ofile << setw(15) << setprecision(16) << vx_ << ", ";
    ofile << setw(15) << setprecision(16) << vy_ << ", ";
    ofile << setw(15) << setprecision(16) << potentialEnergy_ << ", ";
    ofile << setw(15) << setprecision(16) << kineticEnergy_ << ", ";
    ofile << setw(15) << setprecision(16) << AngularMomentum_ << ", ";
    ofile << setw(15) << setprecision(16) << timeUsed_ << ", ";
    ofile << setw(15) << setprecision(16) << log10(timeUsed_) << ", ";
    ofile << setw(15) << setprecision(16) << r_ << endl;
}

