#include "planet.h"

ofstream ofile;

Planet:: Planet () { mass = xPosition = yPosition = xVelocity = yVelocity = 0.0;}

Planet:: Planet (double mass_, double xPosition_, double yPosition_, double xVelocity_, double yVelocity_, string filename_, string planetName_)
{

    mass = mass_;
    xVelocity = xVelocity_;
    yVelocity = yVelocity_;
    xPosition = xPosition_;
    yPosition = yPosition_;
    time      = 0.0;
    radialDistance = sqrt(xPosition*xPosition + yPosition*yPosition);
    //ofstream ofile;
    filename = filename_+ planetName_ + string(".csv");
    ofile.open(filename);
    ofile << "time,x,y,vx/pi,vy/pi,potentialEnergy,kineticEnergy,angularMomentum,timeUsed,logTimeUsed,r" << endl;
    ofile.close();

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
double Planet:: getKineticEnergy()
{
    double velocity2 = xVelocity*xVelocity + yVelocity*yVelocity;
    return 0.5*mass*velocity2;
}
double Planet:: getPotentialEnergy(double r_)
{
    //return 6.67*pow(10,-11)*mass/r_;
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    return FourPi2*mass/r_;
}
double Planet:: getAngularMomentum(double r_)
{
    double velocity = sqrt(xVelocity*xVelocity + yVelocity*yVelocity);
    return r_*mass*velocity;
}
void Planet:: getAcceleration(vector<Planet> planets_, double *accelerationX_, double *accelerationY_, int numberOfPlanets_)
{
    double pi = acos(-1.0);
    double FourPi2 = 4.*pi*pi;
    double r;
    double forceX, forceY;
    forceX = -FourPi2*xPosition*mass/getRadialDistance(planets_[0]);
    forceY = -FourPi2*yPosition*mass/getRadialDistance(planets_[0]);
    for (int planetNumber = 1; planetNumber < numberOfPlanets_; planetNumber++)
    {
        if (planets_[planetNumber].getMass() != mass){
            r = getRadialDistance(planets_[planetNumber]);
            forceX += FourPi2*planets_[planetNumber].getMass()/planets_[0].getMass()*(xPosition - planets_[planetNumber].getXPosition())/pow(r,3);
            forceY += FourPi2*planets_[planetNumber].getMass()/planets_[0].getMass()*(yPosition - planets_[planetNumber].getYPosition())/pow(r,3);
        }
    }
    *accelerationX_ = forceX/mass;
    *accelerationY_ = forceY/mass;
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
void Planet:: setTime(double time_){time = time_;}

void Planet:: writeTofile(double timeUsed_, double r_)
{
    //ofs.open ("test.txt", std::ofstream::out | std::ofstream::app);
    ofile.open(filename, ios::out | ios::app);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(16) << time << ", ";
    ofile << setw(15) << setprecision(16) << xPosition << ", ";
    ofile << setw(15) << setprecision(16) << yPosition << ", ";
    ofile << setw(15) << setprecision(16) << xVelocity << ", ";
    ofile << setw(15) << setprecision(16) << yVelocity << ", ";
    ofile << setw(15) << setprecision(16) << getPotentialEnergy(r_) << ", ";
    ofile << setw(15) << setprecision(16) << getKineticEnergy() << ", ";
    ofile << setw(15) << setprecision(16) << getAngularMomentum(r_) << ", ";
    ofile << setw(15) << setprecision(16) << timeUsed_ << ", ";
    ofile << setw(15) << setprecision(16) << log10(timeUsed_) << ", ";
    ofile << setw(15) << setprecision(16) << r_ << endl;
    /*if (timeUsed_ != NAN)
        ofile.close();*/
    ofile.close();
}
