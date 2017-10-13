#include "planet.h"
#include "solver.h"

using namespace std;

void initialize ( string& outfileName, double& finalTime, int& N, string& solverType, double &initialVy, double &beta, int argc, char** argv);

int main(int argc, char* argv[])
{
    int N;
    double finalTime, initialVy, beta;
    double initialVyJupiter;
    string outfileName, solverType;

    initialize( outfileName, finalTime, N, solverType, initialVy, beta, argc, argv);

    double pi = acos(-1.0);
    Planet sun(1., 0., 0., 0. , 0.);
    Planet earth(0.000003, 1., 0., 0. , initialVy);
    double velocityEarth = 29.7559;
    double velocityJupiter = 13.0697;
    double velocityJupiterToEarth = velocityJupiter/velocityEarth;
    //initialVyJupiter = 2*pi*velocityJupiterToEarth;
    //Planet jupiter(9.5e-4, 5.2, 0., 0., initialVyJupiter);
    Solver solution(N, finalTime, outfileName);
    solution.addPlanet(sun);
    solution.addPlanet(earth);

    //solution.addPlanet(jupiter);
    //cout << earth.getInitialXPosition() << forwardEuler.
    if (solverType == "ForwardEuler")
        solution.forwardEuler();
    else if (solverType == "VelocityVerlet")
        solution.velocityVerlet();
    else
        solution.alternativeForceVelocityVerlet(beta);
    return 0;
}

void initialize (string& outfileName, double& finalTime, int& N, string& solverType, double &initialVy, double &beta, int argc, char** argv)
{
    if( argc<= 1){
      cout << "Insert: outfile-name, number of simulations, amplification factor, start dimension" << endl;
      exit(1);
    }
    else{
      outfileName=argv[1];
    }
    finalTime = atof(argv[2]);
    N = atoi(argv[3]);
    solverType = argv[4];
    initialVy = atof(argv[5]);
    beta = atof(argv[6]);
}
