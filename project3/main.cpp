#include "planet.h"
#include "solver.h"

using namespace std;

void initialize ( string& outfileName, double& finalTime, int& N, string& solverType, int argc, char** argv);

int main(int argc, char* argv[])
{
    int N;
    double finalTime;
    string outfileName, solverType;

    initialize( outfileName, finalTime, N, solverType, argc, argv);

    double pi = acos(-1.0);
    Planet earth(0.000003, 1., 0., 0. ,2*pi);
    Solver solution(N, finalTime, outfileName);
    solution.addPlanet(earth);
    //cout << earth.getInitialXPosition() << forwardEuler.
    solution.forwardEuler();
    return 0;
}

void initialize ( string& outfileName, double& finalTime, int& N, string& solverType, int argc, char** argv)
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
}
