#include <iostream>
#include "solver.h"

using namespace std;
void read_input (string& outfileName, double& dt, double& dx, double& theta, double& T, int argc, char** argv);
mat analyticalU (string outfileName, double dt, double dx, int Nt, int Nx);



int main(int argc, char* argv[])
{
    double dt, dx, theta, T;
    int Nx, Nt;
    double thetaForwardEuler = 0.0;
    double thetaBackwardEuler = 1.0;
    double thetaCranckNicholson = 0.5;
    string outfileName;

    read_input(outfileName, dx, dt, theta, T, argc, argv);

    Nt = (int) (round(T/dt)) + 1;
    Nx = (int) (round(1./dx+1));

    string oufileNameForwardEuler = outfileName + "ForwardEuler";
    string oufileNameBackwardEuler = outfileName + "BackwardEuler";
    string oufileNameCrancNicholson = outfileName + "CrancNicholson";

    Solver forwardEuler =  Solver( dt, dx, thetaForwardEuler, T, Nx, Nt);
    Solver backwardEuler =  Solver( dt, dx, thetaBackwardEuler, T,Nx, Nt);
    Solver cranckNicholson =  Solver( dt, dx, thetaCranckNicholson, T,Nx, Nt);

    forwardEuler.solve(oufileNameForwardEuler);
    backwardEuler.solve(oufileNameBackwardEuler);
    cranckNicholson.solve(oufileNameCrancNicholson);

    analyticalU(outfileName, dt, dx, Nt, Nx);
    return 0;
}

void read_input (string& outfileName, double &dt, double &dx, double &theta, double &T, int argc, char** argv)
{
    if( argc<= 1){
      cout << "Insert: outfile-name, number of simulations, amplification factor, start dimension" << endl;
      exit(1);
    }
    else{
      outfileName=argv[1];
    }
    dx = atof(argv[2]);
    dt = atof(argv[3]);
    theta = atof(argv[4]);
    T = atof(argv[5]);
}

mat analyticalU(string outfileName, double dt, double dx, int Nt, int Nx)
{
   mat analyticalMatrixU = zeros<mat>(Nx,Nt);
   string outfileNameAnalytical;
   double uANalytical;
   for (int t = 1; t < Nt; t++){
        for(int xj = 0; xj < Nx; xj++){
            uANalytical = 0.;
            for (int k = 1; k < Nx*100; k++){
                uANalytical += exp(-pow(k*M_PI, 2)*(dt*t))/M_PI*pow((-1.), k)*sin(k*M_PI*dx*xj);
            }
            uANalytical += dx*xj;
            analyticalMatrixU(xj, t) = uANalytical;
        }
   }
    outfileNameAnalytical = outfileName + "AnalyticalMatrix";
    analyticalMatrixU.save(outfileNameAnalytical , raw_ascii);
    return analyticalMatrixU;

}
