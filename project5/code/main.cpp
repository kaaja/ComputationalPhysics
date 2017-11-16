#include <iostream>
#include "solver.h"

using namespace std;
void read_input (string& outfileName, double& dt, double& dx, double& theta, double& T, int argc, char** argv);
mat analyticalU (string outfileName, double dt, double dx, int Nt, int Nx);
void output_scalars(double computed_error, double dt, double dx);

ofstream ofile1; // File for scalars


int main(int argc, char* argv[])
{
    double dt, dx, theta, T;
    int Nx, Nt;
    double thetaForwardEuler = 0.0;
    double thetaBackwardEuler = 1.0;
    double thetaCranckNicholson = 0.5;
    double computed_error;
    string outfileName, outfileNameNorms;

    read_input(outfileName, dx, dt, theta, T, argc, argv);

    Nt = (int) (round(T/dt)) + 1;
    Nx = (int) (round(1./dx+1));

    mat solutionMatrixFe = zeros<mat>(Nx+1,Nt+1);
    mat solutionMatrixBe = zeros<mat>(Nx+1,Nt+1);
    mat solutionMatrixCn = zeros<mat>(Nx+1,Nt+1);
    mat solutionMatrixExact = zeros<mat>(Nx+1,Nt+1);
    int numberOfSchemes = 3;
    int numberOfVariablesNormMatrix = 3;
    mat normMatrix= zeros<mat>(numberOfSchemes, numberOfVariablesNormMatrix);

    string oufileNameForwardEuler = outfileName + "ForwardEuler";
    string oufileNameBackwardEuler = outfileName + "BackwardEuler";
    string oufileNameCrancNicholson = outfileName + "CrancNicholson";

    Solver forwardEuler =  Solver( dt, dx, thetaForwardEuler, T, Nx, Nt);
    Solver backwardEuler =  Solver( dt, dx, thetaBackwardEuler, T,Nx, Nt);
    Solver cranckNicholson =  Solver( dt, dx, thetaCranckNicholson, T,Nx, Nt);

    solutionMatrixFe  = forwardEuler.solve(oufileNameForwardEuler);
    solutionMatrixBe  = backwardEuler.solve(oufileNameBackwardEuler);
    solutionMatrixCn  = cranckNicholson.solve(oufileNameCrancNicholson);

    solutionMatrixExact = analyticalU(outfileName, dt, dx, Nt, Nx);

    forwardEuler.calculate_error(solutionMatrixFe, solutionMatrixExact, &computed_error, Nx, Nt);
    normMatrix(0,0) = log10(dt);
    normMatrix(0,1) = log10(dx);
    normMatrix(0,2) = computed_error;

    backwardEuler.calculate_error(solutionMatrixBe, solutionMatrixExact, &computed_error, Nx, Nt);
    normMatrix(1,0) = log10(dt);
    normMatrix(1,1) = log10(dx);
    normMatrix(1,2) = computed_error;

    cranckNicholson.calculate_error(solutionMatrixCn, solutionMatrixExact, &computed_error, Nx, Nt);
    normMatrix(2,0) = log10(dt);
    normMatrix(2,1) = log10(dx);
    normMatrix(2,2) = computed_error;

    normMatrix.save(outfileName + ".txt", raw_ascii);
    /*
    outfileNameNorms = outfileName + string("Norms")+string(".csv");
    ofile1.open(outfileNameNorms);
    ofile1 << "dt,dx,computedError" << endl;
    output_scalars(computed_error, dt, dx);
    */
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
   mat analyticalMatrixU = zeros<mat>(Nx+1,Nt+1);
   string outfileNameAnalytical;
   double uANalytical;
   for (int t = 1; t < Nt; t++){
        for(int xj = 0; xj < Nx; xj++){
            uANalytical = 0.;
            for (int k = 1; k < Nx; k++){
                uANalytical += exp(-pow(k*M_PI, 2)*(dt*t))*2./(k*M_PI)*pow((-1.), k)*sin(k*M_PI*dx*xj);
            }
            uANalytical += dx*xj;
            analyticalMatrixU(0, t+1) = t*dt;
            analyticalMatrixU(xj+1, 0) = xj*dx;
            analyticalMatrixU(xj+1, t+1) = uANalytical;
        }
   }
    outfileNameAnalytical = outfileName + "AnalyticalSolutionMatrixU.txt";
    analyticalMatrixU.save(outfileNameAnalytical , raw_ascii);
    return analyticalMatrixU;
}

void output_scalars(double computed_error, double dt, double dx){
  ofile1 << setiosflags(ios::showpoint | ios::uppercase);
  ofile1 << setw(15) << setprecision(16) << log10(dt) << ", ";
  ofile1 << setw(15) << setprecision(16) << log10(dx) << ", ";
  ofile1 << setw(15) << setprecision(16) << computed_error << endl;
}
