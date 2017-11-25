#include <iostream>
#include "solver.h"
#include "twodimensionaldiffusionsolver.h"
#include "lib.h"
#include <omp.h>

using namespace std;
void read_input (string& outfileName, double& dt, double& dx, double& theta, double& T, string& scenario, int& threadNumberFromUser, int argc, char** argv);
mat analyticalU (string outfileName, double dt, double dx, int Nt, int Nx);
void analytical2D(string outfileName, double dt, double dx, double dy, int Nt, int Nx, int Ny, TwoDimensionalDiffusionSolver solver);
void output_scalars(double computed_error, double dt, double dx);
double gaussQuad(int numberOfSummationPoints, int n, int m, int numberOfIntegrationPoints);

ofstream ofile1; // File for scalars


int main(int argc, char* argv[])
{
    double dt, dx, dy, theta, T, wtime, wtime2, wtime3;
    int Nx, Ny, Nt, threadNumberFromUser;
    double thetaForwardEuler = 0.0;
    double thetaBackwardEuler = 1.0;
    double thetaCranckNicholson = 0.5;
    double computed_error = 0.0;
    string outfileName, outfileNameNorms, scenario;

    read_input(outfileName, dx, dt, theta, T, scenario, threadNumberFromUser, argc, argv);

    int thread_num;
    omp_set_num_threads(threadNumberFromUser);
    thread_num = omp_get_max_threads ();
    cout << "  The number of processors available = " << omp_get_num_procs () << endl ;
    cout << "  The number of threads available    = " << thread_num <<  endl;


    Nt = (int) (round(T/dt)) + 1;
    Nx = (int) (round(1./dx+1));
    Ny = Nx;
    dy = dx;

    if (scenario == "1D"){
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

        cout << "dt: " << dt << endl;

        Solver forwardEuler =  Solver( dt, dx, thetaForwardEuler, T, Nx, Nt);
        Solver backwardEuler =  Solver( dt, dx, thetaBackwardEuler, T,Nx, Nt);
        Solver cranckNicholson =  Solver( dt, dx, thetaCranckNicholson, T,Nx, Nt);

        solutionMatrixFe  = forwardEuler.solve(oufileNameForwardEuler);
        solutionMatrixBe  = backwardEuler.solve(oufileNameBackwardEuler);
        solutionMatrixCn  = cranckNicholson.solve(oufileNameCrancNicholson);

        solutionMatrixExact = analyticalU(outfileName, dt, dx, Nt, Nx);

        forwardEuler.calculate_error(solutionMatrixFe, solutionMatrixExact, &computed_error, Nx, Nt);
        normMatrix(0,0) = log2(dt);
        normMatrix(0,1) = log2(dx);
        normMatrix(0,2) = computed_error;

        backwardEuler.calculate_error(solutionMatrixBe, solutionMatrixExact, &computed_error, Nx, Nt);
        normMatrix(1,0) = log2(dt);
        normMatrix(1,1) = log2(dx);
        normMatrix(1,2) = computed_error;

        cranckNicholson.calculate_error(solutionMatrixCn, solutionMatrixExact, &computed_error, Nx, Nt);
        normMatrix(2,0) = log2(dt);
        normMatrix(2,1) = log2(dx);
        normMatrix(2,2) = computed_error;

        normMatrix.save(outfileName + ".txt", raw_ascii);
    }
    else if (scenario == "2D"){
        string outfileName2DExplicit;
        string outfileName2DImplicit;

        TwoDimensionalDiffusionSolver explicit2D = TwoDimensionalDiffusionSolver( dt, dx, dy, thetaForwardEuler, T, Nx, Ny, Nt);
        TwoDimensionalDiffusionSolver implicit2D = TwoDimensionalDiffusionSolver( dt, dx, dy, thetaForwardEuler, T, Nx, Ny, Nt);

        outfileName2DExplicit = outfileName + "Explicit";
        wtime = omp_get_wtime ( );
        explicit2D.solve(outfileName2DExplicit, "Explicit", threadNumberFromUser);
        wtime = omp_get_wtime ( ) - wtime;
        cout << "Time used explicit: " << wtime << endl;

        outfileName2DImplicit = outfileName + "Implicit";
        wtime2 = omp_get_wtime ( );
        implicit2D.solve(outfileName2DImplicit, "Implicit", threadNumberFromUser);
        wtime2 = omp_get_wtime ( ) - wtime2;
        cout << "Time used implicit: " << wtime2 << endl;

        wtime3 = omp_get_wtime ( );
        analytical2D(outfileName, dt, dx, dy, Nt, Nx, Ny, explicit2D);
        wtime3 = omp_get_wtime ( ) - wtime3;
        cout << "Time used analytic: " << wtime3 << endl;


        ofile1.open(outfileName + "Timing.txt");
        ofile1 << "explicit,implicit,analytic" << endl;
        ofile1 << setiosflags(ios::showpoint | ios::uppercase);
        ofile1 << setw(15) << setprecision(16) << wtime<< ", ";
        ofile1 << setw(15) << setprecision(16) << wtime2<< ", ";
        ofile1 << setw(15) << setprecision(16) << wtime3<< endl;
    }
    else if (scenario == "2DJacobiIterations"){
        string outfileName2DImplicit;
        TwoDimensionalDiffusionSolver implicit2D = TwoDimensionalDiffusionSolver( dt, dx, dy, thetaForwardEuler, T, Nx, Ny, Nt);

        outfileName2DImplicit = outfileName ;
        wtime2 = omp_get_wtime ( );
        implicit2D.solve(outfileName2DImplicit, "Implicit", threadNumberFromUser);
        wtime2 = omp_get_wtime ( ) - wtime2;
        cout << "Time used implicit: " << wtime2 << endl;
    }
    return 0;
}

void read_input (string& outfileName, double& dt, double& dx, double& theta, double& T, string& scenario, int& threadNumberFromUser, int argc, char** argv)
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
    scenario = argv[6];
    threadNumberFromUser = atoi(argv[7]);
}

mat analyticalU(string outfileName, double dt, double dx, int Nt, int Nx)
{
   mat analyticalMatrixU = zeros<mat>(Nx+1,Nt+1);
   string outfileNameAnalytical;
   double uANalytical;
   for (int t = 1; t < Nt; t++){
        for(int xj = 0; xj < Nx; xj++){
            uANalytical = 0.;
            for (int k = 1; k < 100; k++){
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
void analytical2D(string outfileName, double dt, double dx, double dy, int Nt, int Nx, int Ny, TwoDimensionalDiffusionSolver solver)
{
    int t=0;int i=0;int j=0;int n=0;int m=0;
    string outfileNameAnalytical = outfileName + "AnalyticalSolutionMatrixU2D";
    double uSs;
    mat analyticalMatrixU2D;
    double fourOverPi = 4.0/M_PI;
    double tempSum = 0.0;
    int counter = 1;
    uSs = 0.0;
    double integral = 0.;
    int integrationPoints = 30;
    int sumLimit = 15;// for sums with integrals
    for (t = 0; t < Nt; t++){
        analyticalMatrixU2D = zeros<mat>(Nx,Ny);
        #pragma omp parallel for default(shared) private(i,j,n,m) reduction(+:tempSum)
        for (i = 0; i < Nx ; i++){//note: changed from 1 to 0
            for (j = 0; j < Ny; j++){
                #pragma omp critical
                tempSum = 0.;
                //#pragma omp parallel for default(shared) private(n,m) reduction(+:tempSum)
                for (n = 1; n < sumLimit; n++){
                    for (m = 1; m < sumLimit; m++){
                        integral = gaussQuad(sumLimit, n, m, integrationPoints);
                        tempSum += integral*sin(n*M_PI*i*dx)*sin(m*M_PI*j*dy)*exp(-pow(M_PI,2)*(n*n + m*m)*t*dt);
                    }
                }
                #pragma omp critical
                uSs = solver.uSteadyState(i*dx, j*dy);
                analyticalMatrixU2D(i,j) = uSs -4.*tempSum;
            }
        }
        analyticalMatrixU2D.save(outfileNameAnalytical+to_string(counter)+ ".txt" , raw_ascii);
        counter +=1;
    }
}

double gaussQuad(int numberOfSummationPoints, int n, int m, int numberOfIntegrationPoints){
    int k=0;
    int integrationCounter=0;
    double leftIntegrationLimit = 0.;
    double rightIntegrationLimit = 1.;
    double *xPoints = new double [numberOfIntegrationPoints];
    double *w = new double [numberOfIntegrationPoints];
    gauleg(leftIntegrationLimit , rightIntegrationLimit,xPoints, w, numberOfIntegrationPoints);
    double integralSum = 0.;
    double beta;
    #pragma omp parallel for default(shared) private(k, integrationCounter) reduction(+:integralSum)
    for ( k = 1;  k < numberOfSummationPoints; k++){
        beta = 2*k-1;
        double intGaussX = 0.;
        double intGaussY = 0.;
        for (integrationCounter =0; integrationCounter < numberOfIntegrationPoints; integrationCounter++){
            intGaussX += w[integrationCounter]*sin(beta*M_PI*xPoints[integrationCounter])*sin(n*M_PI*xPoints[integrationCounter]);
            intGaussY += w[integrationCounter]*sinh(beta*M_PI*xPoints[integrationCounter])*sin(m*M_PI*xPoints[integrationCounter]);
        }
        integralSum += 4./M_PI*intGaussX*intGaussY/(beta*sinh(beta*M_PI));
    }
    return integralSum;
}
