#include "twodimensionaldiffusionsolver.h"

TwoDimensionalDiffusionSolver::TwoDimensionalDiffusionSolver(double dt_, double dx_, double dy_, double theta_, double T_, int Nx_, int Ny_, int Nt_)
{
    dt = dt_;
    dx = dx_;
    dy = dy_;
    theta = theta_;
    T = T_;
    Nt = Nt_;
    Nx = Nx_;
    Ny = Ny_;
}

void TwoDimensionalDiffusionSolver::solve(string outfileName_)
{
    outfileName = outfileName_;
    mat solutionMatrixU = zeros<mat>(Nx,Ny);

    alpha = dt/pow(dx,2);



    double ** u, ** uOld;

    u = new double*[Nx];
    uOld = new double*[Nx];

    for(int i = 0; i <Nx; i++){
        u[i] = new double[Ny];
        uOld[i] = new double[Ny];
    }

    // Initial conditions and declaration arrays
    for (int i =0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            uOld[i][j] = -uSteadyState(i*dx,j*dy);
            u[i][j] = 0.;
        }
    }

    int counter = 1;
    for (int t = 1; t < Nt; t++){
        //generate_right_hand_side( computed_right_hand_side, uOld, offDiagonalRhs, diagonalRhs, Nx);
        //gassianTridiagonalSymmetricSolver( computed_right_hand_side,  u, uOld, offDiagonalLhs, diagonalLhs, Nx);
        // Change of variables
        explicitScheme(u, uOld, Nx, Ny);
        for(int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                u[i][j] += uSteadyState(i*dx,j*dy);
                cout << u[i][j] << ",";
                solutionMatrixU(i,j) = u[i][j];
            }
            cout << endl;
        }
        solutionMatrixU.save(outfileName + "SolutionMatrixUTime" + to_string(counter) + ".txt", raw_ascii);
        counter += 1;
    }
    delete [] u;
    delete [] uOld;
}

double TwoDimensionalDiffusionSolver::uSteadyState(double x, double y)
{
    double uSs = 0.0;
    double fourOverPi = 4.0/M_PI;
    for (int k = 1; k < 100; k++){
        double beta = 2*k-1;
        uSs += sin(beta*M_PI*x)*sinh(beta*M_PI*y)/(beta*sinh(beta*M_PI));
    }
    uSs *= fourOverPi;
    //cout << "u steady in method " << uSs << endl;
    return uSs;
}

void TwoDimensionalDiffusionSolver::explicitScheme(double **u, double **uOld, int Nx, int Ny)
{
    for (int i = 0 ; i < Nx; i++){
        u[i][0] = u[0][i] = u[i][Nx-1] = u[Nx-1][i] = 0.0;
    }

    for (int i = 1; i < Nx-1; i++){
        for (int j = 1; j < Ny-1; j++){
            u[i][j] = uOld[i][j] + alpha*(uOld[i+1][j] + uOld[i-1][j] + uOld[i][j+1] + uOld[i][j-1] - 4*uOld[i][j]);
        }
    }

    //delete [] uOld;
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            uOld[i] = u[i];
        }
    }
}
