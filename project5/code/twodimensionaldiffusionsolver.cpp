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
    int maxIterations = 100000;
    double maxDifference = 0.001;

    alpha = dt/pow(dx,2);

    double ** u, ** uOld;

    uOld = CreateMatrix(Nx, Ny); u = CreateMatrix(Nx, Ny);

    // Initial conditions and declaration arrays
    for (int i =0; i < Nx; i++){
        for (int j = 0; j < Ny; j++)
            uOld[i][j] = -uSteadyState(i*dx,j*dy);
    }

    int counter = 1;
    for (int t = 1; t < Nt; t++){
        //generate_right_hand_side( computed_right_hand_side, uOld, offDiagonalRhs, diagonalRhs, Nx);
        //gassianTridiagonalSymmetricSolver( computed_right_hand_side,  u, uOld, offDiagonalLhs, diagonalLhs, Nx);
        // Change of variables
        //explicitScheme(u, uOld, Nx, Ny);

        backwardEuler(u, uOld, Nx, Ny, maxIterations, maxDifference);
        for(int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                u[i][j] += uSteadyState(i*dx,j*dy);
                solutionMatrixU(i,j) = u[i][j];
            }

        }
        solutionMatrixU.save(outfileName + "SolutionMatrixUTime" + to_string(counter) + ".txt", raw_ascii);
        counter += 1;
    }
    //delete [] u;
    //delete [] uOld;
    DestroyMatrix(u, Nx, Ny);
    DestroyMatrix(uOld, Nx, Ny);
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

    //set uOld = u;
    setMatrixAEqualMatrixB(uOld, u, Nx, Ny);
}

// Backward Euler, 2D, Jacobi
void TwoDimensionalDiffusionSolver::backwardEuler(double **u, double **uOld, int Nx, int Ny,
                                                  int maxIterations, double maxDifference){
    for (int i = 0 ; i < Nx; i++){
        u[i][0] = u[0][i] = u[i][Nx-1] = u[Nx-1][i] = 0.0;
    }
    double **uTemp = CreateMatrix(Nx, Ny);
    double diff = 1.;
    int iterations = 1;
    while ((iterations < maxIterations) && (diff > maxDifference)){
        diff = 0.;
        setMatrixAEqualMatrixB(uTemp, u, Nx, Ny);
        for (int i = 1; i < Nx-1; i++){
            for (int j = 1; j < Ny-1; j++){
                //u[i][j] = uOld[i][j] + alpha*(uTemp[i+1][j] + uTemp[i][j+1] - 4.0*uTemp[i][j] + uTemp[i-1][j] + uTemp[i][j-1]);
                u[i][j] =   1./(1+4.*alpha)*(alpha*(uTemp[i+1][j] + uTemp[i-1][j] + uTemp[i][j+1] + uTemp[i][j-1]) + uOld[i][j]);
                diff += fabs(u[i][j] - uTemp[i][j]);
            }
        }
        diff /= (Nx*Ny);
        iterations += 1;
        cout << diff << "  iterations " << iterations << endl;

    }
    setMatrixAEqualMatrixB(uOld, u, Nx, Ny);
    DestroyMatrix(uTemp, Nx, Ny);
}

double ** TwoDimensionalDiffusionSolver:: CreateMatrix(int m, int n){
  double ** mat;
  mat = new double*[m];
  for(int i=0;i<m;i++){
    mat[i] = new double[n];
    for(int j=0;j<m;j++)
      mat[i][j] = 0.0;
  }
  return mat;
}

void TwoDimensionalDiffusionSolver:: DestroyMatrix(double ** mat, int m, int n){
  for(int i=0;i<m;i++)
    delete[] mat[i];
  delete[] mat;
}

void TwoDimensionalDiffusionSolver::setMatrixAEqualMatrixB(double ** matrixA, double ** matrixB, int m, int n)
{
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            matrixA[i][j] = matrixB[i][j];
        }
    }
}

