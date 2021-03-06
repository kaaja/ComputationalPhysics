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

mat TwoDimensionalDiffusionSolver::solve(string outfileName_, string method_, int thread_num_)
{
    outfileName = outfileName_;

    int thread_num;
    omp_set_num_threads(thread_num_);
    thread_num = omp_get_max_threads ();
    colvec iterationNumberArray= zeros<colvec>(1);
    double iterationNumber = 0.0;
    mat solutionMatrixU = zeros<mat>(Nx,Ny);
    int maxIterations = 100000;
    double maxDifference = 0.000001;

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
        if (method_ =="Explicit") explicitScheme(u, uOld, Nx, Ny);
        else if (method_ == "ImplicitJacobi") backwardEulerJacobi(u, uOld, Nx, Ny, maxIterations, maxDifference, iterationNumber);
        else if (method_ == "ImplicitGaussSeidel") backwardEulerGaussSeidel(u, uOld, Nx, Ny, maxIterations, maxDifference, iterationNumber);
        for(int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                u[i][j] += uSteadyState(i*dx,j*dy);
                solutionMatrixU(i,j) = u[i][j];
            }
        }
        solutionMatrixU.save(outfileName + "SolutionMatrixUTime" + to_string(counter) + ".txt", raw_ascii);
        iterationNumberArray(0) = iterationNumber;
        counter += 1;
    }
    iterationNumberArray.save(outfileName + "IterationNumber.txt", raw_ascii);
    DestroyMatrix(u, Nx, Ny);
    DestroyMatrix(uOld, Nx, Ny);
    return solutionMatrixU;
}

double TwoDimensionalDiffusionSolver::uSteadyState(double x, double y)
{
    int k;
    double uSs = 0.0;
    double fourOverPi = 4.0/M_PI;
    #pragma omp parallel for default(shared) private(k) reduction(+:uSs)
    for (k = 1; k < 15; k++){
        double beta = 2*k-1;
        uSs += sin(beta*M_PI*x)*sinh(beta*M_PI*y)/(beta*sinh(beta*M_PI));
    }
    uSs *= fourOverPi;
    return uSs;
}

void TwoDimensionalDiffusionSolver::explicitScheme(double **u, double **uOld, int Nx, int Ny)
{
    int i; int j;
    for (i = 0 ; i < Nx; i++){
        u[i][0] = u[0][i] = u[i][Nx-1] = u[Nx-1][i] = 0.0;
    }
    #pragma omp parallel for default(shared) private(i,j)
    for (i = 1; i < Nx-1; i++){
        for (j = 1; j < Ny-1; j++){
            u[i][j] = uOld[i][j] + alpha*(uOld[i+1][j] + uOld[i-1][j] + uOld[i][j+1] + uOld[i][j-1] - 4*uOld[i][j]);
        }
    }

    //set uOld = u;
    setMatrixAEqualMatrixB(uOld, u, Nx, Ny);
}

// Backward Euler, 2D, Jacobi
void TwoDimensionalDiffusionSolver::backwardEulerJacobi(double **u, double **uOld, int Nx, int Ny,
                                                  int maxIterations, double maxDifference, double& iterationNumber){
    int i; int j;
    for (i = 0 ; i < Nx; i++){
        u[i][0] = u[0][i] = u[i][Nx-1] = u[Nx-1][i] = 0.0;
    }
    double **uTemp = CreateMatrix(Nx, Ny);
    double diff = 1.;
    int iterations = 1;
    while ((iterations < maxIterations) && (diff > maxDifference)){
        diff = 0.;
        setMatrixAEqualMatrixB(uTemp, u, Nx, Ny);
        #pragma omp parallel for default(shared) private(i,j) reduction(+:diff)
        for (i = 1; i < Nx-1; i++){
            for (j = 1; j < Ny-1; j++){
                u[i][j] =   1./(1+4.*alpha)*(alpha*(uTemp[i+1][j] + uTemp[i-1][j] + uTemp[i][j+1] + uTemp[i][j-1]) + uOld[i][j]);
                diff += fabs(u[i][j] - uTemp[i][j]);
            }
        }
        diff /= (Nx*Ny);
        iterations += 1;
    }
    iterationNumber = iterations/double(Nx);
    if ((iterations == maxIterations) && (diff > maxDifference))
        cout << "Iteration limit reached without convergence. |u^(k+1) - u^k| = " << diff << endl;
    setMatrixAEqualMatrixB(uOld, u, Nx, Ny);
    DestroyMatrix(uTemp, Nx, Ny);
}

void TwoDimensionalDiffusionSolver::backwardEulerGaussSeidel(double **u, double **uOld, int Nx, int Ny,
                                                  int maxIterations, double maxDifference, double& iterationNumber){
    int i; int j;
    for (i = 0 ; i < Nx; i++){
        u[i][0] = u[0][i] = u[i][Nx-1] = u[Nx-1][i] = 0.0;
    }
    double **uTemp = CreateMatrix(Nx, Ny);
    double diff = 1.;
    int iterations = 1;
    while ((iterations < maxIterations) && (diff > maxDifference)){
        diff = 0.;
        setMatrixAEqualMatrixB(uTemp, u, Nx, Ny);
        for (i = 1; i < Nx-1; i++){
            for (j = 1; j < Ny-1; j++){
                u[i][j] =   1./(1+4.*alpha)*(alpha*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) + uOld[i][j]);
                diff += fabs(u[i][j] - uTemp[i][j]);
            }
        }
        diff /= (Nx*Ny);
        iterations += 1;
    }
    iterationNumber = iterations/double(Nx);
    if ((iterations == maxIterations) && (diff > maxDifference))
        cout << "Iteration limit reached without convergence. |u^(k+1) - u^k| = " << diff << endl;
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

void TwoDimensionalDiffusionSolver::calculate_error(mat computed_numerical_solution, mat computed_exact_solution, double *computed_error, int Nx_, int Nt_)
{
    // L2 norm
    double temp_relative_error = 0.0;
    *computed_error = 0.0;
    double sum = 0.;
    for (int i = 0; i < Nx_; i++){
        for (int j = 0; j < Nx_; j++){
            sum += pow((computed_numerical_solution(i, j) - computed_exact_solution(i, j)),2);
        }
    }
    sum = sqrt(dx*dx*sum);
    *computed_error = sum;
}
