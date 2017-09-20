// 1st try on Jacobi algorithm for finding eigenvalues

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

using namespace std;
using namespace arma;

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& rhoMax,int& maxIterations, double& tolerance, string& armadillo, string& interactionRepulsion, double& omega, int argc, char** argv );
void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter, mat& v);
mat get_eigenvecs(mat a, mat v, colvec eigenValues, int N);
double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N);
void rotate(mat &aMatrix, int k, int l, int N, mat& v);
void createEigenvalueVector( mat A, colvec &eigenValues, int N );
void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h, string& interactionRepulsion, double omega);
void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax, double omega, colvec eigenValues);
void output_vectors( double *, int, int, string);
void calculateError(colvec eigenValues, double *computedError);

ofstream ofile1; // File for scalars
ofstream ofile2; // File for vectors

main(int argc, char* argv[]){ 

  if (3 == 3)
      int result = Catch::Session().run( argc, argv );

  double tolerance, computedError, h, timeUsed, omega, rhoMax;
  int N, amplificationFactor, numberOfSimulations, maxIterations, counter;
  string outfileName, outfileNameComputedNumerical, outfileNameScalars, armadillo, interactionRepulsion;
  mat A, v, eigenvectorMatrixSorted3;
  colvec eigenValues;

  clock_t start, finish;

  initialize(outfileName, numberOfSimulations, amplificationFactor,N, rhoMax, maxIterations, tolerance, armadillo, interactionRepulsion, omega, argc, argv );

  if (armadillo == "false"){
      outfileNameScalars = (outfileName) + string("_scalars")+string(".csv");
      outfileNameComputedNumerical = (outfileName) + string("_numerical");
      ofile1.open(outfileNameScalars);
  }
  else{
      outfileNameScalars = (outfileName) + string("Armadillo") + string("_scalars")+string(".csv");
      outfileNameComputedNumerical = (outfileName) + string("_numerical");
      ofile1.open(outfileNameScalars);  
  }
  ofile1 << "rhoMax,omega,h,logH,relError,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter,lambda1,lambda2,lambda3" << endl;

  // Solving for different matrix dimensions
  for (int simulationNumber = 0; simulationNumber < numberOfSimulations; simulationNumber++){
      createTridiagonalMatrix(A, N, rhoMax, 0, &h, interactionRepulsion, omega);
      v = eye<mat>(N,N);
      eigenValues  = zeros<colvec>(N);

      start = clock();
      if (armadillo == "false")
          jacobi(A, eigenValues, tolerance, maxIterations, N, &counter, v);
      else
          eig_sym(eigenValues, A);
      finish = clock();
      timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));

      calculateError(eigenValues, &computedError);
      if(armadillo == "false")
          output_scalars(computedError, h, timeUsed, N, counter, rhoMax, omega, eigenValues);
      else{
          output_scalars(computedError, h, timeUsed, N, 0, rhoMax, omega, eigenValues);
      }

      if (simulationNumber < numberOfSimulations -1)
          N *= amplificationFactor;
      eigenValues.print("EigenValues: ");
  }

  outfileName = (outfileName) + string("_eigenVector.csv");
  eigenvectorMatrixSorted3 = get_eigenvecs(A, v, eigenValues, N);
  eigenvectorMatrixSorted3.save(outfileName, csv_ascii);

  ofile1.close();
  ofile2.close();

  return 0;
}

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& rhoMax,int& maxIterations, double& tolerance, string& armadillo, string& interactionRepulsion, double& omega, int argc, char** argv )
{
    if( argc<= 1){
      cout << "Insert: outfile-name, number of simulations, amplification factor, start dimension" << endl;
      exit(1);
    }
    else{
      outfile_name=argv[1];
    }
    number_of_simulations = atoi(argv[2]);
    amplificationFactor = atoi(argv[3]);
    N = atoi(argv[4]);
    rhoMax = atof(argv[5]);
    maxIterations = atoi(argv[6]);
    tolerance = atof(argv[7]);
    armadillo = argv[8];
    interactionRepulsion = argv[9];
    omega = atof(argv[10]);
}

//ofile1 << "logH,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter" << endl;
void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax, double omega, colvec eigenValues){
  ofile1 << setiosflags(ios::showpoint | ios::uppercase);
  ofile1 << setw(15) << setprecision(16) << rhoMax << ", ";
  ofile1 << setw(15) << setprecision(16) << omega << ", ";
  ofile1 << setw(15) << setprecision(16) << h << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(h) << ", ";
  ofile1 << setw(15) << setprecision(16) << computedError << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(computedError) << ", ";
  ofile1 << setw(15) << setprecision(16) << timeUsed << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(timeUsed) << ", ";
  ofile1 << setw(15) << setprecision(16) << N << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(N) << ", ";
  ofile1 << setw(15) << setprecision(16) << counter << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(counter) << ", ";
  ofile1 << setw(15) << setprecision(16) << eigenValues(0)<< ", ";
  ofile1 << setw(15) << setprecision(16) << eigenValues(1) << ", ";
  ofile1 << setw(15) << setprecision(16) << eigenValues(2) << endl;
}
//"rhoMax,omega, h,logH,relError,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter,lambda1,lambda2,lambda3"

void output_vectors( double *output_array, int simulation_number, int N, string outfile_name){
  string filename = outfile_name + to_string(simulation_number+1)+string(".csv");
  ofile2.open(filename);
  /*
  ofile2 << "simulation_number_" << simulation_number+1 << endl;
  ofile2 << setiosflags(ios::showpoint | ios::uppercase);
  for (int i = 0; i<N; i++){
    ofile2 << setw(15) << setprecision(16) << output_array[i] << endl;;
  }
  ofile2 << setw(15) << setprecision(16) << endl;
  */
  ofile2.close();
}

void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter, mat& v){
    int k,l;
    *counter = 0;

    double aMaxNonDiagonal = findMaxNonDiagonalElement(A, &k, &l, N);

    while ( fabs(aMaxNonDiagonal) > tolerance && *counter < maxIterations ){
        aMaxNonDiagonal  = findMaxNonDiagonalElement(A, &k, &l, N);
        rotate(A, k, l, N, v);
        *counter += 1;
    }
    createEigenvalueVector( A, eigenValues,  N);
}

mat get_eigenvecs(mat a, mat v, colvec eigenValues, int N){
    // Based on Morten's example
    mat vecs(N,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<N;j++){
            if(a(j,j)==eigenValues(i)){
                for(int k=0;k<N;k++){
                      vecs(k,i)=v(k,j);
                }
             }
         }
    }
    return vecs;
}

void rotate(mat &aMatrix, int k, int l, int N, mat& v){
    double c, s;
    if ( aMatrix(k,l) != 0.0 ) {
        double tau, t;
        tau = (aMatrix(l,l) - aMatrix(k,k))/(2.*aMatrix(k,l));
        if ( tau >= 0 ) {
          t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
          t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }

        c = 1/sqrt(1+t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    double aIk, aIl, aKk, aLl, vik, vil;
    for (int i = 0; i < N; i++){
        if ( i != k && i != l ) {
                aIk = aMatrix(i, k);
                aIl = aMatrix(i, l);
                aMatrix(i,k) = aIk*c - aIl*s;
                aMatrix(k, i) = aMatrix(i, k);
                aMatrix(i, l) = aIl*c + aIk*s;
                aMatrix(l, i) = aMatrix(i, l);
        }
        vik=v(i,k);
        vil=v(i,l);
        v(i,k)=c*vik-s*vil;
        v(i,l)=c*vil+s*vik;
    }
    aKk = aMatrix(k,k);
    aLl = aMatrix(l, l);
    aMatrix(k, k) = aKk*c*c - 2.*aMatrix(k, l)*c*s + aLl*s*s ;
    aMatrix(l, l) = aLl*c*c + 2.*aMatrix(k, l)*c*s + aKk*s*s;
    aMatrix(k,l) = 0.0;
    aMatrix(l,k) = 0.0;
}

double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N){
    double aMaxTemp = 0.0;
    for (int row = 0; row < N; row++){
         for (int col = row + 1; col < N; col++){
             if( row != col && fabs(A(row, col)) >  aMaxTemp ) {
                 *k = row;
                 *l = col;
                 aMaxTemp = fabs(A(row, col));
             }
         }
    }
    return aMaxTemp;
}

void createEigenvalueVector( mat A, colvec &eigenValues, int N ){
    for (int row = 0; row < N; row++){
        for (int col= 0; col < N; col++){
            if ( col == row )
                eigenValues(col) = A(row, col);
        }
    }
    eigenValues = sort(eigenValues);
}

void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h, string& interactionRepulsion, double omega){
    double hTemp = (rhoMax - rhoMin)/N;
    //A.zeros(N,N);
    A = zeros<mat>(N,N);
    double offDiagonal = -1.0/(hTemp*hTemp);
    double diagonalFirstTerm = 2.0/(hTemp*hTemp);

    A(0, 1) =  offDiagonal;
    A(N-1, N-2) = offDiagonal;

    if ( interactionRepulsion == "TwoElectronNoCoulomb" )
        A(0, 0) = diagonalFirstTerm + pow(omega, 2)*hTemp*hTemp;
    else if (interactionRepulsion == "TwoElectronCoulomb" )
        A(0, 0) = diagonalFirstTerm + pow(omega, 2)*hTemp*hTemp + 1./hTemp;
    else
        A(0, 0) = diagonalFirstTerm + hTemp*hTemp;

    for (int row = 1; row < N-1; row++){
        if ( interactionRepulsion == "TwoElectronNoCoulomb" )
            A(row, row) = diagonalFirstTerm + pow(omega, 2)*pow((row+1)*hTemp,2);
        else if (interactionRepulsion == "TwoElectronCoulomb" )
            A(row, row) = diagonalFirstTerm + pow(omega, 2)*pow((row+1)*hTemp,2) + 1./((row+1)*hTemp); // Check this! AM avoiding rho=0 on the first step.
        else
            A(row, row) = diagonalFirstTerm + pow((row+1)*hTemp,2);
      A(row, row - 1) = offDiagonal;
      A(row, row + 1) =  offDiagonal;
      }

    if ( interactionRepulsion == "TwoElectronNoCoulomb" )
         A(N-1, N-1) = diagonalFirstTerm + pow(omega, 2)*pow(rhoMax,2);
    else if (interactionRepulsion == "TwoElectronCoulomb" )
        A(N-1, N-1) = diagonalFirstTerm + pow(omega, 2)*pow(rhoMax,2) + 1./rhoMax;
    else
        A(N-1, N-1) = diagonalFirstTerm + pow(rhoMax,2);
    *h = hTemp;
    A.print("A: ");
}

void calculateError(colvec eigenValues, double *computedError){
    // Calculates sup-norm for relative error of 3 lowest eigenvalues
    colvec analyticalSolution  = zeros<colvec>(3);
    analyticalSolution(0) = 3.0;
    analyticalSolution(1) = 7.0;
    analyticalSolution(2) = 11.0;
    *computedError = fabs((analyticalSolution(0) - eigenValues(0))/analyticalSolution(0));;
    double temp;
    for (int i = 1; i < 3; i++) {
        temp = fabs((analyticalSolution(i) - eigenValues(i))/analyticalSolution(i));
        if (temp > *computedError)
            *computedError = temp;
    }
}

TEST_CASE( "Testing maxDiagonal", "[findMaxNonDiagonalElement]" ) {
  // Checks that function finds max abs non-diagonal element in known symmetric matrix
  int k, l, N;
  N = 3;

  mat A(N,N);
  A(0,0) = 10.;
  A(0,1) = 1.;
  A(0,2) = -6.;
  A(1,0) = 1.;
  A(1,1) = 20.;
  A(1,2) = 3.;
  A(2,0) = -6.;
  A(2,1) = 3.;
  A(2,2) = 30.;

  double max = findMaxNonDiagonalElement(A, &k, &l, N);
  REQUIRE( max == 6.);
}


TEST_CASE( "Testing Jacobi eigenvalues", "[jacobi]" ) {
  // Checks that function finds max abs non-diagonal element in known symmetric matrix
  int k, l, N;
  N = 3;
  colvec eigenValues  = zeros<colvec>(N);
  vec EigvalArmadillo(N);
  double tolerance = pow(10,-10); // Tolerance for non-diagonals in Jacobi
  int maxIterations = pow(10,4);
  double tol = pow(10,-1); // Tolerance for difference jacobi and Armadillo
  mat v = eye<mat>(N,N);

  mat A(N,N);
  A(0,0) = 10.;
  A(0,1) = 1.;
  A(0,2) = -6.;
  A(1,0) = 1.;
  A(1,1) = 20.;
  A(1,2) = 3.;
  A(2,0) = -6.;
  A(2,1) = 3.;
  A(2,2) = 30.;
  int counter;
  eig_sym(EigvalArmadillo, A);
  jacobi(A, eigenValues, tolerance, maxIterations, N, &counter, v);
  for (int i = 0; i < N; i++){
          REQUIRE( fabs((eigenValues(i) - EigvalArmadillo(i))/EigvalArmadillo(i)) < tol);
  }
}


/*
TEST_CASE( "Testing jacobi", "[jacobi]" ) {
  mat A;
  int N = pow(10,3)
  createTridiagonalMatrix(A, N, 10, 0);
  colvec eigenValues;
  jacobi(A, eigenValues, pow(10,-7), pow(10,3), N );
  vec Eigval(N);
  eig_sym(Eigval, A);
  REQUIRE( A(0,0)== 2);
}
*/
