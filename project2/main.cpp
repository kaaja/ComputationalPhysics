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

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, int& rhoMax,int& maxIterations, double& tolerance, string& armadillo, int argc, char** argv );
void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter);
double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N);
void rotate(mat &aMatrix, int k, int l, int N);
void createEigenvalueVector( mat A, colvec &eigenValues, int N );
void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h);
void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax);
void output_vectors( double *, int, int, string);
void calculateError(colvec eigenValues, double *computedError);

ofstream ofile1; // File for scalars
ofstream ofile2; // File for vectors

main(int argc, char* argv[]){ 

  if (3 == 3)
      int result = Catch::Session().run( argc, argv );

  double tolerance, computedError, h, timeUsed;
  int N, rhoMax, amplificationFactor, numberOfSimulations, maxIterations, counter;
  string outfileName, outfileNameComputedNumerical, outfileNameScalars, armadillo;
  mat A;
  colvec eigenValues;

  clock_t start, finish;

  initialize(outfileName, numberOfSimulations, amplificationFactor,N, rhoMax, maxIterations, tolerance, armadillo, argc, argv );

  if (armadillo == "false"){
      outfileNameScalars = (outfileName) + string("_scalars")+string(".csv");
      outfileNameComputedNumerical = (outfileName) + string("_numerical");
      ofile1.open(outfileNameScalars);
      ofile1 << "rhoMax,h,logH,relError,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter" << endl;
  }
  else{
      outfileNameScalars = (outfileName) + string("Armadillo") + string("_scalars")+string(".csv");
      outfileNameComputedNumerical = (outfileName) + string("_numerical");
      ofile1.open(outfileNameScalars);
      ofile1 << "rhoMax,h,logH,relError,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter" << endl;

  }

  // Solving for different matrix dimensions
  for (int simulationNumber = 0; simulationNumber < numberOfSimulations; simulationNumber++){
      createTridiagonalMatrix(A, N, rhoMax, 0, &h);
      eigenValues  = zeros<colvec>(N);

      start = clock();
      if (armadillo == "false")
          jacobi(A, eigenValues, tolerance, maxIterations, N, &counter);
      else
          eig_sym(eigenValues, A);

      finish = clock();
      timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));

      calculateError(eigenValues, &computedError);
      if(armadillo == "false")
          output_scalars(computedError, h, timeUsed, N, counter, rhoMax);
      else{
          output_scalars(computedError, h, timeUsed, N, 0, rhoMax);
      }

      if (simulationNumber < numberOfSimulations -1)
          N *= amplificationFactor;
      eigenValues.print("EigenValues: ");
  }

  ofile1.close();
  ofile2.close();


  return 0;
}

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, int& rhoMax,int& maxIterations, double& tolerance, string& armadillo, int argc, char** argv )
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
    rhoMax = atoi(argv[5]);
    maxIterations = atoi(argv[6]);
    tolerance = atof(argv[7]);
    armadillo = argv[8];
}

//ofile1 << "logH,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter" << endl;
void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax){
  ofile1 << setiosflags(ios::showpoint | ios::uppercase);
  ofile1 << setw(15) << setprecision(16) << rhoMax << ", ";
  ofile1 << setw(15) << setprecision(16) << h << ", ";
  ofile1 << setw(15) << setprecision(16) << log10(h) << ", ";
  ofile1 << setw(15) << setprecision(16) << computedError << ", ";
  ofile1 << setw(15) << setprecision(16) << log10(computedError) << ", ";
  ofile1 << setw(15) << setprecision(16) << timeUsed << ", ";
  ofile1 << setw(15) << setprecision(16) << log10(timeUsed) << ", ";
  ofile1 << setw(15) << setprecision(16) << N << ", ";
  ofile1 << setw(15) << setprecision(16) << log10(N) << ", ";
  ofile1 << setw(15) << setprecision(16) << counter << ", ";
  ofile1 << setw(15) << setprecision(16) << log10(counter) << endl;
}

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

void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter){
    int k,l;
    *counter = 0;

    double aMaxNonDiagonal = findMaxNonDiagonalElement(A, &k, &l, N);

    while ( fabs(aMaxNonDiagonal) > tolerance && *counter < maxIterations ){
        aMaxNonDiagonal  = findMaxNonDiagonalElement(A, &k, &l, N);
        rotate(A, k, l, N);
        *counter += 1;
    }
    createEigenvalueVector( A, eigenValues,  N);
    cout << " " << N << " counter: " << *counter << " maxIterations: " << maxIterations << "fabs(aMaxNonDiagonal) " << fabs(aMaxNonDiagonal) << "tolerance " << tolerance << endl;
}

void rotate(mat &aMatrix, int k, int l, int N){
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

    double aIk, aIl, aKk, aLl;
    for (int i = 0; i < N; i++){
        if ( i != k && i != l ) {
                aIk = aMatrix(i, k);
                aIl = aMatrix(i, l);
                aMatrix(i,k) = aIk*c - aIl*s;
                aMatrix(k, i) = aMatrix(i, k);
                aMatrix(i, l) = aIl*c + aIk*s;
                aMatrix(l, i) = aMatrix(i, l);
        }
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

void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h){
    double hTemp = (rhoMax - rhoMin)/N;
    //A.zeros(N,N);
    A = zeros<mat>(N,N);
    double offDiagonal = -1.0/(hTemp*hTemp);
    double diagonalFirstTerm = 2.0/(hTemp*hTemp);
    A(0, 1) =  offDiagonal;
    A(0, 0) = 2.0/(hTemp*hTemp) + hTemp*hTemp; // Check this
    A(N-1, N-2) = offDiagonal;

    for (int row = 1; row < N-1; row++){
      A(row, row) = diagonalFirstTerm + pow(row*hTemp,2);
      A(row, row - 1) = offDiagonal;
      A(row, row + 1) =  offDiagonal;
      }
    A(N-1, N-1) = 2.0/(hTemp*hTemp) + pow((N-1),2);
    *h = hTemp;
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
  jacobi(A, eigenValues, tolerance, maxIterations, N, &counter);
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
