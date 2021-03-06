#include "catch.hpp"
#include "jacobi.h"
#include "eigenvalueBisection.h"
#include "lanczos.h"

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

TEST_CASE( "Testing sturm sequence", "[sturmSeq]" ){
    // Testing example on p 360-361 in Kiusalaas (2014).
    int N = 4;
    double tol = 0.00001;
    mat A(N,N);
    colvec p;
    double lam = 0.5;


    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;


    p = sturmSeq(A, lam, N);
    REQUIRE( abs(p(0) - 1) < tol);
    REQUIRE( abs(p(1) - 1.5) < tol);
    REQUIRE( abs(p(2) - 1.25) < tol);
    REQUIRE( abs(p(3) - 0.375) < tol);
    REQUIRE( abs(p(4) + 0.6875) < tol);
}


TEST_CASE( "Testing numner of lambdas", "[numLambdas]" ){
    // Testing example on p 361 in Kiusalaas (2014).
    // The smallest eigenvalue is supposed to be between 0.25 and 2.25


    int N = 4;
    double tol = 0.00001;
    mat A(N,N);
    colvec p;

    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    double lam = 0.2;
    p = sturmSeq(A, lam, N);
    REQUIRE( numLambdas(p, N) < tol);

    lam = 0.3;
    p = sturmSeq(A, lam, N);
    REQUIRE( numLambdas(p, N) - 1 < tol);

    lam = 10;
    p = sturmSeq(A, lam, N);
    REQUIRE( numLambdas(p, N) >= 1);
}


TEST_CASE( "Testing gerschgorin ", "[gerschgorin]" ){
    // Testing example on p 361 in Kiusalaas (2014).
    // The smallest eigenvalue is supposed to be between 0.25 and 2.25


    int N = 4;
    double tol = 0.00001;
    mat A(N,N);
    colvec eigenvaluesArmadillo;
    eigenvaluesArmadillo = zeros<colvec>(N);
    colvec gerschgorinValues;
    gerschgorinValues = zeros<colvec>(2);

    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    eig_sym(eigenvaluesArmadillo, A);
    gerschgorinValues = gerschgorin(A, N);
    REQUIRE(gerschgorinValues(0) - eigenvaluesArmadillo(0) <= tol);
    REQUIRE(gerschgorinValues(1) - eigenvaluesArmadillo(3) >= tol);
}


TEST_CASE( "Testing lambda range ", "[lamRange]" ){
    // Testing example on p 361 in Kiusalaas (2014).
    // The smallest eigenvalue is supposed to be between 0.25 and 2.25
    // Testing both revised and original version


    int N = 4;
    double tol = 0.00001;
    mat A(N,N);
    colvec eigenvaluesArmadillo;
    eigenvaluesArmadillo = zeros<colvec>(N);
    colvec eigenvalueDomains;
    eigenvalueDomains = zeros<colvec>(N+1);

    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    eig_sym(eigenvaluesArmadillo, A);
    for (int i = 0; i < 2; i++){
        eigenvalueDomains = lamRange(A, N, i);
        for (int eigenvalueNumber = 0; eigenvalueNumber < N; eigenvalueNumber++){
            REQUIRE(eigenvaluesArmadillo(eigenvalueNumber) >= eigenvalueDomains(eigenvalueNumber));
            REQUIRE(eigenvaluesArmadillo(eigenvalueNumber) <= eigenvalueDomains(eigenvalueNumber+1));
        }
    }
}

TEST_CASE( "Testing f ", "[f]" ){
    // Check that f gives last element of Sturm sequence
    int N = 4;
    double tol = 0.00001;
    mat A(N,N);
    colvec p;
    double lam = 0.5;
    double fResult;


    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    p = sturmSeq(A, lam, N);
    fResult = f(A, lam, N);

    REQUIRE(fabs(fResult - p(4)) < tol);
}

TEST_CASE( "Testing bisection ", "bisection" ){
    // Bisection method double bisection(double (*func)(double), double x1, double x2, double xacc, mat &A, int N, int max_iterations)
    double x1 = 0.25; // lower limit minimum eigenvalue
    double x2 = 0.5; //
    double xacc = 0.000001;
    int max_iterations = 1000000;
    double bisectionSolution;
    double tol = 0.00001;

    int N = 4;
    colvec  eigenvaluesArmadillo;
    eigenvaluesArmadillo = zeros<colvec>(N);
    mat A(N,N);

    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    eig_sym(eigenvaluesArmadillo, A);
    bisectionSolution = bisection(f, x1, x2, xacc, A, N, max_iterations);

    REQUIRE(fabs(bisectionSolution- eigenvaluesArmadillo(0)) < tol);
}

TEST_CASE( "Testing eigenvalues3", "eigenvals3" ){
    // eigenvals3(mat &A, int N, double bisectionAccuracy, int max_iterations)
    double xacc = 0.000001;
    int max_iterations = 1000000;
    double tol = 0.00001;

    int N = 4;
    colvec  eigenvaluesArmadillo;
    eigenvaluesArmadillo = zeros<colvec>(N);
    colvec  eigenvaluesBisection;
    eigenvaluesBisection= zeros<colvec>(N);

    mat A(N,N);

    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    eig_sym(eigenvaluesArmadillo, A);

    for ( int revision = 0; revision < 2; revision++){
        eigenvaluesBisection= eigenvals3(A, N, xacc, max_iterations, revision);
        eigenvaluesBisection = sort(eigenvaluesBisection);
        for (int i = 0; i < N; i++){
            REQUIRE(fabs(eigenvaluesArmadillo(i)- eigenvaluesBisection(i)) < tol);
        }
    }
}

TEST_CASE( "NumLambdas revised", "[numLambdasRevised]" ){
    // Using same test case as for numLambdas


    int N = 4;
    double tol = 0.00001;
    mat A(N,N);
    colvec p;

    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    double lam = 0.2;
    p = sturmSeqRevised(A, lam, N);
    REQUIRE( numLambdasRevised(p, N) < tol);

    lam = 0.3;
    p = sturmSeqRevised(A, lam, N);
    REQUIRE( numLambdasRevised(p, N) - 1 < tol);

    lam = 10;
    p = sturmSeqRevised(A, lam, N);
    REQUIRE( numLambdasRevised(p, N) >= 1);
}

TEST_CASE( "lanczos brute force. Armadillo on lanczos tridiagonal matrix", "[lanczos]" ){
    /* Assume symmetric (non-tridiagonal) matrix.
    Use brute force to compute Aq
    Show that Q is orthogonal
    Show that T is trdiagonal
    Test largest eigenvalue against Armadillo
    Test that stopping criterion also gives correct eigenvalue


    */

    int N = 200;
    mat A, Q, T, QtransQ;
    A.randu(N,N);
    A = A*A.t();
    int iterations = 5;
    int stopIteration;

    colvec eigenvaluesArmadillo, eigenvaluesLanczos,alpha, beta, eigenvaluesStopCriterion;
    eigenvaluesArmadillo.zeros(N);
    //eigenvaluesLanczos.zeros(iterations);

    string tridiag = "false";
    string eigenvalueSolver = "armadillo";

    lanczos(eigenvaluesLanczos, A, alpha, beta, Q, N, iterations, tridiag, eigenvalueSolver, &stopIteration);

    QtransQ = trans(Q)*Q;
    T = trans(Q)*A*Q;
    cout << "Non-tridiagonal, brute force. Prespecified iteration limit 5 " << endl;
    QtransQ.print("Prespecified iteration limit 5 trans(Q)*Q: ");
    T.print("Prespecified iteration limit 5 T = trans(Q)*A*Q: ");
    eig_sym(eigenvaluesArmadillo, A);
    eigenvaluesLanczos.print("Prespecified iteration limit 5 EigLanczos brute force: ");

    double tol = 1.e-6;
    REQUIRE( abs(eigenvaluesLanczos(iterations-1)/eigenvaluesArmadillo(N-1) -1) < tol);


    // Test stop iteration requirement
    cout << " Applying stopping critera" << endl;
    iterations = N;
    lanczos(eigenvaluesStopCriterion, A, alpha, beta, Q, N, iterations, tridiag, eigenvalueSolver, &stopIteration);
    //eigenvaluesLanczos.print("Eigenvalues stop set at 5 ");
    eigenvaluesStopCriterion.print("Eigenvalues with stop criterion activesed");
    cout << "Stopping criteria active. Number of iterations: " << stopIteration << endl;
    int lastElementIndex = eigenvaluesStopCriterion.n_elem;
    REQUIRE( abs(eigenvaluesStopCriterion(lastElementIndex-1)-eigenvaluesArmadillo(N-1)) < tol);
    REQUIRE( lastElementIndex  == stopIteration);

}



TEST_CASE( "lanczos tridiagonal not brute force Aq. Armadillo on lanczos tridiagonal matrix", "[lanczos]" ){
    /* Assume symmetric tridiagonal matrix.
    Calculate Aq more efficient than brute force*/

    string tridiag = "true";
    string eigenvalueSolver;
    int N = 5;
    int iterations = 5;
    int stopIteration;
    colvec eigenvaluesArmadillo, eigenvaluesLanczos,alpha, beta;
    mat QtransQ, Q, T;

    mat A(N,N);
    eigenvaluesArmadillo.zeros(N);
    eigenvaluesLanczos.zeros(iterations);

    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;
    A(0,4) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;
    A(1,4) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;
    A(2,4) = 0.;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;
    A(3,4) = -1.;

    A(4,0) = 0.;
    A(4,1) = 0.;
    A(4,2) = 0.;
    A(4,3) = -1.;
    A(4,4) = 2.0;


    double tol = 1.e-2; // 1 per cent
    eig_sym(eigenvaluesArmadillo, A);

    eigenvalueSolver = "armadillo";
    lanczos(eigenvaluesLanczos, A, alpha, beta, Q, N, iterations, tridiag, eigenvalueSolver, &stopIteration);
    QtransQ = trans(Q)*Q;
    T = trans(Q)*A*Q;
    cout << "Tridiagonal not brute force  " << endl;
    QtransQ.print("trans(Q)*Q: ");
    eigenvaluesLanczos.print("EigLanczos not brute force: ");
    REQUIRE( abs(eigenvaluesLanczos(iterations-1)/eigenvaluesArmadillo(N-1)-1) < tol);

}


TEST_CASE( "lanczos tridiagonal not brute force Aq. Jacobi on lanczos tridiagonal matrix", "[lanczos]" ){

    /* Assume symmetric tridiagonal matrix.
    Calculate Aq more efficient than brute force*/

     /*
    string tridiag = "true";
    string eigenvalueSolver;
    int N = 4;
    int iterations = 3;
    int stopIteration;
    colvec eigenvaluesArmadillo, eigenvaluesLanczos,alpha, beta;
    mat QtransQ, Q, T;

    mat A(N,N);
    eigenvaluesArmadillo.zeros(N);
    eigenvaluesLanczos.zeros(iterations);


    A(0,0) = 2.;
    A(0,1) = -1.;
    A(0,2) = 0.;
    A(0,3) = 0.;

    A(1,0) = -1.;
    A(1,1) = 2.;
    A(1,2) = -1.;
    A(1,3) = 0.;

    A(2,0) = 0.;
    A(2,1) = -1.;
    A(2,2) = 2.;
    A(2,3) = -1.0;

    A(3,0) = 0.;
    A(3,1) = 0.;
    A(3,2) = -1.;
    A(3,3) = 2.0;

    double tol = 1.e-2; // 1 per cent
    eig_sym(eigenvaluesArmadillo, A);

    eigenvalueSolver = "jacobi";
    lanczos(eigenvaluesLanczos, A, alpha, beta, Q, N, iterations, tridiag, eigenvalueSolver, stopIteration);
    QtransQ = trans(Q)*Q;
    T = trans(Q)*A*Q;
    cout << "Tridiagonal not brute force  " << endl;
    QtransQ.print("trans(Q)*Q: ");
    REQUIRE( abs(eigenvaluesLanczos(iterations-1)/eigenvaluesArmadillo(N-1)-1) < tol);
    */
}

TEST_CASE( "Armadillo's head and size method" ){
    // Getting to know Armadillo
    colvec vector1, vector2, vectorDiff, vectorDiffExact;
    int N = 6;
    double tolerance = 1.e-7;
    vector1.zeros(N);
    vector2.zeros(N);
    vectorDiffExact.zeros(3);
    for (int i=0; i <N; i++){
        vector1(i) = double(i);
        vector2(i) = double(2*i);
        if (i < 3)
            vectorDiffExact(i) = vector2(i) - vector1(i);
    }

    vectorDiff = vector2.head(3) - vector1.head(3);
    cout << "Armadillo testing " << endl;
    vectorDiff.print("Vector diff with head: ");
    vectorDiffExact.print("Vector diff manual: ");

    for (int i = 0; i < 3; i++){
        REQUIRE( abs(vectorDiff(i)-vectorDiffExact(i)) < tolerance);
    }
    REQUIRE( size(vectorDiff) == size(vectorDiffExact));

    mat A = randn<mat>(5,5);
    mat B = A.cols(0,3);
    A.print("Armadillo test A = ");
    B.print("Armadillo test B = 1st 4 col's of A ");

}

TEST_CASE( "Gram-Schmidt", "[lanczos]" ){
    /* Test Gram-Schmidth
    void gramSmith(mat &QforEigenvalue, int N, int columnNumber)
    Test matrix taken from wikipedia     https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
*/
    int N = 2;
    int columnNumber = 2;
    mat A = randu<mat>(2,2);
     double tolerance = 1.e-12;

    A(0,0) = 3.;
    A(0,1) = 2.;

    A(1,0) = 1.;
    A(1,1) = 2.;

    mat B;
    B.zeros(N,N);
    B(0,0) = 3./sqrt(10);
    B(0,1) = -1./sqrt(10);
    B(1,0) = 1./sqrt(10);
    B(1,1) = 3./sqrt(10);

    gramSmith(A, N, columnNumber);

    A.print("A Gram-Schmidt: ");
    B.print("B (Exact)) Gram-Schmidt: ");

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            REQUIRE( abs(A(i,j)/B(i,j) - 1.) < tolerance);
        }
    }
}
