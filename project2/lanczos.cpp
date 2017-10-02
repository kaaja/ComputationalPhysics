#include "lanczos.h"
#include "jacobi.h"

void lanczos(colvec &eigenvalues, mat &A, colvec &alpha, colvec &beta, mat &QforEigenvalue, int N, int iterationNumber, string tridiag, string eigenvalueSolver, int *stopIteration){
    int k, lanczosIterationCounter, numberOfEigenvaluesConvergenceTest;
    colvec r, rMinusOne, q, qOld,  Aq, eigenvaluesOld, eigenvaluesDifference;
    double betaMinusOne, normMinimumEigenvalues, LanczosIterationTolerance;
    mat T, Q;
    string eigenvaluesConverged;

    LanczosIterationTolerance = 1.e-3;
    Aq = zeros<colvec>(N);
    beta.zeros(iterationNumber);
    alpha.zeros(iterationNumber+1);
    qOld.zeros(N);
    r = zeros<colvec>(N);
    betaMinusOne = 1.0;
    Q = zeros<mat>(N,iterationNumber);
    numberOfEigenvaluesConvergenceTest = 3;

    // Initial step
    rMinusOne.randu(N);
    beta(0) = norm(rMinusOne,2);
    q = rMinusOne/beta(0);
    if (tridiag != "true")
        alpha(1) = as_scalar(trans(q)*A*q);
    else{
        Aq(0) = A(0,0)*q(0) + A(0,1)*q(1);
        Aq(N-1) = A(N-1,N-2)*q(N-2) + A(N-1,N-1)*q(N-1);
        for (int row = 1; row < N-1; row++){
            Aq(row) = A(row,row-1)*q(row-1) + A(row,row)*q(row) + A(row, row+1)*q(row+1);
        }
        alpha(1) = as_scalar(trans(q)*Aq);
    }
    Q.col(0) = q;

    // Proceeding steps
    k = 1;
    lanczosIterationCounter = 1;
    eigenvaluesOld.zeros(5); // Fix if get in variable for how often calculate eigenvalues
    while (beta(k-1) != 0.0 && k < iterationNumber && eigenvaluesConverged != "true"){
        r = A*q -alpha(k)*q - beta(k-1)*qOld;
        qOld = q;
        beta(k) = norm(r,2);
        q = r/beta(k);
        if (tridiag != "true")
            alpha(k+1) = as_scalar(trans(q)*A*q);
        else{
            Aq(0) = A(0,0)*q(0) + A(0,1)*q(1);
            Aq(N-1) = A(N-1,N-2)*q(N-2) + A(N-1,N-1)*q(N-1);
            for (int row = 1; row < N-1; row++){
                Aq(row) = A(row,row-1)*q(row-1) + A(row,row)*q(row) + A(row, row+1)*q(row+1);
            }
            alpha(k+1) = as_scalar(trans(q)*Aq);
        }
        Q.col(k) = q;

        k += 1;

        lanczosIterationCounter += 1;
        if (lanczosIterationCounter == 5){
            lanczosIterationCounter = 0;
            *stopIteration = k;
            QforEigenvalue = Q.cols(0,k-1);
            T = trans(QforEigenvalue)*A*QforEigenvalue; // Could be done more efficient since A tridiagonal
            eigenvalues.zeros(k);
            eig_sym(eigenvalues, T);

            eigenvaluesDifference.zeros(numberOfEigenvaluesConvergenceTest);
            for (int j=0; j < numberOfEigenvaluesConvergenceTest; j++){
                eigenvaluesDifference(j) = eigenvalues(j)/eigenvaluesOld(j)-1.;
            }
            //eigenvaluesDifference = eigenvalues.head(3)-eigenvaluesOld.head(3);
            normMinimumEigenvalues = norm(eigenvaluesDifference, "inf");
            if ( fabs(normMinimumEigenvalues) < LanczosIterationTolerance )
                eigenvaluesConverged = "true";
            //eigenvalues.head(3).print("Eigenvalues during run: ");
            eigenvaluesOld = eigenvalues;
        }
    }
}
