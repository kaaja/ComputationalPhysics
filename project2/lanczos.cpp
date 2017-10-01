#include "lanczos.h"
#include "jacobi.h"

void lanczos(colvec &eigenvalues, mat &A, colvec &alpha, colvec &beta, mat &Q, int N, int iterationNumber, string tridiag, string eigenvalueSolver){
    int k;
    colvec r, rMinusOne, q, qOld,  Aq;
    double betaMinusOne;
    mat T;

    Aq = zeros<colvec>(N);
    beta.zeros(iterationNumber);
    alpha.zeros(iterationNumber+1);
    qOld.zeros(N);
    r = zeros<colvec>(N);
    betaMinusOne = 1.0;
    Q = zeros<mat>(N,iterationNumber);

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
    while (beta(k-1) != 0.0 && k < iterationNumber){
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
    }

    // Eigenvalues
    T = trans(Q)*A*Q; // Inefficient if have eigenvalue solver for tridiagonal matrices that takes diagonal and off-diagonals as input
    if (eigenvalueSolver == "jacobi"){
        double tolerance = 1.e-6; // Take in as function arguement if time.
        int maxIterations = 1e8;
        int counter = 0;
        mat v;
        v = eye<mat>(k,k); // Equal to actual length (takes into account if while loop terminates)
        eigenvalues  = zeros<colvec>(k);
        jacobi(T, eigenvalues, tolerance, maxIterations, N, &counter, v);
    }
    else if (eigenvalueSolver == "armadillo")
        eig_sym(eigenvalues, T);
}
