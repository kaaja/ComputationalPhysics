#include "lanczos.h"

void lanczos(mat &A, colvec &alpha, colvec &beta, mat &Q, int N, int iterationNumber){
    int k;
    colvec r, rMinusOne, qOld, q, Aq;
    mat I, T;
    double betaMinusOne, betaTest;

    T = zeros<mat>(iterationNumber, iterationNumber);
    Aq = zeros<colvec>(N);
    I = eye<mat>(N,N);
    beta.zeros(iterationNumber);
    //beta(0) = 1.0; // To pass first loop iteration
    betaTest = 1.0; // To pass first loop iteration
    alpha.zeros(iterationNumber);
    k = 0;
    qOld = zeros<colvec>(N);
    q.randu(N);
    r = zeros<colvec>(N);
    rMinusOne = q;
    betaMinusOne = 1.0;
    Q = zeros<mat>(N,iterationNumber);
    //Q.col(0) = qOld;
    mat Q2 = zeros<mat>(N,N);
    while (betaTest != 0.0 && k < iterationNumber){
        //cout << " k " << k << " beta(k-1) " << betaTest << endl;
        if (k == 0)
            q = rMinusOne/betaMinusOne;
        else
            q = r/beta(k-1);
        Q.col(k) = q;
        k += 1;
        /*
        Aq(0) = (A(0,0) + A(0,1))*q(0);
        for (int row = 1; row < N-1; row++){
            Aq(row) = (A(row,row-1) + A(row,row) + A(row, row+1))*q(row);
        }
        Aq(N-1) = (A(N-1,N-1) + A(N-1,N-2))*q(N-1);
        */
        Aq = A*q;
        alpha(k-1) = as_scalar(trans(q)*Aq);
        if (k == 1)
            r = Aq - alpha(k-1)*q - betaMinusOne*qOld;
        else
            r = Aq - alpha(k-1)*q - beta(k-2)*qOld;
        qOld = q;
        //beta(k-1) =  norm(r,2);
        beta(k-1) =  sqrt(as_scalar(trans(r)*r));
        betaTest = beta(k-1);
    }
    Q2 = trans(Q)*Q;
    T = trans(Q)*A*Q;
    alpha.print("alpha: ");
    beta.print("beta: ");
    Q2.print("trans(Q)*Q: ");
    T.print("T = trans(Q)*A*Q: ");
}
