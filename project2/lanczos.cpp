#include "lanczos.h"
#include "jacobi.h"

void lanczos(colvec &eigenvalues, mat &A, colvec &alpha, colvec &beta, mat &QforEigenvalue, int N, int iterationNumber, string tridiag, string eigenvalueSolver, int *stopIteration){
    int k, lanczosIterationCounter, numberOfEigenvaluesConvergenceTest;
    colvec r, rMinusOne, q, qOld,  Aq, eigenvaluesOld, eigenvaluesDifference;
    double betaMinusOne, normMinimumEigenvalues, LanczosIterationTolerance;
    mat T, Q;
    string eigenvaluesConverged;

    LanczosIterationTolerance = 1.e-03;
    Aq = zeros<colvec>(N);
    beta.zeros(iterationNumber+1);
    alpha.zeros(iterationNumber+1);
    qOld.zeros(N);
    r = zeros<colvec>(N);
    betaMinusOne = 1.0;
    Q = zeros<mat>(N,iterationNumber+1);
    numberOfEigenvaluesConvergenceTest = 3;

    // Initial step
    //rMinusOne.randu(N);
    //beta(0) = norm(rMinusOne,2);
    //q = rMinusOne/beta(0);

    Q.col(0) = qOld;
    q.randu(N);
    q = q/norm(q,2);
    r = q;
    Q.col(1) = q;

    beta(0)  = 1.;
    k = 0;

    q = r/beta(0);
    k = k+1;
    if (tridiag != "true")
        alpha(k) = as_scalar(trans(q)*A*q);
    else{
        Aq(0) = A(0,0)*q(0) + A(0,1)*q(1);
        Aq(N-1) = A(N-1,N-2)*q(N-2) + A(N-1,N-1)*q(N-1);
        for (int row = 1; row < N-1; row++){
            Aq(row) = A(row,row-1)*q(row-1) + A(row,row)*q(row) + A(row, row+1)*q(row+1);
        }
        alpha(k) = as_scalar(trans(q)*Aq);
    }
    //Q.col(0) = q/norm(q,2);
    r = A*q -alpha(k)*q - beta(k-1)*qOld;
    qOld = q;
    beta(k) = norm(r,2);

    // Proceeding steps
    lanczosIterationCounter = 1;
    eigenvaluesOld.zeros(5); // Fix if get in variable for how often calculate eigenvalues

    while (beta(k-1) != 0.0 && k < iterationNumber && eigenvaluesConverged != "true"){
        q = r/beta(k);
        QforEigenvalue = Q.cols(0,k); //For orthogonalization
        projection(QforEigenvalue, q, N, k); // Orthonormalization of q
        Q.col(k+1) = q;
        QforEigenvalue = Q.cols(0,k+1);
        k+=1;

        if (tridiag != "true")
            alpha(k) = as_scalar(trans(q)*A*q);
        else{
            Aq(0) = A(0,0)*q(0) + A(0,1)*q(1);
            Aq(N-1) = A(N-1,N-2)*q(N-2) + A(N-1,N-1)*q(N-1);
            for (int row = 1; row < N-1; row++){
                Aq(row) = A(row,row-1)*q(row-1) + A(row,row)*q(row) + A(row, row+1)*q(row+1);
            }
            alpha(k) = as_scalar(trans(q)*Aq);
        }

        r = A*q -alpha(k)*q - beta(k-1)*qOld;
        qOld = q;

        beta(k) = norm(r,2);

        lanczosIterationCounter += 1;
        if (lanczosIterationCounter == 5){
            lanczosIterationCounter = 0;
            *stopIteration = k;
            //QforEigenvalue = Q.cols(0,k-1);

            T = trans(QforEigenvalue)*A*QforEigenvalue; // Could be done more efficient since A tridiagonal
            eigenvalues.zeros(k);
            eig_sym(eigenvalues, T);

            eigenvaluesDifference.zeros(numberOfEigenvaluesConvergenceTest);
            for (int j=0; j < numberOfEigenvaluesConvergenceTest; j++){
                eigenvaluesDifference(j) = fabs(eigenvalues(j)/eigenvaluesOld(j)-1.);
            }
            //eigenvaluesDifference = eigenvalues.head(3)-eigenvaluesOld.head(3);
            normMinimumEigenvalues = norm(eigenvaluesDifference, "inf");
            if ( normMinimumEigenvalues < LanczosIterationTolerance )
                eigenvaluesConverged = "true";
            //eigenvalues.head(3).print("Eigenvalues during run: ");
            eigenvaluesOld = eigenvalues;
        }
    }

    //eigenvaluesDifference.print("Eigenvalue difference lanczos: ");

}

void gramSmith(mat &QforEigenvalue, int N, int columnNumber){
    //QforEigenvalue.cols(k)  =
    mat QOrtho = zeros<mat>(N,columnNumber);
    colvec sumVector;
    sumVector.zeros(N);


    for (int vectorNumber = 0; vectorNumber < columnNumber; vectorNumber++){
        for (int k = 0; k < vectorNumber ; k++){
            sumVector = as_scalar(trans(QOrtho.col(k))*QforEigenvalue.col(vectorNumber))/norm(QOrtho.col(k))*QOrtho.col(k);
        }
        QOrtho.col(vectorNumber) = QforEigenvalue.col(vectorNumber) - sumVector;
        QOrtho.col(vectorNumber) = QOrtho.col(vectorNumber)/norm(QOrtho.col(vectorNumber), 2);
    }
    QforEigenvalue = QOrtho;
}

void projection(mat &QforEigenvalue, colvec q, int N, int columnNumber){
    colvec sumVector;
    sumVector.zeros(N);

    for (int k = 0; k < columnNumber ; k++){
        sumVector = as_scalar(trans(QforEigenvalue.col(k))*q)/norm(QforEigenvalue.col(k))*QforEigenvalue.col(k);
    }
    q = q - sumVector;
    q = q/norm(q, 2);
}
