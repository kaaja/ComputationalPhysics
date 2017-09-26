/* This is a rewrite of the python programs in Kiusalaas (2014) "Numerical methods in engineering with python" for finding eigenvalues
 of symmetric tridiagonal matrices  */

#include "eigenvalueBisection.h"

colvec sturmSeq(mat &A, double &lam, int N){
    colvec p = ones<colvec>(N+1);
    p(1) = A(0,0) - lam;
    for (int i = 2; i < N+1; i++ ){
        if (A(i-1, i-2) == 0.0)
            A(i-1,i-2) = 1.0e-12;
        p(i) = (A(i-1,i-1)-lam)*p(i-1) - (A(i-1,i-2)*A(i-1,i-2))*p(i-2);
    }
    return p;
}

int numLambdas(colvec &p, int N){
    int sign;
    int signOld = 1;
    int numLam = 0.0;
    for (int i = 1; i < N+1; i++){
        if (p(i) > 0.0)
            sign = 1;
        else if (p(i) < 0.0)
            sign = -1;
        else
            sign = - signOld;
        if (sign*signOld < 0)
            numLam = numLam + 1;
        signOld = sign;
    }
    return numLam;
}

colvec gerschgorin(mat &A, int N){
    double lamMin, lamMax, lam;
    colvec lambdaExtremes;
    lambdaExtremes = zeros<colvec>(2);

    lamMin = A(0,0) - fabs(A(0,1));
    lamMax = A(0,0) + fabs(A(0,1));
    for (int i = 1; i < N-1; i++){
        lam = A(i, i) - fabs(A(i, i-1)) - fabs(A(i,i-1));
        if (lam < lamMin)
            lamMin = lam;
        lam = A(i,i) + fabs(A(i,i-1)) + fabs(A(i,i-1));
        if (lam > lamMax)
            lamMax = lam;
    }
    lam = A(N-1,N-1) - fabs(A(N-2,N-1));
    if (lam < lamMin)
        lamMin = lam;
    lam = A(N-1, N-1) + fabs(A(N-2, N-1));
    if (lam > lamMax)
        lamMax = lam;
    lambdaExtremes(0) = lamMin;
    lambdaExtremes(1) = lamMax;

    return lambdaExtremes;
}

colvec lamRange(mat &A,int N){
    colvec lambdaExtrmes, r, p;
    double lamMin,lamMax, lam, h;
    int numLam;
    lambdaExtrmes = zeros<colvec>(2);
    r = ones<colvec>(N+1);
    p = ones<colvec>(N+1);
    lambdaExtrmes = gerschgorin(A,N);
    lamMin = lambdaExtrmes(0); lamMax = lambdaExtrmes(1);
    r(0) = lamMin;
    // Search for eigenvalues in descending order
    for (int k = N; k > 0; k--){
        // First bisection of interval(lamMin,lamMax)
        lam = (lamMax + lamMin)/2.0;
        h = (lamMax - lamMin)/2.0;
        for (int i = 0; i < 1000; i++){
            // Find number of eigenvalues less than lam
            p = sturmSeq(A, lam, N);
            numLam = numLambdas(p,N);
            // Bisect again & find the half containing lam
            h = h/2.0;
            if (numLam < k)
                lam = lam + h;
            else if (numLam > k)
                lam = lam - h;
            else
                break;
        }
        // If eigenvalue located, change the upper limit
        // of search and record it in [r]
        lamMax = lam;
        r(k) = lam;
    }
    return r;
}

double f(mat &A, double eigenvalueGuess, int N){
    colvec p;
    p = zeros<colvec>(N+1);
    p = sturmSeq(A, eigenvalueGuess, N);
    return p(N);
}

colvec eigenvals3(mat &A, int N, double bisectionAccuracy, int max_iterations){
    colvec lam, r;
    lam = zeros<colvec>(N);
    r = zeros<colvec>(N+1);
    r = lamRange(A,N);
    for (int i = 0; i < N; i++){
        lam(i) = bisection(f, r(i), r(i+1), bisectionAccuracy, A, N, max_iterations );
    }
    return lam;
}
double bisection(double (*func)(mat &A, double eigenvalueGuess, int N), double x1, double x2, double xacc, mat &A, int N, int max_iterations){
    // From Morten HJ lecture notes
    int j;
    double dx, f, fmid, xmid, rtb;
    f = (*func)(A, x1, N);
    fmid = (*func)(A, x2, N);
    if(f*fmid >= 0.0) {
        cout << "\n\nError in function bisection():" << endl;
        cout << "\nroot in function must be within" << endl;
        cout << "x1 =" << x1 << "and x2 " << x2 << endl;
        exit(1);
    }
    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    for(j = 0; j < max_iterations; j++) {
        fmid = (*func)(A, xmid = rtb + (dx *= 0.5), N);
        if (fmid <= 0.0)
            rtb=xmid;
        if(fabs(dx) < xacc || fmid == 0.0)
            return rtb;
    }
    cout << "Error in the bisection:" << endl; // should never reach this point
    cout << "Too many iterations!" << endl;
}

