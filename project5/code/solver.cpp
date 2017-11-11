#include "solver.h"

Solver::Solver(double dt_, double dx_, double theta_, double T_)
{
    dt = dt_;
    dx = dx_;
    theta = theta_;
    T = T_;
    Nt = (int) (round(T/dt));
    Nx = (int) (round(1./dx+1));
}

void Solver::solve(){

    alpha = dt/pow(dx,2);
    offDiagonal = 2.*alpha*(1. - theta);
    diagonal = 2.*(1 - 2.*alpha*(1 - theta));

    double * u, * computed_right_hand_side, *uOld;

    u = new double[Nx];
    uOld = new double[Nx];
    computed_right_hand_side = new double[Nx];


    // IC
    for (int i =0; i < Nx; i++){
        uOld[i] = -dx*i;
        u[i] = 0.;
        computed_right_hand_side[i] = 0.;
    }

    for (int t = 1; t <= Nt; t++){
        generate_right_hand_side( computed_right_hand_side, uOld, offDiagonal, diagonal, Nx);
        gassianTridiagonalSymmetricSolver( computed_right_hand_side,  u, uOld, offDiagonal, diagonal, Nx);
        // Change of variables
        double uANalytical;
        for(int i = 0; i < Nx; i++){
            u[i] += +dx*i;
            for (int j = 0; j < Nx*100; j++){
                uANalytical += exp(-pow(j*3.14, 2)*0.1)*2./3.14*pow((-1.), j)*sin(j*3.14*dx*i) + dx*i;
            }
            cout << u[i] << " analytic " << uANalytical << endl;
        }
    }
    delete [] computed_right_hand_side;
    delete [] u;
    delete [] uOld;
}

void Solver::generate_right_hand_side(double *computed_right_hand_side , double *uOld, double offDiagonal, double diagonal, int Nx){
    computed_right_hand_side[1] = offDiagonal*uOld[0] + offDiagonal*uOld[1];
    for ( int i=2; i < Nx-2; i++ ){
        computed_right_hand_side[i] = offDiagonal*uOld[i-1] + diagonal*uOld[i] + offDiagonal*uOld[i+1];
    }
    computed_right_hand_side[Nx-2] = offDiagonal*uOld[Nx-3] + diagonal*uOld[Nx-2];
    computed_right_hand_side[0] = 0.;
    computed_right_hand_side[Nx-1] = 0.;
}

void Solver::gassianTridiagonalSymmetricSolver( double * computed_right_hand_side, double * u, double * uOld, double offDiagonal, double diagonal, int Nx){
    double * f; // Temererary right hand side for forward substitution
    f = new double[Nx];
    f[0] = 0.;
    f[Nx-1] = 0.;
    f[1] = computed_right_hand_side[1];
    double * diagonalTilde; // Temererary right hand side for forward substitution
    diagonalTilde = new double[Nx];
    diagonalTilde[0] = 0.;
    diagonalTilde[Nx-1] = 0.;
    diagonalTilde[1] = diagonal;
    for (int i = 2; i < Nx-1; i++){
        diagonalTilde[i] = diagonal - pow((-offDiagonal), 2)/diagonalTilde[i-1];
        f[i] = computed_right_hand_side[i] -offDiagonal*f[i-1]/diagonalTilde[i-1];
    }

    // Backward substitution
    u[0] = 0.;
    u[Nx-1] = 0.;
    u[Nx-2] = f[Nx-2]/diagonalTilde[Nx-2];
    for (int i = Nx-2; i >1; i--){
        u[i-1] = f[i-1] - offDiagonal*u[i]/diagonalTilde[i-1];
    }

    delete [] f;
    delete [] diagonalTilde;
    //delete [] uOld;
    for (int i = 0; i < Nx; i++) uOld[i] = u[i];
}
