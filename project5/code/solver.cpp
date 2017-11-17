#include "solver.h"

Solver::Solver(double dt_, double dx_, double theta_, double T_, int Nx_, int Nt_)
{
    dt = dt_;
    dx = dx_;
    theta = theta_;
    T = T_;
    Nt = Nt_;
    Nx = Nx_;
}

mat Solver::solve(string outfileName_)
{
    outfileName = outfileName_;
    mat solutionMatrixU = zeros<mat>(Nx+1,Nt+1);

    alpha = dt/pow(dx,2);
    offDiagonalLhs = -2.*alpha*theta;
    diagonalLhs = 2.*(2.*alpha*theta + 1);
    offDiagonalRhs = 2.*alpha*(1. - theta);
    diagonalRhs = 2.*(1.-2.*alpha*(1.-theta));


    double * u, * computed_right_hand_side, *uOld;

    u = new double[Nx];
    uOld = new double[Nx];
    computed_right_hand_side = new double[Nx];


    // Initial conditions and declaration arrays
    for (int i =0; i < Nx; i++){
        uOld[i] = -dx*i;
        u[i] = 0.;
        computed_right_hand_side[i] = 0.;
    }

    for (int t = 1; t < Nt; t++){
        generate_right_hand_side( computed_right_hand_side, uOld, offDiagonalRhs, diagonalRhs, Nx);
        gassianTridiagonalSymmetricSolver( computed_right_hand_side,  u, uOld, offDiagonalLhs, diagonalLhs, Nx);
        // Change of variables
        double uANalytical;
        for(int i = 0; i < Nx; i++){
            u[i] += dx*i;
            solutionMatrixU(0,t+1) = t*dt;
            solutionMatrixU(i+1,0) = i*dx;
            solutionMatrixU(i+1,t+1) = u[i];
        }
    }
    solutionMatrixU.save(outfileName + "SolutionMatrixU.txt", raw_ascii);
    delete [] computed_right_hand_side;
    delete [] u;
    delete [] uOld;

    return solutionMatrixU;
}

void Solver::generate_right_hand_side(double *computed_right_hand_side , double *uOld, double offDiagonalRhs, double diagonalRhs, int Nx){
    computed_right_hand_side[1] = offDiagonalRhs*uOld[0] + offDiagonalRhs*uOld[1];
    for ( int i=2; i < Nx-2; i++ ){
        computed_right_hand_side[i] = offDiagonalRhs*uOld[i-1] + diagonalRhs*uOld[i] + offDiagonalRhs*uOld[i+1];
    }
    computed_right_hand_side[Nx-2] = offDiagonalRhs*uOld[Nx-3] + diagonalRhs*uOld[Nx-2];
    computed_right_hand_side[0] = 0.;
    computed_right_hand_side[Nx-1] = 0.;
}

void Solver::gassianTridiagonalSymmetricSolver(double * computed_right_hand_side, double * u, double * uOld, double offDiagonalLhs, double diagonalLhs, int Nx){
    double * f; // Temererary right hand side for forward substitution
    f = new double[Nx];
    f[0] = 0.;
    f[Nx-1] = 0.;
    f[1] = computed_right_hand_side[1];
    double * diagonalTilde; // Temererary right hand side for forward substitution
    diagonalTilde = new double[Nx];
    diagonalTilde[0] = 0.;
    diagonalTilde[Nx-1] = 0.;
    diagonalTilde[1] = diagonalLhs;
    for (int i = 2; i < Nx-1; i++){
        diagonalTilde[i] = diagonalLhs - pow(offDiagonalLhs, 2)/diagonalTilde[i-1];
        f[i] = computed_right_hand_side[i] -offDiagonalLhs*f[i-1]/diagonalTilde[i-1];
    }

    // Backward substitution
    u[0] = 0.;
    u[Nx-1] = 0.;
    u[Nx-2] = f[Nx-2]/diagonalTilde[Nx-2];
    for (int i = Nx-2; i >1; i--){
        u[i-1] = (f[i-1] - offDiagonalLhs*u[i])/diagonalTilde[i-1];
    }

    delete [] f;
    delete [] diagonalTilde;
    //delete [] uOld;
    for (int i = 0; i < Nx; i++) uOld[i] = u[i];
}

void Solver::calculate_error(mat computed_numerical_solution, mat computed_exact_solution, double *computed_error, int Nx_, int Nt_){
    // Calculates sup-norm for relative error
    double temp_relative_error;
    *computed_error = log2(fabs((computed_numerical_solution(2,2) - computed_exact_solution(2,2))
                                  /computed_exact_solution(2,2)));
    for (int tj = 3; tj < Nt_; tj++){
        for (int i = 3; i < Nx_; i++){
            temp_relative_error = log2(fabs((computed_numerical_solution(i, tj) - computed_exact_solution(i, tj))
                                            /computed_exact_solution(i, tj)));
            if(temp_relative_error > *computed_error) *computed_error = temp_relative_error;
        }
    }
}
