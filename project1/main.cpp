/* Project1 fys4150 2017.
Peter Even Killingstad and Karl Jacobsen
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
#include "lib.h"

using namespace std;


void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& lowerDiagonal, double& diagonal, double& upperDiagonal, int , char** argv);
void generate_tridiagonal_matrix(int,double, double, double, double **);
void generate_right_hand_side( int, double *, double *);
void gaussianTridiagonalSolver(double ** computed_tridiagonal_matrix, double * computed_right_hand_side, double * computed_numerical_solution, int N);
void gassianTridiagonalSymmetricSolver( double * computed_right_hand_side, double * computed_numerical_solution, int N);
void generate_exact_solution(int, double *);
void calculate_error(double *, double *, double *, double *, int);
void output_scalars( double, double, double, double, double);
void output_vectors( double *, int, int, string);


ofstream ofile1; // File for scalars
ofstream ofile2; // File for vectors


int main(int argc, char *argv[]){
  string outfile_name;
  int number_of_simulations, N, amplificationFactor, *indxLu;
  double lowerDiagonal, diagonal,  upperDiagonal, h_step, computed_error, time_used, L2Norm, dLu;
  double ** computed_tridiagonal_matrix, ** LU_Lower, ** LU_Upper;
  double * computed_numerical_solution, * computed_right_hand_side, * computed_exact_solution;
  string outfile_name_scalars, outfile_name_computed_numerical,outfile_name_computed_exact;
  clock_t start, finish;

  initialize(outfile_name, number_of_simulations, amplificationFactor, N, lowerDiagonal, diagonal,  upperDiagonal, argc, argv);
  outfile_name_scalars = (outfile_name) + string("_scalars")+string(".csv");
  outfile_name_computed_numerical = (outfile_name) + string("_numerical");
  outfile_name_computed_exact = (outfile_name) + string("_exact");
  ofile1.open(outfile_name_scalars);
  ofile1 << "log_h,log_rel_error,l2norm,time_used,logTimeUsed" << endl;

  // Running solvers for different mesh sizes, N
  for (int simulation_number = 0; simulation_number < number_of_simulations; simulation_number++){
    computed_tridiagonal_matrix = new double*[N];
    computed_right_hand_side = new double[N];
    computed_numerical_solution = new double[N];
    computed_exact_solution = new double[N];
    indxLu = new int[N];

    for (int i = 0; i < N; i++ ){
      computed_tridiagonal_matrix[i] = new double[N];
    }

    generate_tridiagonal_matrix(N, lowerDiagonal, diagonal,  upperDiagonal, computed_tridiagonal_matrix);
    generate_right_hand_side(N, computed_right_hand_side, &h_step);

    start = clock();
    if (outfile_name == "gaussianTridiagonal")
        gaussianTridiagonalSolver(computed_tridiagonal_matrix, computed_right_hand_side, computed_numerical_solution, N);
    else if (outfile_name == "gaussianTridiagonalSymmetric")
        gassianTridiagonalSymmetricSolver(computed_right_hand_side, computed_numerical_solution, N);
    else if (outfile_name=="luLib" ){
            ludcmp(computed_tridiagonal_matrix, N, indxLu, &dLu);
            for ( int i = 0; i < N; i++) computed_numerical_solution[i] = computed_right_hand_side[i];
            lubksb(computed_tridiagonal_matrix, N, indxLu, computed_numerical_solution);
    }
    else {
        cout << "define solver type. " << endl;
        exit(0);
    }

    finish = clock();
    time_used = (double)((finish - start)/double(CLOCKS_PER_SEC));
    double logTimeUsed = log10(time_used);

    generate_exact_solution(N, computed_exact_solution);
    calculate_error(computed_numerical_solution, computed_exact_solution, &computed_error, &L2Norm, N);
    output_scalars( L2Norm, computed_error, h_step, time_used, logTimeUsed);
    output_vectors( computed_exact_solution, simulation_number, N, outfile_name_computed_exact);
    output_vectors( computed_numerical_solution, simulation_number, N, outfile_name_computed_numerical);
    if (simulation_number < number_of_simulations -1)
        N *= amplificationFactor;
    }

  ofile1.close();
  ofile2.close();

  for (int i = 0; i < N; i++ )
    delete[] computed_tridiagonal_matrix[i];

  delete [] computed_tridiagonal_matrix;
  delete [] computed_numerical_solution;
  delete [] computed_right_hand_side;
  delete [] computed_exact_solution;
  delete [] indxLu;

  return 0;
}

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& lowerDiagonal, double& diagonal, double& upperDiagonal, int argc, char** argv )
{
    if( argc<= 1){
      cout << "Insert: outfile-name, number of simulations, amplification factor, start dimension, tridiagonal elements" << endl;
      exit(1);
    }
    else{
      outfile_name=argv[1];
    }
    number_of_simulations = atoi(argv[2]);
    amplificationFactor = atoi(argv[3]);
    N = atoi(argv[4]);
    lowerDiagonal= atof(argv[5]); // Lower diagonal
    diagonal= atof(argv[6]); // Digonal
    upperDiagonal= atof(argv[7]); // Upper diagonal
}

void generate_tridiagonal_matrix(int N,double lowerDiagonal, double diagonal, double  upperDiagonal, double **computed_tridiagonal_matrix ){
  computed_tridiagonal_matrix[0][1] =  upperDiagonal;
  computed_tridiagonal_matrix[0][0] = diagonal;
  computed_tridiagonal_matrix[N-1][N-2] = lowerDiagonal;
  computed_tridiagonal_matrix[N-1][N-1] = diagonal;

  for (int k = 1; k < N-1; k++){
    computed_tridiagonal_matrix[k][k] = diagonal;
    computed_tridiagonal_matrix[k][k-1] = lowerDiagonal;
    computed_tridiagonal_matrix[k][k+1] =  upperDiagonal;
    }
}

void generate_right_hand_side( int N, double *computed_right_hand_side, double *h_step){
    // Adjusted right hand side with h^2
    double h = 1./(N+1);
    *h_step = h;
    for ( int i=0; i < N; i++ ){
        double x = (i+1)*h;
        double source_term = 100.*exp(-10.*x);
        computed_right_hand_side[i] = h*h*source_term;
    }
}

void gaussianTridiagonalSolver(double ** computed_tridiagonal_matrix, double * computed_right_hand_side, double * computed_numerical_solution, int N){
    double multiplicationFactor;
    for (int i = 1; i < N; i++){
        multiplicationFactor = computed_tridiagonal_matrix[i][i-1]/computed_tridiagonal_matrix[i-1][i-1];
        computed_tridiagonal_matrix[i][i] += - multiplicationFactor*computed_tridiagonal_matrix[i-1][i];
        computed_right_hand_side[i] += - multiplicationFactor*computed_right_hand_side[i-1];
    }
    for(int i = N-1; i > -1; i--){
        if(i == N-1)
            computed_numerical_solution[i] = computed_right_hand_side[i]/computed_tridiagonal_matrix[i][i];
        else
            computed_numerical_solution[i] = (computed_right_hand_side[i] - computed_tridiagonal_matrix[i][i+1]*computed_numerical_solution[i+1])/computed_tridiagonal_matrix[i][i];
    }
}

void gassianTridiagonalSymmetricSolver( double * computed_right_hand_side, double * computed_numerical_solution, int N){
    double diagonal;
    double * f; // Temererary right hand side for forward substitution
    f = new double[N];
    f[0] = computed_right_hand_side[0];
    diagonal = 2.;
    for (int i=1; i < N; i++){
        f[i] = computed_right_hand_side[i] + double(i)/(i+1.)*f[i-1];
    }
    diagonal = ((N-1)+2.)/((N-1)+1.);
    computed_numerical_solution[N-1] = f[N-1]/diagonal;
    for (int i = N-1; i >0; i--){
        computed_numerical_solution[i-1] = double(i)/(i+1.)*(f[i-1] + computed_numerical_solution[i]);
    }
    delete [] f;
}

void generate_exact_solution(int N, double *computed_exact_solution){
    for ( int i=0; i < N; i++ ){
        double h = 1./(N+1);
        double x = (i+1)*h;
        computed_exact_solution[i] = 1. - (1. - exp(-10.))*x - exp(-10.*x);
    }
}

void calculate_error(double *computed_numerical_solution, double *computed_exact_solution, double *computed_error, double * L2Norm, int N){
    // Calculates sup-norm for relative error
    double temp_relative_error;
    *computed_error = log10(fabs((computed_numerical_solution[0] - computed_exact_solution[0])/computed_exact_solution[0]));
    for (int i = 1; i < N; i++){
        temp_relative_error = log10(fabs((computed_numerical_solution[i] - computed_exact_solution[i])/computed_exact_solution[i]));
        if(temp_relative_error > *computed_error){
            *computed_error = temp_relative_error;
        }
    }
    double L2NormTemp;
    for (int i = 0; i < N; i++){
        L2NormTemp = (computed_numerical_solution - computed_exact_solution)*(computed_numerical_solution - computed_exact_solution);
        *L2Norm += + L2NormTemp;
    }
    *L2Norm = sqrt(*L2Norm/double(N));
}

void output_scalars( double L2Norm, double computed_error, double h_step, double time_used, double logTimeUsed){
  ofile1 << setiosflags(ios::showpoint | ios::uppercase);
  ofile1 << setw(15) << setprecision(16) << log10(h_step) << ", ";
  ofile1 << setw(15) << setprecision(16) << computed_error << ", ";
  ofile1 << setw(15) << setprecision(16) << L2Norm << ", ";
  ofile1 << setw(15) << setprecision(16) << time_used << ", ";
  ofile1 << setw(15) << setprecision(16) << logTimeUsed << endl;
}

void output_vectors( double *output_array, int simulation_number, int N, string outfile_name){
  string filename = outfile_name + to_string(simulation_number+1)+string(".csv");
  ofile2.open(filename);
  ofile2 << "simulation_number_" << simulation_number+1 << endl;
  ofile2 << setiosflags(ios::showpoint | ios::uppercase);
  for (int i = 0; i<N; i++){
    ofile2 << setw(15) << setprecision(16) << output_array[i] << endl;;
  }
  ofile2 << setw(15) << setprecision(16) << endl;
  ofile2.close();
}
