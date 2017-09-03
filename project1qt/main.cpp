/* Project1 fys4150 2017.


*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"

using namespace std;


void initialize(string& outfile_name, int& number_of_simulations, int& N, double& a, double& b, double& c, int , char** argv);
void generate_tridiagonal_matrix(int,double, double, double, double **);
void generate_LU(int, double **, double **, double **);
void generate_right_hand_side( int, double *, double *);
void LU_tridiagonal_solver(double **, double **, double *, double *, int );
void gaussianTridiagonalSolver(double ** computed_tridiagonal_matrix, double * computed_right_hand_side, double * computed_numerical_solution, int N);
void gassianTridiagonalSymmetricSolver(double ** computed_tridiagonal_matrix, double * computed_right_hand_side, double * computed_numerical_solution, int N);
void generate_exact_solution(int, double *);
void calculate_error(double *, double *, double *, double *, int);
void output_scalars( double, double, double, double);
void output_vectors( double *, int, int, string);


ofstream ofile1;
ofstream ofile2;


int main(int argc, char *argv[]){
  // Declaration variables
  string outfile_name;
  int number_of_simulations, N;
  double a, b, c, h_step, computed_error, time_used, L2Norm;
  double ** computed_tridiagonal_matrix, ** LU_Lower, ** LU_Upper;
  double * computed_numerical_solution, * computed_right_hand_side, * computed_exact_solution;
  string outfile_name_scalars, outfile_name_computed_numerical,outfile_name_computed_exact;
  clock_t start, finish;

  // Read variables from command line
  initialize(outfile_name, number_of_simulations, N, a, b, c, argc, argv);
  outfile_name_scalars = (outfile_name) + string("_scalars");
  outfile_name_computed_numerical = (outfile_name) + string("_numerical");
  outfile_name_computed_exact = (outfile_name) + string("_exact");
  ofile1.open(outfile_name_scalars);
  ofile1 << "log_h,log_rel_error, l2norm, time_used" << endl; // DO NOT USE WHITESPACE BETWEEN VAR-NAMES. Pandas does not handle it.
  //ofile2.open(outfile_name_vectors);

  for (int simulation_number = 0; simulation_number < number_of_simulations; simulation_number++){
    computed_tridiagonal_matrix = new double*[N];
    LU_Lower = new double*[N];
    LU_Upper = new double*[N];
    computed_right_hand_side = new double[N];
    computed_numerical_solution = new double[N];
    computed_exact_solution = new double[N];

    for (int i = 0; i < N; i++ ){
      computed_tridiagonal_matrix[i] = new double[N];
      LU_Lower[i] = new double[N];
      LU_Upper[i] = new double[N];
    }

    generate_tridiagonal_matrix(N, a, b, c, computed_tridiagonal_matrix);
    generate_LU(N, computed_tridiagonal_matrix, LU_Lower, LU_Upper);
    generate_right_hand_side(N, computed_right_hand_side, &h_step);

    start = clock();
    if(outfile_name=="LU")
        LU_tridiagonal_solver(LU_Lower, LU_Upper, computed_numerical_solution, computed_right_hand_side, N);
    else if (outfile_name == "gaussianTridiagonal")
        gaussianTridiagonalSolver(computed_tridiagonal_matrix, computed_right_hand_side, computed_numerical_solution, N);
    else if (outfile_name == "gaussianTridiagonalSymmetric")
        gassianTridiagonalSymmetricSolver(computed_tridiagonal_matrix, computed_right_hand_side, computed_numerical_solution, N);
    else {
        cout << "define solver type. LU or gaussian" << endl;
        exit(0);
    }

    finish = clock();
    time_used = (double)((finish - start)/double(CLOCKS_PER_SEC));

    //gassianTridiagonalSymmetricSolver(computed_tridiagonal_matrix, computed_right_hand_side, computed_numerical_solution, N);

    //gaussianTridiagonalSolver(computed_tridiagonal_matrix, computed_right_hand_side, computed_numerical_solution, N);

    generate_exact_solution(N, computed_exact_solution);
    calculate_error(computed_numerical_solution, computed_exact_solution, &computed_error, &L2Norm, N);
    output_scalars( L2Norm, computed_error, h_step, time_used);
    output_vectors( computed_exact_solution, simulation_number, N, outfile_name_computed_exact);
    output_vectors( computed_numerical_solution, simulation_number, N, outfile_name_computed_numerical);
    if (simulation_number < number_of_simulations -1)
        N *= 2;
    }
  ofile1.close();
  ofile2.close();
  for (int i = 0; i < N; i++ ){
    delete[] computed_tridiagonal_matrix[i];
    delete[] LU_Lower[i];
    delete[] LU_Upper[i];
  }
  delete [] computed_tridiagonal_matrix;
  delete [] LU_Lower;
  delete [] LU_Upper;
  delete [] computed_numerical_solution;
  delete [] computed_right_hand_side;
  delete [] computed_exact_solution;

  return 0;
}

void initialize(string& outfile_name, int& number_of_simulations, int& N, double& a, double& b, double& c, int argc, char** argv )
{
    if( argc <= 1){
      cout << "Insert: outfile-name, number of simulations, start dimension, tridiagonal elements" << endl;
      exit(1);
    }
    else{
      outfile_name=argv[1];
    }
    number_of_simulations = atoi(argv[2]);
    N = atoi(argv[3]);
    a = atof(argv[4]);
    b = atof(argv[5]);
    c = atof(argv[6]);
}

void generate_tridiagonal_matrix(int N,double a, double b, double c, double **computed_tridiagonal_matrix ){
  computed_tridiagonal_matrix[0][1] = c;
  computed_tridiagonal_matrix[0][0] = b;
  computed_tridiagonal_matrix[N-1][N-2] = a;
  computed_tridiagonal_matrix[N-1][N-1] = b;

  for (int k = 1; k < N-1; k++){
    computed_tridiagonal_matrix[k][k] = b;
    computed_tridiagonal_matrix[k][k-1] = a;
    computed_tridiagonal_matrix[k][k+1] = c;
    }
  return;
}

void generate_LU(int N, double **computed_tridiagonal_matrix, double **LU_Lower, double **LU_Upper){
  for ( int i = 0; i< N; i++ ){
    LU_Lower[i][i] = 1;
  }

  LU_Upper[0][0] = computed_tridiagonal_matrix[0][0];
  //U[row-1][row-1] = A[row-1][row-1];
  for (int i = 0; i < N-1; i++ ){
    if ( i != 0){
      LU_Upper[i][i] = computed_tridiagonal_matrix[i][i] -LU_Lower[i][i-1]*LU_Upper[i-1][i];
    }
    LU_Lower[i+1][i] = computed_tridiagonal_matrix[i+1][i]/LU_Upper[i][i];
    LU_Upper[i][i+1] = computed_tridiagonal_matrix[i][i+1];
    }
  LU_Upper[N-1][N-1] = computed_tridiagonal_matrix[N-1][N-1] - LU_Lower[N-1][N-2]*LU_Upper[N-2][N-1];
}

void generate_right_hand_side( int N, double *computed_right_hand_side, double *h_step){
    double h = 1./(N+1);
    *h_step = h;
    for ( int i=0; i < N; i++ ){
        double x = (i+1)*h;
        double source_term = 100.*exp(-10.*x);
        computed_right_hand_side[i] = h*h*source_term;
    }
}

void LU_tridiagonal_solver(double **LU_Lower, double **LU_Upper, double *computed_numerical_solution, double *computed_right_hand_side, int N){
    double *y;
    y = new double[N];
    for(int i =0; i < N; i++){
        if(i==0)
            y[i] = computed_right_hand_side[i]/LU_Lower[i][i];
        else
            y[i] = computed_right_hand_side[i] - LU_Lower[i][i-1]*y[i-1];
    }

    for(int i = N-1; i > -1; i--){
        if(i == N-1)
            computed_numerical_solution[i] = y[i];
        else
            computed_numerical_solution[i] = y[i+1] - LU_Upper[i][i+1]*computed_numerical_solution[i+1]/LU_Upper[i][i];
    }
}

void gaussianTridiagonalSolver(double ** computed_tridiagonal_matrix, double * computed_right_hand_side, double * computed_numerical_solution, int N){
    double multiplicationFactor;
    for (int i = 1; i < N; i++){
        multiplicationFactor = computed_tridiagonal_matrix[i][i-1]/computed_tridiagonal_matrix[i-1][i-1];
        //computed_tridiagonal_matrix[i][i-1] = 0;
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

void gassianTridiagonalSymmetricSolver(double ** computed_tridiagonal_matrix, double * computed_right_hand_side, double * computed_numerical_solution, int N){
    double multiplicationFactor;
    for (int i = 1; i < N; i++){
        multiplicationFactor = computed_tridiagonal_matrix[i][i-1]/computed_tridiagonal_matrix[i-1][i-1];
        //computed_tridiagonal_matrix[i][i-1] = 0;
        computed_tridiagonal_matrix[i][i] += - multiplicationFactor*computed_tridiagonal_matrix[i-1][i];
        computed_right_hand_side[i] += - multiplicationFactor*computed_right_hand_side[i-1];
    }
    for(int i = N-1; i > -1; i--){
        if(i == N-1)
            computed_numerical_solution[i] = computed_right_hand_side[i]/computed_tridiagonal_matrix[i][i];
        else
            computed_numerical_solution[i] = (computed_right_hand_side[i] - computed_tridiagonal_matrix[0][1]*computed_numerical_solution[i+1])/computed_tridiagonal_matrix[i][i];
    }

}

void generate_exact_solution(int N, double *computed_exact_solution){
    for ( int i=0; i < N; i++ ){
        double h = 1./(N+1);
        double x = (i+1)*h;
        computed_exact_solution[i] = 1. - (1. - exp(-10.))*x - exp(-10.*x);
    }
}

void calculate_error(double *computed_numerical_solution, double *computed_exact_solution, double *computed_error, double * L2Norm, int N){
    double temp_relative_error;
    *computed_error = log10(fabs((computed_numerical_solution[0] - computed_exact_solution[0])/computed_exact_solution[0]));
    for (int i = 1; i < N; i++){
        temp_relative_error = log10(fabs((computed_numerical_solution[i] - computed_exact_solution[i])/computed_exact_solution[i]));
        if(temp_relative_error > *computed_error){
            *computed_error = temp_relative_error;
            cout << i << endl;
        }
    }
    double L2NormTemp;
    for (int i = 0; i < N; i++){
        L2NormTemp = (computed_numerical_solution - computed_exact_solution)*(computed_numerical_solution - computed_exact_solution);
        *L2Norm += + L2NormTemp;
    }
    *L2Norm = sqrt(*L2Norm/double(N));
}

void output_scalars( double L2Norm, double computed_error, double h_step, double time_used){
  ofile1 << setiosflags(ios::showpoint | ios::uppercase);
  ofile1 << setw(15) << setprecision(16) << log10(h_step) << ", ";
  ofile1 << setw(15) << setprecision(16) << computed_error << ", ";
  ofile1 << setw(15) << setprecision(16) << L2Norm << ", ";
  ofile1 << setw(15) << setprecision(16) << time_used << endl;
}

void output_vectors( double *output_array, int simulation_number, int N, string outfile_name){
  string filename = outfile_name + to_string(simulation_number+1);
  ofile2.open(filename);
  ofile2 << "simulation_number_" << simulation_number+1 << endl; // DO NOT USE WHITESPACE BETWEEN VAR-NAMES. Pandas does not handle it.
  ofile2 << setiosflags(ios::showpoint | ios::uppercase);
  for (int i = 0; i<N; i++){
    ofile2 << setw(15) << setprecision(16) << output_array[i] << endl;;
    //ofile2 << setw(15) << setprecision(8) << computed_numerical_solution[i] << endl;
  }
  ofile2 << setw(15) << setprecision(16) << endl;
  ofile2.close();
}
