/* Project1 fys4150 2017. 
    

*/

using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

void initialize( int *,int *, double *, double *, double *);
void generate_tridiagonal_matrix(int,double, double, double, double **);
void generate_LU(int, double **, double **, double **);
void generate_right_hand_side( int, double *, double *);
void LU_tridiagonal_solver(double **, double **, double *, double *, int );
void generate_exact_solution(int, double *);
void calculate_error(double *, double *, double *, int);
void output_scalars( double, double);
void output_vectors( double *, double *, int, int);


ofstream ofile1;
ofstream ofile2;


int main(int argc, char *argv[]){
  // Declaration variables
  char *outfile_name;
  int number_of_simulations, N;
  double a, b, c, h_step, computed_error;
  double ** computed_tridiagonal_matrix, ** LU_Lower, ** LU_Upper;
  double * computed_numerical_solution, * computed_right_hand_side, * computed_exact_solution;
  string outfile_name_scalars, outfile_name_vectors;

  if( argc <= 1){
    cout << "Give a name to the output_scalars file after the program name" << endl;
    exit(1);
  }
  else{
    outfile_name=argv[1];
  }

  outfile_name_scalars = string(outfile_name) + string("_scalars");
  outfile_name_vectors = string(outfile_name) + string("_vectors");

  // Read variables from command line
  initialize(&number_of_simulations, &N, &a, &b, &c);

  // Pointer/arrays
  //h_step = new double[number_of_simulations];

  ofile1.open(outfile_name_scalars);
  ofile1 << "log_h,log_rel_error" << endl; // DO NOT USE WHITESPACE BETWEEN VAR-NAMES. Pandas does not handle it.
  ofile2.open(outfile_name_vectors);

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
    LU_tridiagonal_solver(LU_Lower, LU_Upper, computed_numerical_solution, computed_right_hand_side, N);
    generate_exact_solution(N, computed_exact_solution);
    calculate_error(computed_numerical_solution, computed_exact_solution, &computed_error, N);
    output_scalars( computed_error, h_step);
    output_vectors( computed_exact_solution, computed_numerical_solution, simulation_number, N);
    //output_scalars(h_step, computed_derivative, x, number_of_steps, outfile_name_scalars);
    if (simulation_number < number_of_simulations -1)
        N *= 10;
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

void initialize( int *number_of_simulations, int *N, double *a, double *b, double *c )
{
  printf("Read in from screen:number of simulations, matrix dimension, lower diagonal, diagonal, upperdiagonal \n");
  scanf(" %d %d %lf %lf %lf",number_of_simulations, N, a, b, c);
  return;
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
        double source_term = 100*exp(-10*x);
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

void generate_exact_solution(int N, double *computed_exact_solution){
    for ( int i=0; i < N; i++ ){
        double h = 1./(N+1);
        double x = (i+1)*h;
        computed_exact_solution[i] = 1 - (1 - exp(-10))*x - exp(-10*x);
    }
}

void calculate_error(double *computed_numerical_solution, double *computed_exact_solution, double *computed_error, int N){
    double temp_relative_error;
    *computed_error = log10(fabs(computed_numerical_solution[0] - computed_exact_solution[0])/fabs(computed_exact_solution[0]));
    for (int i = 1; i < N; i++){
        temp_relative_error = log10(fabs(computed_numerical_solution[i] - computed_exact_solution[i])/computed_exact_solution[i]);
        if(temp_relative_error > *computed_error)
            *computed_error = temp_relative_error;
    }
}

void output_scalars( double computed_error, double h_step){
  ofile1 << setiosflags(ios::showpoint | ios::uppercase);
  ofile1 << setw(15) << setprecision(8) << log10(h_step) << ", ";
  ofile1 << setw(15) << setprecision(8) << computed_error << endl;
}

void output_vectors( double *computed_exact_solution, double *computed_numerical_solution, int simulation_number, int N){
  ofile2 << "exact_simulation_number_" << simulation_number+1 << ", " << "numerical_simulation_number_" <<simulation_number+1 << endl; // DO NOT USE WHITESPACE BETWEEN VAR-NAMES. Pandas does not handle it.
  ofile2 << setiosflags(ios::showpoint | ios::uppercase);
  for (int i = 0; i<N; i++){
    ofile2 << setw(15) << setprecision(8) << computed_exact_solution[i] << ", ";
    ofile2 << setw(15) << setprecision(8) << computed_numerical_solution[i] << endl;
  }
}


