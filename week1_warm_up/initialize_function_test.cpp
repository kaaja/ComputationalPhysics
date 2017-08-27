// Testing initialize function in warm-up exercise

using namespace std;

#include <iostream>

void initialize (double *, double *, int *);

int main()
{
  int number_of_steps;
  double x, initial_step;

  initialize( &x, &initial_step, &number_of_steps ); /* Need "&" since the input-variables are not defined as pointers */
  cout << x << initial_step << number_of_steps << endl;
  return 0;
}

void initialize( double *x, double *initial_step, int *number_of_steps )
{
  printf("Read in from screen:x,  initial step,  and number of steps\n");
  scanf("%lf %lf %d",x,initial_step, number_of_steps);
  return;
}
  
