#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"


using namespace std;

int main(){
	double * a, * b;
	int N = 10;
	a = new double[N];
	b = new double[N]; 
	for ( int i = 0; i < N; i++){
		a[i] = 0;
		b[i] = i;
	}
	//a = b;
	for ( int i = 0; i < N; i++){
		a[i] = b[i];
	}

	for ( int i = 0; i < N; i++){
		b[i] = 2*i;
	}
	

	for ( int i = 0; i < N; i++)
		cout << a[i] << endl;

}
