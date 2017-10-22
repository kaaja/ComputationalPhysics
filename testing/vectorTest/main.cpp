#include <iostream>
#include <vector>

using namespace std;

int main()
{
    vector<double> acc;
    acc.push_back(1.0);
    cout << "acc[0] " << acc[0] <<  endl;
    acc.push_back(2.0);
    cout << "acc[1] " << acc[1] <<  endl;
    cout << "Hello " << endl;
    double a = 3.;
    double b = 3.;
    a += -1.;
    b -= 1.;
    cout << "3 += -1 = " << a << endl;
    cout << "3 -= 1 = " << b << endl;
    return 0;
}
