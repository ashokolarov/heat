#include "linalg.h"
#include "heat.h"
#include <iostream>

int main(){
    Heat<1, double> H(0.3125, 25, 0.25);
    Vector<double> exact = H.exact(10.5);
    Vector<double> approx = H.solve(10.5);
    Vector<double> diff = exact-approx;
    double norm = mag(diff);
    std::cout << norm;
}

