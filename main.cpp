#include "linalg.h"
#include "heat.h"

int main(){
    Heat<2, double> H(0.15, 35, 1, 0.01);

    Vector<double> exact = H.exact();
    Vector<double> approx = H.solve();
    Vector<double> diff = exact-approx;

}

