#include "linalg.h"
#include "heat.h"

int main(){
    Heat<1, double> H(0.05, 100, 5, 0.005);

    Vector<double> exact = H.exact();
    Vector<double> approx = H.solve();
    Vector<double> diff = exact-approx;

}

