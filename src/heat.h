#ifndef HEAT_HEAT_H
#define HEAT_HEAT_H

#include "linalg.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

static const std::string dir = "/home/alexshokolarov/Documents/Projects/heat/plot/";

static inline double initial(double x){
    return 10*std::sin(x);
}

template<typename T>
int cg(const Sparse<T>& A,
       const Vector<T>& b,
       Vector<T>&       x,
       T                tol     = (T)1e-8,
       int              maxiter = 100)
{

    Vector<T> rk = b - A * x;
    Vector<T> rk_1(x.size());
    auto pk = rk;

    for (int i=0; i < maxiter; i++){
        auto alpha = dot(rk, rk) / dot(pk, (A * pk));
        x = x + alpha * pk;
        rk_1 = rk - alpha * (A * pk);

        if (dot(rk_1, rk_1) < (tol * tol)){
            return i;
        }

        auto beta = dot(rk_1, rk_1) / dot(rk, rk);
        pk = rk_1 + beta * pk;
        rk = rk_1;
    }

    return -1;
}

template <int n, typename T>
class Heat
{
public:
    const T alpha;
    const int m;
    const T Ttot;
    const T dt;
    Sparse<T> M;
    Vector<T> w0;

    std::ofstream temp;
    std::ofstream xcor;
    std::ofstream time;

    Heat(T _alpha, int _m, T Ttot, T _dt) :
            alpha(_alpha), m(_m), Ttot(Ttot), dt(_dt), M(std::pow(_m, n), std::pow(_m,n)), w0(std::pow(_m,n))
    {
        int totp = std::pow(m, n);
        T dx = 1.0 / (1.0 + m);
        T coeff = alpha * (dt/dx/dx);

        time.open(dir + "time.txt", std::ios::ate);
        temp.open(dir + "temp.txt", std::ios::ate);
        xcor.open(dir + "xcor.txt", std::ios::ate);

        for (int i=0; i<totp; i++)
        {
            for (int j=0; j<totp; j++)
            {
                int D = 0;

                for (int k=0; k<n; k++)
                {
                    if (i==j)
                    {
                        D -= 2;
                    }

                    else if (std::abs(i-j) == std::pow(m,k))
                    {
                        D += 1;
                    }
                }

                if (i==j)
                {
                    M[{i,j}] = 1 - coeff * D;
                }
                else
                {
                    M[{i,j}] = - coeff * D;
                }

            }
        }

        w0 = Vector<T>::ones(w0.size());

        for (int i=0; i<totp; i++)
        {
            for (int k=0; k<n; k++)
            {
                T coor = (std::fmod(i / std::pow(m,k), m) + 1) / (m + 1);
                w0[i] *= initial(M_PI * coor);

                xcor << coor << " ";
                temp << w0[i] << " ";
            }

        }

        xcor.close();
        temp << '\n';
    }

    Vector<T> exact()
    {
        Vector<T> sol(std::pow(m,n));
        for (int i=0; i<std::pow(m,n); i++)
        {
            sol[i] = std::exp(-n*M_PI*M_PI*alpha*Ttot) * w0[i];
        }

        return sol;
    }

    Vector<T> solve()
    {
        Vector<T> sol(std::pow(m,n));
        Vector<T> interm(std::pow(m,n));
        interm = w0;
        int N = static_cast<int>(Ttot/dt);
        time << '0' << ' ';

        for (int i=1; i<=N; i++)
        {
            time << i*dt << ' ';
            int cg_status = cg<T>(M, interm, sol);
            if (cg_status == -1)
            {
                throw std::runtime_error("Max iteration of cg reached");
            }

            for (int j=0; j<sol.size(); j++){
                temp << sol[j] << ' ';
            }

            temp << "\n";
            interm = sol;
        }

        temp.close();
        time.close();
        return sol;
    }
};

template<typename T>
[[maybe_unused]] double mag(const Vector<T> V)
{
    double mag = 0;
    for (int i=0; i<V.size(); i++)
    {
        mag += V[i] * V[i];
    }
    return std::pow(mag, 0.5);
}

#endif //HEAT_HEAT_H
