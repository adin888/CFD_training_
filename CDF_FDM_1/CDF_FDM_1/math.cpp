#include"cfd.h"

using namespace std;

/*
Calculate L - 2 norm of a vector
*/
double compute_l2norm(int nx, vector<double> erro)
{
    double rms = 0.0;
    for (int i = 1; i < nx; i++)
    {
        rms = rms + pow(erro[i], 2);
    }
    rms = sqrt(rms / (nx - 1));
    return rms;
}
/*
Calculate right hand term of the inviscid Burgers equation
r = -udu/dx
*/
vector<double> rhs(int nx, double dx, double dt, vector<double> u, vector<double> r, double alpha)
{
    for (int i = 1; i < nx; i++)
    {
        r[i] = alpha * dt * (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (dx * dx);
    }
    return r;
}
/*
tridiagonal matrix algorithm(Thomas algorithm)
*/
vector<double> tdms(vector<double> a, vector<double> b, vector<double> c, vector<double> r, vector<double> u, int start, int end)
{

}