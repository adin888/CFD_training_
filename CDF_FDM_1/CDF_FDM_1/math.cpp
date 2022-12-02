#include"cfd.h"

using namespace std;

/*
* Calculate L - 2 norm of a vector
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
* Tridiagonal matrix algorithm(Thomas algorithm) Au=d, A is a matrix containing a_i,b_i,c_i
*/
vector<double> tdms(vector<double> a, vector<double> b, vector<double> c, vector<double> d, int start, int end)
{
    vector<double> q(end);                    //Storage of superdiagonal array
    vector<double> u(end + 1, 0.0);
    double b_tem = b[start];                  //Temporary storage of b
    u[start] = d[start] / b_tem;
    /*
    * Forward elimination
    */
    for (int i = start + 1; i < end; i++)
    {
        q[i] = c[i - 1] / b_tem;
        b_tem = b[i] - a[i] * q[i];
        u[i] = (d[i] - a[i] * u[i - 1]) / b_tem;
    }
    /*
    * Back substitution
    */
    for (int i = end - 2; i >= start; i--)
    {
        u[i] -= q[i + 1] * u[i + 1];
    }
    return u;
}