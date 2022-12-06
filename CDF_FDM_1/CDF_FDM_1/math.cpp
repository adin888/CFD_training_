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
* -Solution to tridigonal system using cyclic Thomas algorithm(with Sherman-Morrison formula)
*/
vector<double> ctdms(vector<double> a, vector<double> b, vector<double> c, double alpha, double beta,
                    vector<double> d, vector<double> u,int start, int end)
{
    vector<double> bb(end + 1);
    vector<double> p(end + 1, 0.0);
    vector<double> y(end + 1);
    vector<double> z(end + 1);

    double gamma = -b[start];
    bb[start] = b[start] - gamma;
    bb[end] = b[end] - alpha * beta / gamma;

    for (int i = start+1; i < end; i++)
    {
        bb[i] = b[i];
    }

    y = tdms(a, bb, c, d, y, start, end);

    p[start] = gamma;
    p[end] = alpha;

    z = tdms(a, bb, c, p, z, start, end);

    double f = (y[start] + beta * y[end] / gamma) / (1.0 + z[start] + beta * z[end] / gamma);

    for (int i = start; i <= end; i++)
    {
        u[i] = y[i] - f * z[i];
    }
    return u;
}
/*
* -Tridiagonal matrix algorithm(Thomas algorithm) Au=d, A is a matrix containing a_i,b_i,c_i
*/
vector<double> tdms(vector<double> a, vector<double> b, vector<double> c, vector<double> d, vector<double> u, int start, int end)
{
    vector<double> q(end + 1);                    //Storage of superdiagonal array
    double b_tem = b[start];                  //Temporary storage of b
    u[start] = d[start] / b_tem;
    /*
    * Forward elimination
    */
    for (int i = start + 1; i <= end; i++)
    {
        q[i] = c[i - 1] / b_tem;
        b_tem = b[i] - a[i] * q[i];
        u[i] = (d[i] - a[i] * u[i - 1]) / b_tem;
    }
    /*
    * Back substitution
    */
    for (int i = end - 1; i >= start; i--)
    {
        u[i] -= q[i + 1] * u[i + 1];
    }
    return u;
}