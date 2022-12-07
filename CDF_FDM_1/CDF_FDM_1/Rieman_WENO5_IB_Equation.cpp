#include"cfd.h"

/*
* -Using Riemann solvers to compute the flux at the interface
* -Using WENO-5 Scheme to reconstruct the left and right side fluxes at the interface
* -Using Runge-Kutta-3 Scheme for time integration
* -Data is saved every 0.025s
*/
using namespace std;

/*
* - Riemann solver: Rusanov
*/
vector<double> rusanov(int nx, vector<double>u, vector<double>uL, vector<double>uR, vector<double>f, vector<double>fL, vector<double>fR)
{
    vector<double> u_a(nx + 1);

    for (int i = 1; i < nx; i++)
    {
        u_a[i] = max(abs(u[i]), abs(u[i - 1]));
    }
    u_a[0] = max(abs(u[0]), abs(u[nx - 1]));
    u_a[nx] = max(abs(u[0]), abs(u[nx - 1]));

    /*Using Rusanov to compute Interface fluxes*/
    for (int i = 0; i < nx + 1; i++)
    {
        f[i] = 0.5 * (fR[i] + fL[i]) - 0.5 * u_a[i] * (uR[i] - uL[i]);
    }
    return f;
}
/*
* Calculate right hand term of the inviscid Burgers equation
* r = -udu/dx
*/
vector<double> rhs_rweno5(int nx, double dx, vector<double> u, vector<double> r)
{
    vector<double> f(nx + 1);            //flux computed at i+1/2,i-1/2
    vector<double> fL(nx + 1);           //flux component using the left reconstructed state
    vector<double> fR(nx + 1);           //flux component using the right reconstructed state
    vector<double> uL(nx + 1);          
    vector<double> uR(nx + 1);

    uL = fweno5L(nx, u, uL);
    uR = fweno5R(nx, u, uR);

    for (int i = 0; i < nx + 1; i++)
    {
        fL[i] = 0.5 * uL[i] * uL[i];
        fR[i] = 0.5 * uR[i] * uR[i];
    }

    f = rusanov(nx, u, uL, uR, f, fL, fR);

    for (int i = 0; i < nx; i++)
    {
        r[i] = -(f[i + 1] - f[i]) / dx;
    }

    return r;
}

vector< vector<double> > numerical_rweno5(int nx, int ns, int nt, double dx, double dt, vector<double> x, vector< vector<double> > u_n)
{
    vector<double> u_nn(nx);
    vector<double> u_nt(nx);
    vector<double> r(nx);

    int freq = ceil(nt / ns);

    for (int i = 0; i < nx; i++)
    {
        u_nn[i] = sin(2.0 * PI * x[i]);
        u_n[0][i] = u_nn[i];
    }

    for (int j = 1; j < nt + 1; j++)
    {
        r = rhs_rweno5(nx, dx, u_nn, r);
        for (int i = 0; i < nx; i++)
        {
            u_nt[i] = u_nn[i] + dt * r[i];
        }

        r = rhs_rweno5(nx, dx, u_nt, r);
        for (int i = 0; i < nx; i++)
        {
            u_nt[i] = 0.75 * u_nn[i] + 0.25 * u_nt[i] + 0.25 * dt * r[i];
        }

        r = rhs_rweno5(nx, dx, u_nt, r);
        for (int i = 0; i < nx; i++)
        {
            u_nn[i] = (1.0 / 3.0) * u_nn[i] + (2.0 / 3.0) * u_nt[i] + (2.0 / 3.0) * dt * r[i];
        }

        if (j % freq == 0)
        {
            for (int i = 0; i < nx; i++)
            {
                u_n[j / freq][i] = u_nn[i];
            }
        }

    }
    return u_n;
}

void Rieman_WENO5_IB_Equation()
{
    double x_l = 0.0;
    double x_r = 1.0;
    int nx = 200;
    double dx = (x_r - x_l) / nx;

    double t = 0.25;
    double dt = 0.0001;
    int nt = ceil(t / dt);

    int ns = 10;   //Save ten sets of data results
    double ds = t / ns;          //The time interval for saving data

    vector<double> x(nx);
    vector< vector<double> > u_n(ns + 1, vector<double>(nx));

    for (int i = 0; i < nx; i++)
    {
        x[i] = x_l + dx * i + 0.5 * dx;  //Assign node locations
    }

    u_n = numerical_rweno5(nx, ns, nt, dx, dt, x, u_n);

    ofstream field;
    field.open("field_final.csv", ios::out | ios::trunc);
    field << "x" << ",";
    for (int i = 0; i < ns + 1; i++)
    {
        field << i * ds << ",";
    }
    field << endl;

    for (int i = 0; i < nx; i++)
    {
        field << setprecision(16) << x[i] << ",";
        for (int j = 0; j < ns + 1; j++)
        {
            field << setprecision(16) << u_n[j][i] << ",";
        }
        field << endl;
    }

    field.close();

}