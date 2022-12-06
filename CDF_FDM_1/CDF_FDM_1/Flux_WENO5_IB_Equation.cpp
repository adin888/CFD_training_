#include"cfd.h"

/*
* -Using Lax-Friedrichs flux splitting to compute positive and negative components of the flux
* -Using WENO-5 Scheme to reconstruct the left and right side fluxes at the interface
* -Using Runge-Kutta-3 Scheme for time integration
* -Data is saved every 0.025s
*/
using namespace std;

/*
* WENO reconstruction for upwind methode
*/
vector<double> fweno5L(int nx, vector<double>fP, vector<double>fL)
{
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    /* Ghost points with periodic condition*/
    int i = -1;
    u1 = fP[nx - 3];
    u2 = fP[nx - 2];
    u3 = fP[nx - 1];
    u4 = fP[i + 1];
    u5 = fP[i + 2];
    fL[i + 1] = wL(u1, u2, u3, u4, u5);

    i = 0;
    u1 = fP[nx - 2];
    u2 = fP[nx - 1];
    u3 = fP[i];
    u4 = fP[i + 1];
    u5 = fP[i + 2];
    fL[i + 1] = wL(u1, u2, u3, u4, u5);

    i = 1;
    u1 = fP[nx - 1];
    u2 = fP[i - 1];
    u3 = fP[i];
    u4 = fP[i + 1];
    u5 = fP[i + 2];
    fL[i + 1] = wL(u1, u2, u3, u4, u5);

    for (int i = 2; i < nx - 2; i++)
    {
        u1 = fP[i - 2];
        u2 = fP[i - 1];
        u3 = fP[i];
        u4 = fP[i + 1];
        u5 = fP[i + 2];
        fL[i + 1] = wL(u1, u2, u3, u4, u5);
    }

    i = nx - 2;
    u1 = fP[i - 2];
    u2 = fP[i - 1];
    u3 = fP[i];
    u4 = fP[i + 1];
    u5 = fP[0];
    fL[i + 1] = wL(u1, u2, u3, u4, u5);

    i = nx - 1;
    u1 = fP[i - 2];
    u2 = fP[i - 1];
    u3 = fP[i];
    u4 = fP[0];
    u5 = fP[1];
    fL[i + 1] = wL(u1, u2, u3, u4, u5);
    return fL;
}

/*
* WENO reconstruction for downwind methode
*/
vector<double> fweno5R(int nx, vector<double>fN, vector<double>fR)
{
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    int i = 0;
    u1 = fN[nx - 2];
    u2 = fN[nx - 1];
    u3 = fN[i];
    u4 = fN[i + 1];
    u5 = fN[i + 2];
    fR[i] = wR(u1, u2, u3, u4, u5);

    i = 1;
    u1 = fN[nx - 1];
    u2 = fN[i - 1];
    u3 = fN[i];
    u4 = fN[i + 1];
    u5 = fN[i + 2];
    fR[i] = wR(u1, u2, u3, u4, u5);

    for (int i = 2; i < nx - 2; i++)
    {
        u1 = fN[i - 2];
        u2 = fN[i - 1];
        u3 = fN[i];
        u4 = fN[i + 1];
        u5 = fN[i + 2];
        fR[i] = wR(u1, u2, u3, u4, u5);
    }

    i = nx - 2;
    u1 = fN[i - 2];
    u2 = fN[i - 1];
    u3 = fN[i];
    u4 = fN[i + 1];
    u5 = fN[0];
    fR[i]= wR(u1, u2, u3, u4, u5);
         
    i = nx - 1;
    u1 = fN[i - 2];
    u2 = fN[i - 1];
    u3 = fN[i];
    u4 = fN[0];
    u5 = fN[1];
    fR[i] = wR(u1, u2, u3, u4, u5);

    i = nx;
    u1 = fN[i - 2];
    u2 = fN[i - 1];
    u3 = fN[0];
    u4 = fN[1];
    u5 = fN[2];
    fR[i] = wR(u1, u2, u3, u4, u5);
    return fR;
}

vector<double> wavespeed(int nx, vector<double>u, vector<double>u_a)
{
    for (int i = 2; i < nx - 2; i++)
    {
        u_a[i] = max({ abs(u[i - 2]), abs(u[i - 1]), abs(u[i]), abs(u[i + 1]), abs(u[i + 2]) });
    }
    /*periodic condition*/
    u_a[0] = max({ abs(u[nx - 2]), abs(u[nx - 1]), abs(u[0]), abs(u[1]), abs(u[2]) });
    u_a[1] = max({ abs(u[nx - 1]), abs(u[0]), abs(u[1]), abs(u[2]), abs(u[3]) });
    u_a[nx - 2] = max({ abs(u[nx - 4]), abs(u[nx - 3]), abs(u[nx - 2]), abs(u[nx - 1]), abs(u[0]) });
    u_a[nx - 1] = max({ abs(u[nx - 3]), abs(u[nx - 2]), abs(u[nx - 1]), abs(u[0]), abs(u[1]) });
    
    return u_a;
}
/*
* Calculate right hand term of the inviscid Burgers equation
* r = -udu/dx
*/
vector<double> rhs_fweno5(int nx, double dx, vector<double> u, vector<double> r)
{
    vector<double> f(nx);            //flux computed at nodal points and positive and negative splitting
    vector<double> fP(nx);
    vector<double> fN(nx);
    vector<double> u_a(nx);          //wave speed at nodal points

    vector<double> fL(nx+1);         //left and right side fluxes at the interface
    vector<double> fR(nx + 1);

    for (int i = 0; i < nx; i++)
    {
        f[i] = 0.5 * u[i] * u[i];
    }

    u_a = wavespeed(nx, u, u_a);

    for (int i = 0; i < nx; i++)
    {
        fP[i] = 0.5 * (f[i] + u_a[i] * u[i]);
        fN[i] = 0.5 * (f[i] - u_a[i] * u[i]);
    }

    fL = fweno5L(nx, fP, fL);
    fR = fweno5R(nx, fN, fR);

    for (int i = 1; i < nx; i++)
    {
        r[i] = -(fL[i + 1] - fL[i]) / dx - (fR[i + 1] - fR[i]) / dx;
    }

    return r;
}

vector< vector<double> > numerical_fweno5(int nx, int ns, int nt, double dx, double dt, vector<double> x, vector< vector<double> > u_n)
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
        r = rhs_fweno5(nx, dx, u_nn, r);
        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = u_nn[i] + dt * r[i];
        }

        r = rhs_fweno5(nx, dx, u_nt, r);
        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = 0.75 * u_nn[i] + 0.25 * u_nt[i] + 0.25 * dt * r[i];
        }

        r = rhs_fweno5(nx, dx, u_nt, r);
        for (int i = 1; i < nx; i++)
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

void Flux_WENO5_IB_Equation()
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
        x[i] = x_l + dx * i;  //Assign node locations
    }

    u_n = numerical_fweno5(nx, ns, nt, dx, dt, x, u_n);

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
