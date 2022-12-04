#include"cfd.h"

/*
* -Using WENO-5 Scheme for spatial terms
* -Using Runge-Kutta-3 Scheme for time integration
* -Data is saved every 0.025s
*/
using namespace std;

/*
* nonlinear weights for upwind direction
*/
float wL(float u1, float u2, float u3, float u4, float u5)
{
    double eps = 1.0e-6;
    double u_L;

    /* smoothness indicators */
    double s1 = pow((13.0 / 12.0) * (u1 - 2.0 * u2 + u3), 2) + pow(0.25 * (u1 - 4.0 * u2 + 3.0 * u3), 2);
    double s2 = pow((13.0 / 12.0) * (u2 - 2.0 * u3 + u4), 2) + pow(0.25 * (u2 - u4), 2);
    double s3 = pow((13.0 / 12.0) * (u3 - 2.0 * u4 + u5), 2) + pow(0.25 * (3.0 * u3 - 4.0 * u4 + u5), 2);

    double d1 = 1.0 / 10.0;
    double d2 = 3.0 / 5.0;
    double d3 = 3.0 / 10.0;

    double c1 = d1 / pow(eps + s1, 2);
    double c2 = d2 / pow(eps + s2, 2);
    double c3 = d3 / pow(eps + s3, 2);

    double w1 = c1 / (c1 + c2 + c3);
    double w2 = c2 / (c1 + c2 + c3);
    double w3 = c3 / (c1 + c2 + c3);

    double b1 = u1 / 3.0 - 7.0 / 6.0 * u2 + 11.0 / 6.0 * u3;
    double b2 = -u2 / 6.0 + 5.0 / 6.0 * u3 + u4 / 3.0;
    double b3 = u3 / 3.0 + 5.0 / 6.0 * u4 - u5 / 6.0;

    return u_L = w1 * b1 + w2 * b2 + w3 * b3;
}

/*
* nonlinear weights for downwind direction
*/
float wR(float u1, float u2, float u3, float u4, float u5)
{
    double eps = 1.0e-6;
    double u_R;

    /* smoothness indicators */
    double s1 = pow((13.0 / 12.0) * (u1 - 2.0 * u2 + u3), 2) + pow(0.25 * (u1 - 4.0 * u2 + 3.0 * u3), 2);
    double s2 = pow((13.0 / 12.0) * (u2 - 2.0 * u3 + u4), 2) + pow(0.25 * (u2 - u4), 2);
    double s3 = pow((13.0 / 12.0) * (u3 - 2.0 * u4 + u5), 2) + pow(0.25 * (3.0 * u3 - 4.0 * u4 + u5), 2);

    double d1 = 3.0 / 10.0;
    double d2 = 3.0 / 5.0;
    double d3 = 1.0 / 10.0;

    double c1 = d1 / pow(eps + s1, 2);
    double c2 = d2 / pow(eps + s2, 2);
    double c3 = d3 / pow(eps + s3, 2);

    double w1 = c1 / (c1 + c2 + c3);
    double w2 = c2 / (c1 + c2 + c3);
    double w3 = c3 / (c1 + c2 + c3);

    double b1 = -u1 / 6.0 + 5.0 / 6.0 * u2 + u3 / 3.0;
    double b2 = u2 / 3.0 + 5.0 / 6.0 * u3 - u4 / 6.0;
    double b3 = 11.0 / 6.0 * u3 - 7.0 / 6.0 * u4 + u5 / 3.0;

    return u_R = w1 * b1 + w2 * b2 + w3 * b3;
}

/*
* WENO reconstruction for upwind methode
*/
vector<double> weno5L(int nx, vector<double>u, vector<double>uL)
{
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    /* Ghost points with Dirichlet condition*/
    int i = 0;
    u1 = 3.0 * u[i] - 2.0 * u[i + 1];          //periodic condition = u[nx-2]
    u2 = 2.0 * u[i] - u[i + 1];                //periodic condition = u[nx-1]
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[i + 2];
    uL[i] = wL(u1, u2, u3, u4, u5);

    i = 1;
    u1 = 2.0 * u[i - 1] - u[i];                //periodic condition = u[nx-1]
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[i + 2];
    uL[i] = wL(u1, u2, u3, u4, u5);

    for (int i = 2; i < nx - 1; i++)
    {
        u1 = u[i - 2];
        u2 = u[i - 1];
        u3 = u[i];
        u4 = u[i + 1];
        u5 = u[i + 2];
        uL[i] = wL(u1, u2, u3, u4, u5);
    }

    i = nx - 1;
    u1 = u[i - 2];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = 2.0 * u[i + 1] - u[i];                //periodic condition = u[1]
    uL[i] = wL(u1, u2, u3, u4, u5);
    return uL;
}

/*
* WENO reconstruction for downwind methode
*/
vector<double> weno5R(int nx, vector<double>u, vector<double>uR)
{
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    int i = 1;
    u1 = 2.0 * u[i - 1] - u[i];                //periodic condition = u[nx-1]
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[i + 2];
    uR[i] = wR(u1, u2, u3, u4, u5);

    for (int i = 2; i < nx - 1; i++)
    {
        u1 = u[i - 2];
        u2 = u[i - 1];
        u3 = u[i];
        u4 = u[i + 1];
        u5 = u[i + 2];
        uR[i] = wR(u1, u2, u3, u4, u5);
    }

    i = nx - 1;
    u1 = u[i - 2];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = 2 * u[i + 1] - u[i];                //periodic condition = u[1]
    uR[i] = wR(u1, u2, u3, u4, u5);

    i = nx;
    u1 = u[i - 2];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = 2.0 * u[i] - u[i - 1];                //periodic condition = u[1]
    u5 = 3.0 * u[i] - 2.0 * u[i - 1];                //periodic condition = u[2]
    uR[i] = wR(u1, u2, u3, u4, u5);
    return uR;
}

/*
* Calculate right hand term of the inviscid Burgers equation
* r = -udu/dx
*/
vector<double> rhs_weno5(int nx, double dx, vector<double> u, vector<double> r)
{
    vector<double> u_L(nx);
    vector<double>u_R(nx + 1);

    u_L = weno5L(nx, u, u_L);
    u_R = weno5R(nx, u, u_R);

    for (int i = 1; i < nx; i++)
    {
        if (u[i]>=0.0)
        {
            r[i] = -u[i] * (u_L[i] - u_L[i - 1]) / dx;
        }
        else
        {
            r[i] = -u[i] * (u_R[i+1] - u_R[i]) / dx;
        }
    }
    /* periodic condition
    if (u[0]>=0.0)
        {
            r[0] = -u[0] * (u_L[0] - u_L[nx - 1]) / dx;
        }
        else
        {
            r[0] = -u[0] * (u_R[1] - u_R[nx]) / dx;
        }
    */
    return r;
}

vector< vector<double> > numerical_weno5(int nx, int ns, int nt, double dx, double dt, vector<double> x, vector< vector<double> > u_n)
{
    vector<double> u_nn(nx + 1);
    vector<double> u_nt(nx + 1);
    vector<double> r(nx);

    int freq = ceil(nt / ns);

    for (int i = 0; i < nx + 1; i++)
    {
        u_nn[i] = sin(2.0 * PI * x[i]);
        u_n[0][i] = u_nn[i];
    }

    /* Dirichlet boundary condition*/
    u_nn[0], u_nn[nx] = 0.0, 0.0;
    u_nt[0], u_nt[nx] = 0.0, 0.0;
    u_n[0][0], u_n[0][nx] = 0.0, 0.0;

    for (int j = 1; j < nt + 1; j++)
    {
        r = rhs_weno5(nx, dx, u_nn, r);
        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = u_nn[i] + dt * r[i];
        }

        //u_nt[nx]=u_nt[0];   //periodic condition

        r = rhs_weno5(nx, dx, u_nt, r);
        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = 0.75 * u_nn[i] + 0.25 * u_nt[i] + 0.25 * dt * r[i];
        }

        //u_nt[nx]=u_nt[0];   //periodic condition

        r = rhs_weno5(nx, dx, u_nt, r);
        for (int i = 1; i < nx; i++)
        {
            u_nn[i] = (1.0 / 3.0) * u_nn[i] + (2.0 / 3.0) * u_nt[i] + (2.0 / 3.0) * dt * r[i];
        }

        //u_nt[nx]=u_nt[0];   //periodic condition

        if (j%freq == 0)
        {
            for (int i = 0; i < nx + 1; i++)
            {
                u_n[j/freq][i] = u_nn[i];
            }
        }

    }
    return u_n;
}

void WENO5_IB_Equation()
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

    vector<double> x(nx + 1);
    vector< vector<double> > u_n(ns + 1, vector<double>(nx + 1));

    for (int i = 0; i < nx + 1; i++)
    {
        x[i] = x_l + dx * i;  //Assign node locations
    }

    u_n = numerical_weno5(nx, ns, nt, dx, dt, x, u_n);

    ofstream field;
    field.open("field_final.csv", ios::out | ios::trunc);
    field << "x" << ",";
    for (int i = 0; i < ns + 1; i++)
    {
        field << i * ds << ",";
    }
    field << endl;

    for (int i = 0; i < nx + 1; i++)
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