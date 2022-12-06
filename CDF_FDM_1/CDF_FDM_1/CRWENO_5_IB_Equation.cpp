#include"cfd.h"

/*
* -Using CRWENO-5 Scheme for spatial terms
* -Using Runge-Kutta-3 Scheme for time integration
* -Data is saved every 0.025s
* -periodic condition
*/
using namespace std;

/*
* nonlinear weights for upwind direction
*/
vector<double> crwL(float u1, float u2, float u3, float u4, float u5, vector<double> coeff)
{
    double eps = 1.0e-6;

    /* smoothness indicators */
    double s1 = pow((13.0 / 12.0) * (u1 - 2.0 * u2 + u3), 2) + pow(0.25 * (u1 - 4.0 * u2 + 3.0 * u3), 2);
    double s2 = pow((13.0 / 12.0) * (u2 - 2.0 * u3 + u4), 2) + pow(0.25 * (u2 - u4), 2);
    double s3 = pow((13.0 / 12.0) * (u3 - 2.0 * u4 + u5), 2) + pow(0.25 * (3.0 * u3 - 4.0 * u4 + u5), 2);

    double d1 = 1.0 / 5.0;
    double d2 = 1.0 / 2.0;
    double d3 = 3.0 / 10.0;

    double c1 = d1 / pow(eps + s1, 2);
    double c2 = d2 / pow(eps + s2, 2);
    double c3 = d3 / pow(eps + s3, 2);

    double w1 = c1 / (c1 + c2 + c3);
    double w2 = c2 / (c1 + c2 + c3);
    double w3 = c3 / (c1 + c2 + c3);

    coeff[0] = (2.0 * w1 + w2) / 3.0;
    coeff[1] = (w1 + 2.0 * w2 + 2.0 * w3) / 3.0;
    coeff[2] = w3 / 3.0;

    coeff[3] = w1 / 6.0;
    coeff[4] = (5.0 * w1 + 5.0 * w2 + w3) / 6.0;
    coeff[5] = (w2 + 5.0 * w3) / 6.0;

    return coeff;
}

/*
* nonlinear weights for downwind direction
*/
vector<double> crwR(float u1, float u2, float u3, float u4, float u5, vector<double> coeff)
{
    double eps = 1.0e-6;

    /* smoothness indicators */
    double s1 = pow((13.0 / 12.0) * (u1 - 2.0 * u2 + u3), 2) + pow(0.25 * (u1 - 4.0 * u2 + 3.0 * u3), 2);
    double s2 = pow((13.0 / 12.0) * (u2 - 2.0 * u3 + u4), 2) + pow(0.25 * (u2 - u4), 2);
    double s3 = pow((13.0 / 12.0) * (u3 - 2.0 * u4 + u5), 2) + pow(0.25 * (3.0 * u3 - 4.0 * u4 + u5), 2);

    double d1 = 3.0 / 10.0;
    double d2 = 1.0 / 2.0;
    double d3 = 1.0 / 5.0;

    double c1 = d1 / pow(eps + s1, 2);
    double c2 = d2 / pow(eps + s2, 2);
    double c3 = d3 / pow(eps + s3, 2);

    double w1 = c1 / (c1 + c2 + c3);
    double w2 = c2 / (c1 + c2 + c3);
    double w3 = c3 / (c1 + c2 + c3);

    /**/
    coeff[0] = w1 / 3.0;
    coeff[1] = (w3 + 2.0 * w2 + 2.0 * w1) / 3.0;
    coeff[2] = (2.0 * w3 + w2) / 3.0;

    coeff[3] = (w2 + 5.0 * w1) / 6.0;
    coeff[4] = (5.0 * w3 + 5.0 * w2 + w1) / 6.0;
    coeff[5] = w3 / 6.0;
    

    /*
    coeff[0] = w3 / 3.0;
    coeff[1] = (w1 + 2.0 * w2 + 2.0 * w3) / 3.0;
    coeff[2] = (2.0 * w1 + w2) / 3.0;

    coeff[3] = w1 / 6.0;
    coeff[4] = (5.0 * w1 + 5.0 * w2 + w3) / 6.0;
    coeff[5] = (w2 + 5.0 * w3) / 6.0;
    */
    return coeff;
}

/*
* WENO reconstruction for upwind methode
*/
vector<double> crweno5L(int nx, vector<double>u, vector<double>uL)
{
    vector<double> a(nx);                    //subdiagonal array of tridiagonal matrix
    vector<double> b(nx);                    //diagonal array of tridiagonal matrix
    vector<double> c(nx);                    //superdiagonal array of tridiagonal matrix
    vector<double> r(nx);
    vector<double> coeff(6);
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    /* Ghost points with periodic condition*/
    int i = 0;
    u1 = u[nx - 2];
    u2 = u[nx - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[i + 2];

    coeff = crwL(u1, u2, u3, u4, u5, coeff);
    a[i] = coeff[0];
    b[i] = coeff[1];
    c[i] = coeff[2];
    r[i] = coeff[3] * u[nx - 1] + coeff[4] * u[i] + coeff[5] * u[i + 1];

    i = 1;
    u1 = u[nx - 1];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[i + 2];

    coeff = crwL(u1, u2, u3, u4, u5, coeff);
    a[i] = coeff[0];
    b[i] = coeff[1];
    c[i] = coeff[2];
    r[i] = coeff[3] * u[i - 1] + coeff[4] * u[i] + coeff[5] * u[i + 1];

    for (int i = 2; i < nx - 1; i++)
    {
        u1 = u[i - 2];
        u2 = u[i - 1];
        u3 = u[i];
        u4 = u[i + 1];
        u5 = u[i + 2];
        
        coeff = crwL(u1, u2, u3, u4, u5, coeff);
        a[i] = coeff[0];
        b[i] = coeff[1];
        c[i] = coeff[2];
        r[i] = coeff[3] * u[i - 1] + coeff[4] * u[i] + coeff[5] * u[i + 1];
    }

    i = nx - 1;
    u1 = u[i - 2];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[1];
    
    coeff = crwL(u1, u2, u3, u4, u5, coeff);
    a[i] = coeff[0];
    b[i] = coeff[1];
    c[i] = coeff[2];
    r[i] = coeff[3] * u[i - 1] + coeff[4] * u[i] + coeff[5] * u[i + 1];

    double alpha = c[nx - 1];
    double beta = a[0];

    uL = ctdms(a, b, c, alpha, beta, r, uL, 0, nx - 1);

    return uL;
}

/*
* WENO reconstruction for downwind methode
*/
vector<double> crweno5R(int nx, vector<double>u, vector<double>uR)
{
    vector<double> a(nx + 1);                    //subdiagonal array of tridiagonal matrix
    vector<double> b(nx + 1);                    //diagonal array of tridiagonal matrix
    vector<double> c(nx + 1);                    //superdiagonal array of tridiagonal matrix
    vector<double> r(nx + 1);
    vector<double> coeff(6);
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    int i = 1;
    u1 = u[nx - 1];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[i + 2];

    coeff = crwL(u1, u2, u3, u4, u5, coeff);
    a[i] = coeff[0];
    b[i] = coeff[1];
    c[i] = coeff[2];
    r[i] = coeff[3] * u[i - 1] + coeff[4] * u[i] + coeff[5] * u[i + 1];

    for (int i = 2; i < nx - 1; i++)
    {
        u1 = u[i - 2];
        u2 = u[i - 1];
        u3 = u[i];
        u4 = u[i + 1];
        u5 = u[i + 2];
        
        coeff = crwL(u1, u2, u3, u4, u5, coeff);
        a[i] = coeff[0];
        b[i] = coeff[1];
        c[i] = coeff[2];
        r[i] = coeff[3] * u[i - 1] + coeff[4] * u[i] + coeff[5] * u[i + 1];
    }

    i = nx - 1;
    u1 = u[i - 2];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[i + 1];
    u5 = u[1];
    
    coeff = crwL(u1, u2, u3, u4, u5, coeff);
    a[i] = coeff[0];
    b[i] = coeff[1];
    c[i] = coeff[2];
    r[i] = coeff[3] * u[i - 1] + coeff[4] * u[i] + coeff[5] * u[i + 1];

    i = nx;
    u1 = u[i - 2];
    u2 = u[i - 1];
    u3 = u[i];
    u4 = u[1];
    u5 = u[2];
    
    coeff = crwL(u1, u2, u3, u4, u5, coeff);
    a[i] = coeff[0];
    b[i] = coeff[1];
    c[i] = coeff[2];
    r[i] = coeff[3] * u[i - 1] + coeff[4] * u[i] + coeff[5] * u[1];

    double alpha = c[nx];
    double beta = a[1];

    uR = ctdms(a, b, c, alpha, beta, r, uR, 1, nx);

    return uR;
}

/*
* Calculate right hand term of the inviscid Burgers equation
* r = -udu/dx
*/
vector<double> rhs_crweno5(int nx, double dx, vector<double> u, vector<double> r)
{
    vector<double> u_L(nx);
    vector<double> u_R(nx + 1);

    u_L = crweno5L(nx, u, u_L);
    u_R = crweno5R(nx, u, u_R);

    for (int i = 1; i < nx; i++)
    {
        if (u[i] >= 0.0)
        {
            r[i] = -u[i] * (u_L[i] - u_L[i - 1]) / dx;
        }
        else
        {
            r[i] = -u[i] * (u_R[i + 1] - u_R[i]) / dx;
        }
    }
    /* 
    * - periodic condition
    */
    if (u[0]>=0.0)
        {
            r[0] = -u[0] * (u_L[0] - u_L[nx - 1]) / dx;
        }
    else
        {
            r[0] = -u[0] * (u_R[1] - u_R[nx]) / dx;
        }
    return r;
}

vector< vector<double> > numerical_crweno5(int nx, int ns, int nt, double dx, double dt, vector<double> x, vector< vector<double> > u_n)
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

    for (int j = 1; j < nt + 1; j++)
    {
        r = rhs_crweno5(nx, dx, u_nn, r);

        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = u_nn[i] + dt * r[i];
        }

        u_nt[nx] = u_nt[0];   //periodic condition

        r = rhs_crweno5(nx, dx, u_nt, r);

        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = 0.75 * u_nn[i] + 0.25 * u_nt[i] + 0.25 * dt * r[i];
        }

        u_nt[nx] = u_nt[0];   //periodic condition

        r = rhs_crweno5(nx, dx, u_nt, r);
        for (int i = 1; i < nx; i++)
        {
            u_nn[i] = (1.0 / 3.0) * u_nn[i] + (2.0 / 3.0) * u_nt[i] + (2.0 / 3.0) * dt * r[i];
        }

        u_nn[nx] = u_nn[0];   //periodic condition

        if (j % freq == 0)
        {
            for (int i = 0; i < nx + 1; i++)
            {
                u_n[j / freq][i] = u_nn[i];
            }
        }

    }
    return u_n;
}

void CRWENO5_IB_Equation()
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

    u_n = numerical_crweno5(nx, ns, nt, dx, dt, x, u_n);

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