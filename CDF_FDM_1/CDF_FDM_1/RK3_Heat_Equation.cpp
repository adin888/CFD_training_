#include"cfd.h"

/*
With Runge-Kutta Numerical Scheme
*/

using namespace std;
vector< vector<double> > numerical(int nx, int nt, double dx, double dt, vector<double> x, vector< vector<double> > u_n, double alpha)
{
    vector<double> u_nn(nx + 1);
    vector<double> u_nt(nx + 1);
    vector<double> r(nx);

    for (int i = 0; i < nx + 1; i++)
    {
        u_nn[i] = -sin(PI * x[i]);           //initial condition for u_n and u_nn
        u_n[0][i] = u_nn[i];
    }

    u_nn[0] = 0.0;
    u_nn[nx] = 0.0;
    u_nt[0] = 0.0;
    u_nt[nx] = 0.0;

    for (int j = 1; j < nt + 1; j++)
    {
        r = rhs(nx, dx, dt, u_nn, r, alpha);
        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = u_nn[i] + r[i];
        }

        r = rhs(nx, dx, dt, u_nt, r, alpha);
        for (int i = 1; i < nx; i++)
        {
            u_nt[i] = 0.75 * u_nn[i] + 0.25 * u_nt[i] + 0.25 * r[i];
        }

        r = rhs(nx, dx, dt, u_nt, r, alpha);
        for (int i = 1; i < nx; i++)
        {
            u_nn[i] = (1.0 / 3.0) * u_nn[i] + (2.0 / 3.0) * u_nt[i] + (2.0 / 3.0) * r[i];
        }

        for (int i = 0; i < nx + 1; i++)
        {
            u_n[j][i] = u_nn[i];
        }
    }
    return u_n;
}
void Runge_Heat_Equation()
{
    double x_l = -1.0;
    double x_r = 1.0;
    double dx = 0.025;
    int nx = ceil((x_r - x_l) / dx);

    double t = 1.0;
    double dt = 0.0025;
    int nt = ceil(t / dt);

    double alpha = 1 / (PI * PI);

    vector<double> x(nx + 1);
    vector<double> u_a(nx + 1);
    vector< vector<double> > u_n(nt + 1, vector<double>(nx + 1));
    vector<double> u_erro(nx + 1, 0);

    for (int i = 0; i < nx + 1; i++)
    {
        x[i] = x_l + dx * i;  //Assign node locations
        u_a[i] = -exp(-t) * sin(PI * x[i]);  ////Analytical solution in the last step
    }

    u_n = numerical(nx, nt, dx, dt, x, u_n, alpha);

    for (int i = 0; i < nx + 1; i++)
    {
        u_erro[i] = u_n[nt][i] - u_a[i];
    }

    double rms_error = compute_l2norm(nx, u_erro);

    vector<double> u_erro_abs(nx + 1);
    for (int i = 0; i < nx + 1; i++)
    {
        u_erro_abs[i] = abs(u_erro[i]);
    }

    double max_erro = *max_element(u_erro_abs.begin(), u_erro_abs.end());

    ofstream outfile;
    outfile.open("result.txt", ios::out | ios::trunc);
    outfile << "Erro details" << endl;
    outfile << "L-2 Norm = " << rms_error << endl;
    outfile << "Maximum Norm = " << max_erro << endl;

    ofstream field;
    field.open("field_final.csv", ios::out | ios::trunc);
    field << "x" << "," << "u_a" << "," << "u_n" << "," << "u_erro" << endl;
    for (int i = 0; i < nx + 1; i++)
    {
        field << setprecision(16) << x[i] << "," << u_a[i] << "," << u_n[nt][i] << "," << u_erro_abs[i] << endl;
    }
    outfile.close();
    field.close();

}