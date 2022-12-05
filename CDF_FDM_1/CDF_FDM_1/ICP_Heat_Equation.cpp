#include"cfd.h"

/*
* With Implicit Compact Pade (ICP) Scheme and Thomas algorithm for solving the tridiagonal matrix.
*/

using namespace std;

void ICP_Heat_Equation()
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
    vector<double> u_tem(nx + 1, 0.0);
    vector<double> u_erro(nx + 1, 0.0);
    vector<double> a(nx + 1, 0.0);
    vector<double> b(nx + 1, 0.0);
    vector<double> c(nx + 1, 0.0);
    vector<double> d(nx + 1, 0.0);

    for (int i = 0; i < nx + 1; i++)
    {
        x[i] = x_l + dx * i;  // Assign node locations
        u_n[0][i] = -sin(PI * x[i]);
        u_a[i] = -exp(-t) * sin(PI * x[i]);  // Analytical solution in the last step
    }

    /* boundary conditions */
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = 0.0;
    d[0] = 0.0;
    u_n[0][0] = 0.0;

    a[nx] = 0.0;
    b[nx] = 1.0;
    c[nx] = 0.0;
    d[nx] = 0.0;
    u_n[0][nx] = 0.0;

    for (int i = 1; i < nx; i++)
    {
        a[i] = 12.0 / (dx * dx) - 2.0 / (alpha * dt);
        b[i] = -24.0 / (dx * dx) - 20.0 / (alpha * dt);
        c[i] = 12.0 / (dx * dx) - 2.0 / (alpha * dt);
    }

    int start = 0;
    int end = nx;

    for (int j = 1; j < nt + 1; j++)
    {

        for (int i = 1; i < nx; i++)
        {
            d[i] = -2.0 / (alpha * dt) * (u_n[j - 1][i + 1] + 10.0 * u_n[j - 1][i] + u_n[j - 1][i - 1])
                - 12.0 / (dx * dx) * (u_n[j - 1][i + 1] - 2.0 * u_n[j - 1][i] + u_n[j - 1][i - 1]);
        }

        u_tem = tdms(a, b, c, d, u_tem, start, end);

        for (int i = 0; i < nx + 1; i++)
        {
            u_n[j][i] = u_tem[i];
        }
    }

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