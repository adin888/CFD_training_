#include"cfd.h"

/*With Forward Time Central Space Scheme*/

#define PI acos(-1)
using namespace std;

float compute_l2norm(int nx, vector<float> erro)
{
    float rms = 0.0;
    for (int i = 1; i < nx; i++)
    {
        rms = rms + pow(erro[i], 2);
    }
    rms = sqrt(rms / (nx - 1));
    return rms;
}
void FTCS_Heat_Equation_1()
{
    float x_l = -1.0;
    float x_r = 1.0;
    float dx = 0.025;
    int nx = ceil((x_r - x_l) / dx);

    float t = 1.0;
    float dt = 0.0025;
    int nt = ceil(t / dt);

    float alpha = 1 / (PI * PI);

    vector<float> x(nx + 1);
    vector<float> u_a(nx + 1);
    vector< vector<float> > u_n(nt + 1, vector<float>(nx + 1));
    vector<float> u_erro(nx + 1);
    /*float* x = new float[nx + 1];
    float *u_a = new float[nx + 1];     //Analytical solution
    float *u_n = new float[nt + 1, nx + 1];  //Matrix for storing numerical solution at every time step
    float *u_erro = new float[nx + 1];*/


    for (int i = 0; i < nx + 1; i++)
    {
        x[i] = x_l + dx * i;  //Assign node locations
        u_n[0][i] = -sin(PI * x[i]);   //initial condition for u_n
        u_a[i] = -exp(-t) * sin(PI * x[i]);  //initial condition for u_a
    }

    u_n[0][0] = 0.0;
    u_n[0][nx] = 0.0;

    float beta = alpha * dt / (dx * dx);   

    for (int i = 1; i < nt + 1; i++)
    {
        for (int j = 1; j < nx; j++)
        {
            u_n[i][j] = u_n[i - 1][j] + beta * (u_n[i - 1][j + 1] - 2.0 * u_n[i - 1][j] + u_n[i - 1][j - 1]);
        }
        u_n[i][0] = 0.0;
        u_n[i][nx] = 0.0;
    }

    for (int i = 0; i < nx + 1; i++)
    {
        u_erro[i] = u_n[nt][i] - u_a[i];
    }
    float rms_error = compute_l2norm(nx, u_erro);
    vector<float> u_erro_abs(nx + 1);
    for (int i = 0; i < nx + 1; i++)
    {
        u_erro_abs[i] = abs(u_erro[i]);
    }
    //float u_erro_abs = fabs(*u_erro);
    float max_erro = *max_element(u_erro_abs.begin(), u_erro_abs.end());
    //cout << max_element(u_erro_abs, u_erro_abs + nx) << endl;

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