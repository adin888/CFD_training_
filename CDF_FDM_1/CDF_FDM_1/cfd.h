#pragma once
#include <iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<algorithm>
#include<vector>

#define PI acos(-1)

void FTCS_Heat_Equation();
void Runge_Heat_Equation();
void CN_Heat_Equation();
void ICP_Heat_Equation();
void WENO5_IB_Equation();
void CRWENO5_IB_Equation();
void Flux_WENO5_IB_Equation();
void Rieman_WENO5_IB_Equation();
double compute_l2norm(int nx, std::vector<double> erro);
std::vector<double> tdms(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, std::vector<double> u, int start, int end);
std::vector<double> ctdms(std::vector<double> a, std::vector<double> b, std::vector<double> c, double alpha, double beta,
    std::vector<double> d, std::vector<double> u, int start, int end);
float wL(float u1, float u2, float u3, float u4, float u5);
float wR(float u1, float u2, float u3, float u4, float u5);
std::vector<double> fweno5L(int nx, std::vector<double>fP, std::vector<double>fL);
std::vector<double> fweno5R(int nx, std::vector<double>fN, std::vector<double>fR);