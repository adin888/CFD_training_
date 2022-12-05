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
double compute_l2norm(int nx, std::vector<double> erro);
std::vector<double> tdms(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, std::vector<double> u, int start, int end);
std::vector<double> ctdms(std::vector<double> a, std::vector<double> b, std::vector<double> c, double alpha, double beta,
    std::vector<double> d, std::vector<double> u, int start, int end);