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
double compute_l2norm(int nx, std::vector<double> erro);
std::vector<double> rhs(int nx, double dx, double dt, std::vector<double> u, std::vector<double> r, double alpha);
std::vector<double> tdms(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, int start, int end);