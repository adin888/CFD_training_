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
double compute_l2norm(int nx, std::vector<double> erro);
std::vector<double> rhs(int nx, double dx, double dt, std::vector<double> u, std::vector<double> r, double alpha);