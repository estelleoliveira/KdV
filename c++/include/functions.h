#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>

std::vector<double> linspace(double start, double end, int num);
double max_abs(const std::vector<std::vector<double>>& vec);
void periodicite(int j, int M, int& jm, int& jp, int& jv, int& jw);
double compute_residu(int M, int N, const std::vector<std::vector<double>>& u_num, const std::vector<std::vector<double>>& u_exact);

#endif // FUNCTIONS_H
