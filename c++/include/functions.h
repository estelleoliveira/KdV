#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>

std::vector<double> linspace(double start, double end, int num);
double max_abs(const std::vector<std::vector<double>>& vec);
void periodicite(int j, int M, int& jm, int& jp, int& jv, int& jw);
double compute_residu(int M, int N, const std::vector<std::vector<double>>& u_num, const std::vector<std::vector<double>>& u_exact);

void plotResults(const std::string& wang_file, const std::string& zk_file, const std::string& sixp_file, const std::string& exact_file);
void saveToFile(const std::vector<std::vector<double>>& u, const std::vector<double>& x, const std::string& filename);

#endif // FUNCTIONS_H
