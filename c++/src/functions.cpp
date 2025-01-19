#include <cmath>
#include <vector>

#include "../include/functions.h"

//creating the linspace function
std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + step * i;
    }
    return result;
}

double max_abs(const std::vector<std::vector<double>>& vec) {
    double max_val = 0.0;
    for (const auto& row : vec) {
        for (double val : row) {
            max_val = std::max(max_val, std::abs(val));
        }
    }
    return max_val;
}

//periodic boundary conditions
void periodicite(int j, int M, int& jm, int& jp, int& jv, int& jw) {
    jm = (j - 1 + M) % M;   //using modulo 
    jp = (j + 1) % M;
    jv = (j - 2 + M) % M;
    jw = (j + 2) % M;
}

//compute residu
double compute_residu(int M, int N, const std::vector<std::vector<double>>& u_num, const std::vector<std::vector<double>>& u_exact) {
    double residu = 0.0;
    
    for (int n = 0; n < N; n++) {
        for (int j = 0; j < M; j++) {
            residu += pow(u_num[j][n] - u_exact[j][n], 2);
        }
    }
    
    return std::sqrt(residu);   //L2 norm
}