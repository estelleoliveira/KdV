#include <iostream>
#include <vector>
#include <cmath> //pow et sqrt

#include "../include/fex.h"
#include "../include/functions.h"

std::vector<std::vector<double>> fex(int M, int N, double T, double L, double delta, double U0, double Uinf, double x0) {
    double DELTA = delta / std::sqrt((U0 - Uinf) / 12);
    double c = Uinf + (U0 - Uinf) / 3;
    double dx = L / M;  //spatial step
    auto x = linspace(0, L - dx, M);
    auto t = linspace(0, T, N + 1);

    std::vector<std::vector<double>> fex(M, std::vector<double>(N + 1, 0.0));  // Solution exacte

    for (int n = 0; n < N + 1; n++) {
        for (int j = 0; j < M; j++) {
            double T1 = t[n];
            double X = x[j];
            double arg = (X - x0 - c * T1 - std::round((X - x0 - c * T1) / L) * L) / DELTA;
            fex[j][n] = U0 / std::pow(std::cosh(arg), 2);
        }
    }
    
    return fex;
}