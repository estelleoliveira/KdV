#include <iostream>
#include <vector>
#include <cmath> //pow et sqrt

#include "../include/sixP.h"
#include "../include/functions.h"
#include "../include/fex.h"


// KdV 6P solver
std::vector<std::vector<double>> sixP_schema(double L, double T, int M, int N, double eta, double mu, double x0, double Uinf, double U0) {
    std::cout << "Calculating sixP scheme..." << std::endl;

    double dx = L / M;
    double dt = T / N;
    
    auto x = linspace(0, L - dx, M);
    auto t = linspace(0, T, N + 1);
    
    std::vector<std::vector<double>> u(M, std::vector<double>(N + 1, 0.0));
    
    double alpha = (eta * dt) / (3 * dx);
    double S = (mu * mu * dt) / (pow(dx, 3));
    double delta = 0.022;
    
    //print stability parameters
    double rnum = dt/dx;
    double rtheoric = 1 / (max_abs(fex(M, N, T, L, delta, U0, Uinf, x0)) + (4 * (pow(mu, 2))) / (pow(dx, 2)));
    std::cout << "rZK numérique: " << rnum << std::endl;
    std::cout << "rZK thérorique: " << rtheoric << std::endl;
    
    // Initialize first column
    for (int j = 0; j < M; j++) {
        u[j][0] = U0 / pow(cosh((x[j] - x0) * sqrt(U0/12) / mu), 2);
    }
    
    // First time step
    for (int j = 0; j < M; j++) {
        int jm, jp, jv, jw;
        periodicite(j, M, jm, jp, jv, jw);
        
        double A = u[j][0];
        double B = -((alpha/2) * (u[jp][0] + u[j][0] + u[jm][0]) * 
                    (u[jp][0] - u[jm][0]));
        double C = -((S/2) * (u[jw][0] - 2*u[jp][0] + 2*u[jm][0] - u[jv][0]));
        u[j][1] = A + B + C;
    }
    
    // Main time loop
    bool unstable = false;
    for (int n = 1; n < N && !unstable; n++) {
        for (int j = 0; j < M; j++) {
            int jm, jp, jv, jw;
            periodicite(j, M, jm, jp, jv, jw);
            
            double D = -3 * alpha * (u[jp][n] + u[j][n]) * (u[jp][n] - u[j][n]);
            double E = -2 * S * (u[jw][n] - 3*u[jp][n] + 3*u[j][n] - u[jm][n]);
            double F = -u[jp][n] + u[j][n] + u[jp][n-1];
            u[j][n+1] = D + E + F;
        }
        
        // Check for numerical instability
        double umax = 0;
        for (int j = 0; j < M; j++) {
            umax = std::max(umax, std::abs(u[j][n+1]));
        }
        if (umax > 1e6) {
            std::cout << "Instabilité numérique 6P: boum au pas " << n << std::endl;
            unstable = true;
        }
        if ((n+1) % 500 == 0) {
            std::cout << "Pas " << n+1 << " sur " << N << std::endl;
        }
        if ((n+1)==N){
        printf("Dernier pas atteint 6P\n");
        }
    }
    
    return u;
}