#include <iostream>
#include <vector>
#include <cmath> //pow et sqrt
#include <chrono>

#include "../include/Zk.h"
#include "../include/functions.h"
#include "../include/fex.h"


//KdV_ZK scheme solver
std::vector<std::vector<double>> ZK_schema(double L, double T, int M, int N, double eta, double mu, double x0, double Uinf, double U0) {
    #ifdef DEBUG
    std::cout << "Calculating ZK scheme..." << std::endl;
    #endif
    auto start_timeZK = std::chrono::high_resolution_clock::now();
    
    double dx = L / M;  //spatial step
    double dt = T / N;  //temporal step
    
    auto x = linspace(0, L - dx, M);
    auto t = linspace(0, T, N + 1);
    
    //initialisation of u, vectorisation
    std::vector<std::vector<double>> u(M, std::vector<double>(N + 1, 0.0));
    
    //parameters
    double alpha = (eta*dt) / (3*dx);
    double S = ((pow(mu, 2))*dt) / (pow(dx, 3));
    
    #ifdef DEBUG
    double delta = 0.022;
    //print stability parameters
    double rnum = dt/dx;
    double rtheoric = 1 / (max_abs(fex(M, N, T, L, delta, U0, Uinf, x0)) + (4 * (pow(mu, 2))) / (pow(dx, 2)));
    std::cout << "rZK numérique: " << rnum << std::endl;
    std::cout << "rZK thérorique: " << rtheoric << std::endl;
    #endif 

    //initialise first column u[j,0]
    for (int j = 0; j < M; j++) {
        u[j][0] = U0 / pow(cosh((x[j] - x0) * sqrt(U0/12) / mu), 2);
    }
    
    //calculating u[j,1] using u[j,0]
    for (int j = 0; j < M; j++) {
        int jm, jp, jv, jw;
        periodicite(j, M, jm, jp, jv, jw);
        
        double A = u[j][0];
        double B = -((alpha/2) * (u[jp][0] + u[j][0] + u[jm][0]) * (u[jp][0] - u[jm][0]));
        double C = -((S/2) * (u[jw][0] - 2*u[jp][0] + 2*u[jm][0] - u[jv][0]));
        u[j][1] = A + B + C;
    }
    
    //main time loop
    bool unstable = false;
    for (int n = 1; n < N && !unstable; n++) {
        for (int j = 0; j < M; j++) {
            int jm, jp, jv, jw;
            periodicite(j, M, jm, jp, jv, jw);
            
            // Calcul de u[jm,n+1]
            double D = -alpha * (u[jp][n] + u[j][n] + u[jm][n]) * (u[jp][n] - u[jm][n]);
            double E = -S * (u[jw][n] - 2*u[jp][n] + 2*u[jm][n] - u[jv][n]);
            double F = u[j][n-1];
            u[j][n+1] = D + E + F;
        }
        
        // Check for numerical instability
        double umax = 0;
        for (int j = 0; j < M; j++) {
            umax = std::max(umax, std::abs(u[j][n+1]));
        }
        if (umax > 1e6) {
            std::cout << "Instabilité numérique Wang: boum au pas " << n << std::endl;
            unstable = true;
        }
        #ifdef DEBUG
        if ((n+1) % 500 == 0) {
            std::cout << "Pas " << n+1 << " sur " << N << std::endl;
        }
        #endif
        #ifdef DEBUG
        if ((n+1)==N){
        printf("Dernier pas atteint ZK\n");
        }
        #endif
    }

    auto end_timeZK = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_timeZK - start_timeZK;
    std::cout << "Temps de calcul modèle ZK = " << duration.count() << " millisecondes" << std::endl;
    
    return u;
}