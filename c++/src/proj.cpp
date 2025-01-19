#include <iostream>
#include <vector>
#include <cmath> //pow et sqrt


//creating the linspace function
std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + step * i;
    }
    return result;
}

//periodic boundary conditions
void periodicite(int j, int M, int& jm, int& jp, int& jv, int& jw) {
    jm = (j - 1 + M) % M;   //using modulo 
    jp = (j + 1) % M;
    jv = (j - 2 + M) % M;
    jw = (j + 2) % M;
}

//KdV_Wang scheme solver
std::vector<std::vector<double>> W_schema(double L, double T, int M, int N, double eta, double mu, double x0, double Uinf, double U0) {
    std::cout << "Calculating Wang scheme..." << std::endl;
    
    double dx = L / M;  //spatial step
    double dt = T / N;  //temporal step
    
    auto x = linspace(0, L - dx, M);
    auto t = linspace(0, T, N + 1);
    
    //initialisation of u, vectorisation
    std::vector<std::vector<double>> u(M, std::vector<double>(N + 1, 0.0));
    
    //parameters
    double alpha = (eta*dt) / (3*dx);
    double S = ((pow(mu, 2))*dt) / (pow(dx, 3));
    double delta = 0.022;
    double DELTA = delta / (sqrt((U0 - Uinf)/12));
    
    //print stability parameters
    double rnum = dt/dx;
    std::cout << "rW numérique: " << rnum << std::endl;
    //std::cout << "rW thérorique: " << rtheoric << std::endl;

    //initialise first column u[j,0]
    for (int j = 0; j < M; j++) {
        u[j][0] = U0 / pow(cosh((x[j] - x0)/DELTA), 2);
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
            double F = u[jm][n] - u[jp][n] + u[jp][n-1];
            u[jm][n+1] = D + E + F;
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
        if ((n+1) % 500 == 0) {
            std::cout << "Pas " << n+1 << " sur " << N << std::endl;
        }
    }
    
    return u;
}


int main() {
    //numerical parameters
    double eta = 1;
    double mu = 0.022;
    double L = 2;
    double T = 10;
    double Uinf = 0;
    double U0 = 1;
    int M = 100;
    int N = 3000;
    double x0 = L/2;
    
    //run solver
    auto u_w = W_schema(L, T, M, N, eta, mu, x0, Uinf, U0);

    return 0;
}