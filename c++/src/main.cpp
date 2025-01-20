#include <iostream>
#include <vector>
#include <cmath> //pow et sqrt

#include "../include/Wang.h"
#include "../include/Zk.h"
#include "../include/functions.h"
#include "../include/fex.h"

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
    
    //compute exact solution
    auto u_exact = fex(M, N, T, L, 0.022, U0, Uinf, x0);

    //run Wang solver
    auto u_w = W_schema(L, T, M, N, eta, mu, x0, Uinf, U0);
    //residu
    double residu = compute_residu(M, N, u_w, u_exact);
    std::cout << "Résidu avec Wang : " << residu << "\n" << std::endl;

    //run ZK solver
    auto u_zk = ZK_schema(L, T, M, N, eta, mu, x0, Uinf, U0);
    //residu
    residu = compute_residu(M, N, u_zk, u_exact);
    std::cout << "Résidu avec ZK : " << residu << "\n" << std::endl;


    return 0;
}