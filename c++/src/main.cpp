#include <iostream>
#include <vector>
#include <cmath> //pow et sqrt

#include "../include/Wang.h"
#include "../include/Zk.h"
#include "../include/sixP.h"
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

    //run 6P solver
    auto u_sixp = sixP_schema(L, T, M, N, eta, mu, x0, Uinf, U0);
    //residu
    residu = compute_residu(M, N, u_sixp, u_exact);
    std::cout << "Résidu avec 6P : " << residu << "\n" << std::endl;


    // Save results to CSV files
    auto x = linspace(0, L - L/M, M);

    std::cout << "\nSaving results..." << std::endl;
    saveToFile(u_w, x, "kdv_w_results.txt");
    saveToFile(u_zk, x, "kdv_zk_results.txt");
    saveToFile(u_sixp, x, "kdv_6p_results.txt");
    saveToFile(u_exact, x, "kdv_exact_results.txt");
    std::cout << "Results saved to kdv_w_results.txt, kdv_zk_results.txt, kdv_6p_results.txt and kdv_exact_results.txt" << std::endl;

    std::cout << "\nCréation des visualisations..." << std::endl;
    plotResults("kdv_w_results.txt", "kdv_zk_results.txt", "kdv_6p_results.txt", "kdv_exact_results.txt");

    return 0;
}