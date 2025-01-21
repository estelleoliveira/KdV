#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

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


// class GnuplotPipe for graphics
class GnuplotPipe {
private:
    FILE* pipe;
    
public:
    GnuplotPipe() {
        #ifdef _WIN32
            pipe = _popen("gnuplot -persist", "w");
        #else
            pipe = popen("gnuplot -persist", "w");
        #endif
        if (!pipe) {
            throw std::runtime_error("Gnuplot not found!");
        }
    }
    
    ~GnuplotPipe() {
        if (pipe) {
            #ifdef _WIN32
                _pclose(pipe);
            #else
                pclose(pipe);
            #endif
        }
    }
    
    void operator()(const std::string& command) {
        fprintf(pipe, "%s\n", command.c_str());
        fflush(pipe);
    }
};

//function to create a gnuplot script
void plotResults(const std::string& wang_file, const std::string& zk_file, const std::string& sixp_file, const std::string& exact_file) {
    try {
        GnuplotPipe gp;

        gp("set terminal pngcairo size 800,600 enhanced font 'Verdana,10'");
        
        // configuration
        gp("set terminal wxt");
        gp("set grid lt 1 lw 0.5 lc rgb '#cccccc'"); //smooth grid
        gp("set title 'Solutions KdV'");
        gp("set xlabel 'x'");
        gp("set ylabel 'u'");
        gp("set samples 2000");

        // Customise the line styles
        gp("set style line 1 linecolor rgb '#FF0000' linetype 1 linewidth 1.2");  // Red line, solid, width 1
        gp("set style line 2 linecolor rgb '#FFD858' linetype 2 linewidth 1");  // Yellow dashed line, width 1
        gp("set style line 3 linecolor rgb '#0000FF' linetype 1 linewidth 1.2");  // Blue line, solid, width 1
        gp("set style line 4 linecolor rgb '#32CD32' linetype 1 linewidth 1.2");  // Lime line, solid, width 1

        // Animation loop to generate frames
        for (int i = 0; i <= 150; ++i) {
            std::string plot_cmd = "plot '" + wang_file + "' using 1:(column(" + std::to_string(i + 2) + ")) title 'Wang' with lines linestyle 1, "
                                   "'" + zk_file + "' using 1:(column(" + std::to_string(i + 2) + ")) title 'ZK' with lines linestyle 3, "
                                   "'" + sixp_file + "' using 1:(column(" + std::to_string(i + 2) + ")) title '6P' with lines linestyle 4, "
                                   "'" + exact_file + "' using 1:(column(" + std::to_string(i + 2) + ")) title 'exact' with lines linestyle 2";
            gp(plot_cmd);
        }

        #ifdef DEBUG
        gp("set print 'gnuplot.log'");
        gp("show plot");
        #endif
        
    } catch (const std::exception& e) {
        std::cerr << "Erreur de plotting: " << e.what() << std::endl;
        std::cerr << "Assurez-vous que gnuplot est installé sur votre système." << std::endl;
    }
}

// Function to save results to a file
void saveToFile(const std::vector<std::vector<double>>& u, const std::vector<double>& x, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    // for each point
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i];
        for (size_t n = 0; n < u[0].size(); ++n) {
            file << " " << u[i][n];
        }
        file << "\n";
    }
    file.close();
}