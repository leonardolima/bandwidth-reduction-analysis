/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <string>
#include <chrono>
#include "bandwidth_minimization.h"

std::string generate_csv_line (const Eigen::VectorXf& U)
{
    std::string line;

    for (int i = 0; i < U.size(); ++i) {    
        if (i < U.size()-1) line += std::to_string(U[i]) + ", ";
        if (i == U.size()-1) line += std::to_string(U[i]);
    }
    
    return line;
}

/* Function: diffusion_1d
 *
 * Inputs: N - number of grid points (assuming a square n x n grid)
 *         L - size of the grid
 *         dt - time step
 *         nsteps - number of time steps
 *
 * Finite Difference Method for numerically solving the
 * Heat Equation in 1 dimension
 */
void diffusion_1d (int N, float L, float dt, int nsteps, std::ofstream& f)
{
    // Numerical parameters, assuming heat coefficient = 1
    float dx = L/(N-1); // Grid spacing
    float z = dt/pow(dx, 2);

    // Initializing matrices and vectors
    Eigen::MatrixXf T = Eigen::MatrixXf::Zero(N-2, N-2);
    Eigen::MatrixXf I = Eigen::MatrixXf::Identity(N-2, N-2);
    Eigen::VectorXf U = Eigen::VectorXf::Zero(N-2);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(N-2);
    // Auxiliary matrices
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N-2, N-2);
    Eigen::MatrixXf B = Eigen::MatrixXf::Zero(N-2, N-2);

    // Defining T
    // First row
    T(0, 0) = 2;
    T(0, 1) = -1;
    // Last row
    T(N-3, N-3) = 2;
    T(N-3, N-4) = -1;
    for (int i = 1; i < N-3; ++i)
    {
        T(i, i-1) = -1;
        T(i, i) = 2;
        T(i, i+1) = -1;
    }

    // Initial condition, assuming u(0, x) = 2x for 0 <= x <= 1
    for (int i = 0; i < N-2; ++i) U[i] = 2*dx;

    // std::cout << "T = " << std::endl;
    // std::cout << T << std::endl << std::endl;

    // Defining auxiliary matrices A and B
    A = (I + (z/2)*T);
    B = (I - (z/2)*T);

    // std::cout << "A = " << std::endl;
    // std::cout << A << std::endl << std::endl;
    // std::cout << "B = " << std::endl;
    // std::cout << B << std::endl;
    
    for (int m = 0; m < nsteps; ++m)
    {
        // b(i) = ( I - (z/2)*T) * U(:,i) 
        for (int i = 0; i < N-2; ++i)
        {
            b[i] = B.row(i)*U;
        }

        // Solving the system of eq. for m+1
        // ( I + (z/2)*T) * U = b
        U = A.colPivHouseholderQr().solve(b);

        // CSV output
        std::string line = generate_csv_line(U);

        f << line << std::endl;
    }
}

void diffusion_1d_results_to_csv (int N, float L, float dt, int nsteps)
{
    // Initialize output file
    std::ofstream f("out.csv");

    if (!f.is_open()) return;

    diffusion_1d(N, L, dt, nsteps, f);

    f.close();
}

/* Function: diffusion_2d
 *
 * Inputs: N - number of grid points (assuming a square n x n grid)
 *         L - size of the grid
 *         dt - time step
 *         nsteps - number of time steps
 *
 * Finite Difference Method for numerically solving the
 * Heat Equation in 2 dimensions
 */
void diffusion_2d (int N, float L, float dt, int nsteps, std::ofstream& f, bool apply_cuthill_mckee)
{
    // Numerical parameters, assuming heat coefficient = 1
    int n = (int)sqrt(N-1);
    float dx = L/(N-1); // Grid spacing, assuming a square grid (i.e. dx = dy)
    float a = 1 + dt/pow(dx, 2);
    float c = -dt/(2*pow(dx, 2));
    float d = 1 - (2*dt/pow(dx, 2));
    float e = dt/(2*pow(dx, 2));

    // Initializing matrice and vectors
    Eigen::MatrixXf H = Eigen::MatrixXf::Zero(N-2, N-2);
    Eigen::VectorXf U = Eigen::VectorXf::Zero(N-2);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(N-2);
    // Auxiliary matrix
    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(N-2, N-2);

    // Initial condition, assuming u(0, x, y) = x + y
    for (int i = 0; i < N-2; ++i) U[i] = dx + dt;

    // Defining H
    H = a*Eigen::MatrixXf::Identity(N-2, N-2);
    H.topRightCorner(N-3, N-3) += c*Eigen::MatrixXf::Identity(N-3, N-3);
    H.bottomLeftCorner(N-3, N-3) += c*Eigen::MatrixXf::Identity(N-3, N-3);
    H.topRightCorner(N-2-n, N-2-n) += c*Eigen::MatrixXf::Identity(N-2-n, N-2-n);
    H.bottomLeftCorner(N-2-n, N-2-n) += c*Eigen::MatrixXf::Identity(N-2-n, N-2-n);

    // Defining M
    M = d*Eigen::MatrixXf::Identity(N-2, N-2);
    M.topRightCorner(N-3, N-3) += e*Eigen::MatrixXf::Identity(N-3, N-3);
    M.bottomLeftCorner(N-3, N-3) += e*Eigen::MatrixXf::Identity(N-3, N-3);
    M.topRightCorner(N-2-n, N-2-n) += e*Eigen::MatrixXf::Identity(N-2-n, N-2-n);
    M.bottomLeftCorner(N-2-n, N-2-n) += e*Eigen::MatrixXf::Identity(N-2-n, N-2-n);

    // Removing elements from the offdiagonals at positions 
    // (n, n-1), (n-1, n) and so on
    for (int i = 0; i < N-2; ++i)
    {
        if ((i % n) == 0 && (i > 0))
        {
            H(i, i-1) = 0;
            H(i-1, i) = 0;
            M(i, i-1) = 0;
            M(i-1, i) = 0;
        }
    }

    std::cout << "H = " << std::endl;
    std::cout << H << std::endl;

    // Applying the Cuthill-McKee algorithm to H
    if (apply_cuthill_mckee)
    {
        Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N-2, N-2);
        Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N-2, N-2);
        run_algorithm(H, P, R);

        std::cout << "R = " << std::endl;
        std::cout << R << std::endl;
        H = R;
    }
    
    for (int m = 0; m < nsteps; ++m)
    {
        for (int i = 0; i < N-2; ++i)
        {
            b[i] = M.row(i)*U;
        }
 
        // Solving the system of eq. for m+1
        U = H.fullPivLu().solve(b);
 
        // CSV output
        std::string line = generate_csv_line(U);

        f << line << std::endl;
    }
}

void diffusion_2d_results_to_csv (int N, float L, float dt, int nsteps)
{
    // Initialize output file
    std::ofstream f("out.csv");

    if (!f.is_open()) return;

    diffusion_2d(N, L, dt, nsteps, f, false);

    f.close();
}

void diffusion_2d_compare (int N, float L, float dt, int nsteps)
{
    std::ofstream f("out.csv");

    auto start = std::chrono::high_resolution_clock::now();

    diffusion_2d(N, L, dt, nsteps, f, false);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Execution duration = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;

    start = std::chrono::high_resolution_clock::now();

    diffusion_2d(N, L, dt, nsteps, f, true);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(stop - start);
    duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Applying the Cuthill-McKee algorithm: " << std::endl;
    std::cout << "Execution duration = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;
}
