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
void diffusion_1d (int N, float L, float dt, int nsteps)
{
    // Initialize output file
    std::ofstream f("out.csv");

    if (!f.is_open()) return;

    // Numerical parameters, assuming heat coefficient = 1
    float dx = L/(N-1); // Grid spacing
    float z = dt/pow(dx, 2);

    // Initializing matrice and vectors
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
        std::string line;
        for (int i = 0; i < N-2; ++i) {

            if (i < N-3) line += std::to_string(U[i]) + ", ";
            if (i == N-3) line += std::to_string(U[i]);
        }

        f << line << std::endl;
    }

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
void diffusion_2d (int N, float L, float dt, int nsteps)
{
    // // Numerical parameters, assuming heat coefficient = 1
    // int n = (int)sqrt(N);
    // float dx = L/(N-1); // Grid spacing, assuming a square grid (i.e. dx = dy)
    // float a = 1 + dt/pow(dx, 2);
    // float c = -k/2*pow(dx, 2);

    // // Initializing matrice and vectors
    // Eigen::MatrixXf H = Eigen::MatrixXf::Zero(N-2, N-2);
    // Eigen::VectorXf U = Eigen::VectorXf::Zero(N-2);
    // Eigen::VectorXf b = Eigen::VectorXf::Zero(N-2);

    // // Initial condition, assuming u(0, x, y) = x + y
    // for (int i = 0; i < N-2; ++i) U[i] = dx + dt;

    // // Defining matrice and vectors
}
