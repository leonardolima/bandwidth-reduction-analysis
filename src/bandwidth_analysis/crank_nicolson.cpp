#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include <map>
#include <Eigen/Dense>
#include "crank_nicolson.h"
#include "cuthill_mckee.h"
#include "gauss_jordan.h"
#include "lu_decomposition.h"
#include "../io.h"

/*******************************************************************************
 * Finite Difference Method for numerically solving the Heat Equation in 1
 * dimension.
 *
 *
 * @param N Number of grid points.
 * @param L Size of the grid.
 * @param dt Time step unit.
 * @param nsteps Number of steps.
 * @param R - Matrix used to store solution at each one of the steps.
 ******************************************************************************/
void diffusion_1d (int N, float L, float dt, int nsteps, Eigen::MatrixXf& R)
{
    // Assuming a square N x N grid
    // Numerical parameters, assuming heat coefficient = 1
    float dx = L/(N-1);      // Grid spacing
    float z = dt/pow(dx, 2); // Auxiliary variable

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
    // [1, N-4] rows
    for (int i = 1; i < N-3; ++i)
    {
        T(i, i-1) = -1;
        T(i, i) = 2;
        T(i, i+1) = -1;
    }

    // Initial condition, assuming u(0, x) = 2x for 0 <= x <= 1,
    // Boundary conditions, assuming u(0, 0) = 0 and u(0, N-3) = 0
    U[0] = 0;
    U[N-3] = 0;
    for (int i = 1; i < N-3; ++i) U[i] = 2*dx;

    // Defining auxiliary matrices A and B
    A = (I + (z/2)*T);
    B = (I - (z/2)*T);

    for (int m = 0; m < nsteps; ++m)
    {
        // b(i) = ( I - (z/2)*T) * U(:,i)
        for (int i = 0; i < N-2; ++i)
        {
            // Solving b for m
            b[i] = B.row(i)*U;
        }

        // Solving the system of eq. for m+1
        // ( I + (z/2)*T) * U = b
        U = A.colPivHouseholderQr().solve(b);
        // Boundary conditions
        U[0] = 0;
        U[N-3] = 0;

        R.col(m) = U;
    }
}

/*******************************************************************************
 * Save solutions (1 dimension case) for each time step as a CSV file.
 *
 *
 * @param N Number of grid points.
 * @param L Size of the grid.
 * @param dt Time step unit.
 * @param nsteps Number of time steps.
 ******************************************************************************/
void diffusion_1d_results_to_csv (int N, float L, float dt, int nsteps)
{
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N-2, nsteps);
    diffusion_1d(N, L, dt, nsteps, R);
    matrix_to_csv(R);
}

/*******************************************************************************
 * Finite Difference Method for numerically solving the Heat Equation in 2
 * dimensions.
 *
 *
 * @param N Number of grid points.
 * @param L Size of the grid.
 * @param dt Time step unit.
 * @param nsteps Number of steps.
 * @param R Matrix used to store solution at each one of the steps.
 ******************************************************************************/
void diffusion_2d (int N, float L, float dt, int nsteps, Eigen::MatrixXf& R)
{
    // Numerical parameters, assuming heat coefficient = 1
    int n = (int)sqrt(N-1);
    float dx = L/(N-1);                // Grid spacing, assuming a square grid
    float a = 1 + ((2*dt)/pow(dx, 2)); // Diagonal coefficient
    float c = -dt/(2*pow(dx, 2));      // Offdiagonal coefficient
    float d = 1 - (2*(dt/pow(dx, 2))); // Auxiliary coefficient used to solve b

    // Initializing matrice and vectors
    Eigen::MatrixXf H = Eigen::MatrixXf::Zero(N-2, N-2);
    Eigen::VectorXf U = Eigen::VectorXf::Zero(N-2);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(N-2);
    // Auxiliary matrix
    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(N-2, N-2);

    for (int i = 0; i < N-2; ++i)
    {
        // Boundary conditions, assuming homogeneous dirichlet conditions
        // u(0, y, t) = u(n, y, t) = 0 and
        // u(x, 0, t) = u(x, n, t) = 0
        if ((i < n) || (i > N-2-n)) U[i] = 0;
        else {
            if ((i % n) == 0)
            {
                U[i-1] = 0;
                U[i] = 0;
            } else {
                // Initial condition, assuming u(x, y, 0) = 2x + 2y
                U[i] = dx + dt;
            }
        }
    }

    // Defining H
    H = a*Eigen::MatrixXf::Identity(N-2, N-2);
    H.topRightCorner(N-3, N-3) += c*Eigen::MatrixXf::Identity(N-3, N-3);
    H.bottomLeftCorner(N-3, N-3) += c*Eigen::MatrixXf::Identity(N-3, N-3);
    H.topRightCorner(N-2-n, N-2-n) += c*Eigen::MatrixXf::Identity(N-2-n, N-2-n);
    H.bottomLeftCorner(N-2-n, N-2-n) += c*Eigen::MatrixXf::Identity(N-2-n, N-2-n);

    // Defining M
    M = d*Eigen::MatrixXf::Identity(N-2, N-2);
    M.topRightCorner(N-3, N-3) += -c*Eigen::MatrixXf::Identity(N-3, N-3);
    M.bottomLeftCorner(N-3, N-3) += -c*Eigen::MatrixXf::Identity(N-3, N-3);
    M.topRightCorner(N-2-n, N-2-n) += -c*Eigen::MatrixXf::Identity(N-2-n, N-2-n);
    M.bottomLeftCorner(N-2-n, N-2-n) += -c*Eigen::MatrixXf::Identity(N-2-n, N-2-n);

    // Removing some elements from the offdiagonals at positions
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

    // std::cout << "H = " << std::endl;
    // std::cout << H << std::endl;
    // std::cout << "M = " << std::endl;
    // std::cout << M << std::endl;

    for (int m = 0; m < nsteps; ++m)
    {
        // std::cout << "U[m] = " << std::endl;
        // std::cout << U << std::endl;

        // Mb = U (solving b for U at step m)
        b = M.partialPivLu().solve(U);

        // std::cout << "b[m] = " << std::endl;
        // std::cout << b << std::endl;

        // Solving the system of eq. for m+1
        // HU = b
        U = H.partialPivLu().solve(b);

        // for (int i = 0; i < N-2; ++i)
        // {
        //     // Solving b for m
        //     U[i] = H.row(i)*b;
        // }

        // Boundary conditions
        for (int i = 0; i < N-2; ++i)
        {
            if ((i < n) || (i > N-2-n)) U[i] = 0;
            else {
                if ((i % n) == 0)
                {
                    U[i-1] = 0;
                    U[i] = 0;
                }
            }
        }

        R.col(m) = U;
    }
}

/*******************************************************************************
 * Save solutions (2 dimensions case) for each time step as a CSV file.
 *
 *
 * @param N Number of grid points.
 * @param L Size of the grid.
 * @param dt Time step unit.
 * @param nsteps Number of time steps.
 ******************************************************************************/
void diffusion_2d_results_to_csv (int N, float L, float dt, int nsteps)
{
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N-2, nsteps);
    diffusion_2d(N, L, dt, nsteps, R);
    matrix_to_csv(R);
}

/*******************************************************************************
 * Finite Difference Method for numerically solving the Heat Equation in 2
 * dimensions.
 *
 *
 * @param X Number of grid points in the X direction.
 * @param Y Number of grid points in the Y direction.
 * @param m Mapping of grid points on the solutions vector.
 ******************************************************************************/
void map_to_U (int X, int Y, std::map<std::pair<int, int>, int>& m)
{
    int index_U = 0;

    // First submatrix (bottom part)
    for (int k = 0; k < X; ++k)
    {
        if (k % 2 == 0)
        {
            for (int l = 0; l < (int)Y/3; ++l)
            {
                m[{k, l}] = index_U++;
            }

        } else {
            for (int l = (int)((Y/3)-1); l >= 0; --l)
            {
                m[{k, l}] = index_U++;
            }
        }
    }

    // Second submatrix (middle part)
    for (int l = (int)(Y/3); l < (int)((2*(Y/3))+1); ++l)
    {
        if (l % 2 == 0)
        {
            for (int k = X-1; k >= (int)(X/2); --k)
            {
                m[{k, l}] = index_U++;
            }

        } else {
            for (int k = (int)(X/2); k < X; ++k)
            {
                m[{k, l}] = index_U++;
            }
        }
    }

    // Third submatrix (top part)
    for (int k = X-1; k >= 0; --k)
    {
        if (k % 2 == 0)
        {
            for (int l = (int)((2*(Y/3))+1); l < Y; ++l)
            {
                m[{k, l}] = index_U++;
            }

        } else {
            for (int l = Y-1; l >= (int)((2*(Y/3))+1); --l)
            {
                m[{k, l}] = index_U++;
            }
        }
    }
}

/*******************************************************************************
 * Finite Difference Method for numerically solving the Heat Equation in 2
 * dimensions on an irregular shape object. Satisfying the constraint Y = X + 3,
 * where X = {9, 11, 13, ...}.
 *
 *  C inversed shape:
 *   _____________
 *   |           |
 *   |           |
 *   |_____      |
 *         |     |
 *         |     |
 *   ______|     |
 *   |           |
 *   |           |
 *   |___________|
 *
 *
 * @param X Number of grid points in the X direction.
 * @param Y Number of grid points in the Y direction.
 * @param L Size of the grid.
 * @param dt Time step unit.
 * @param nsteps Number of steps.
 * @param R Matrix used to store solution at each one of the steps.
 * @param apply_cuthill_mckee If true applies the Cuthill-McKee algorithm to
 *                            matrices in order to improve performance.
 ******************************************************************************/
void diffusion_2d_irregular (int X, int Y, float L, float dt, int nsteps,
                             Eigen::MatrixXf& R, bool apply_cuthill_mckee)
{
    // Numerical parameters, assuming heat coefficient = 1
    int N = ((2*X*Y)/3)+((X/2)*(Y/3));           // Total number of grid points considered
    float dx = L/X;                              // Grid spacing, assuming dy = dx
    float a = 1 + ((2*dt)/pow(dx, 2));           // Diagonal coefficient
    float c = -dt/(2*pow(dx, 2));                // Offdiagonal coefficient
    float d = 1 - (2*(dt/pow(dx, 2)));           // Auxiliary coefficient used to solve b

    // Initializing matrice and vectors
    Eigen::MatrixXf H = Eigen::MatrixXf::Zero(N, N);
    Eigen::VectorXf U = Eigen::VectorXf::Zero(N);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(N);

    // Auxiliary matrix
    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(N, N);

    std::map<std::pair<int, int>, int> U_map;
    map_to_U(X, Y, U_map);

    // Defining H and M
    H = a*Eigen::MatrixXf::Identity(N, N);
    M = d*Eigen::MatrixXf::Identity(N, N);

    for (auto it = U_map.begin(); it != U_map.end(); ++it)
    {
        std::pair<int, int> p = it->first;
        int row = it->second;

        std::pair<int, int> p1;
        std::pair<int, int> p2;
        std::pair<int, int> p3;
        std::pair<int, int> p4;

        // TODO: Rewrite the if's below
        if (p.first+1 < X)
        {
            p1 = std::make_pair(p.first+1, p.second);
            auto p1_it = U_map.find(p1);
            if (p1_it != U_map.end())
            {
                int U1 = p1_it->second;
                H(row, U1) = c;
                M(row, U1) = -c;
            }
        }

        if (p.second+1 < Y)
        {
            p2 = std::make_pair(p.first, p.second+1);
            auto p2_it = U_map.find(p2);
            if (p2_it != U_map.end())
            {
                int U2 = p2_it->second;
                H(row, U2) = c;
                M(row, U2) = -c;
            }
        }

        if (p.first > 0)
        {
            p3 = std::make_pair(p.first-1, p.second);
            auto p3_it = U_map.find(p3);
            if (p3_it != U_map.end())
            {
                int U3 = p3_it->second;
                H(row, U3) = c;
                M(row, U3) = -c;
            }
        }
        if (p.second > 0)
        {
            p4 = std::make_pair(p.first, p.second-1);
            auto p4_it = U_map.find(p4);
            if (p4_it != U_map.end())
            {
                int U4 = p4_it->second;
                H(row, U4) = c;
                M(row, U4) = -c;
            }
        }
    }

    // Applying the Cuthill-McKee algorithm to H
    if (apply_cuthill_mckee)
    {
        Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N, N);
        Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N, N);
        run_algorithm(H, P, R);

        H = R;
        M = (P*M*P.transpose());
    }

    for (int m = 0; m < nsteps; ++m)
    {
        // std::cout << "U[m] = " << std::endl;
        // std::cout << U << std::endl;

        // Mb = U (solving b for U at step m)
        if (apply_cuthill_mckee) band_solver(M, b, U);
        else b = M.partialPivLu().solve(U);

        // std::cout << "b[m] = " << std::endl;
        // std::cout << b << std::endl;

        // Solving the system of eq. for m+1
        // HU = b
        if (apply_cuthill_mckee) band_solver(H, U, b);
        else U = H.partialPivLu().solve(b);

        R.col(m) = U;
    }
}

/*******************************************************************************
 * Compares time performance of the Crank-Nicolson approach for solving the 2D
 * diffusion problem before and after applying the Cuthill-McKee algorithm to
 * the matrices. Reducing their bandwidth improves the performance of each time
 * step when solving the linear system. An irregular shape object is assumed, in
 * particular a C inverted shape.
 *
 *
 * @param X Number of grid points in the X direction.
 * @param Y Number of grid points in the Y direction.
 * @param L Size of the grid.
 * @param dt Time step unit.
 * @param nsteps Number of steps.
 ******************************************************************************/
void diffusion_2d_irregular_compare (int X, int Y, float L, float dt, int nsteps)
{
    int N = ((2*X*Y)/3)+((X/2)*(Y/3));
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N, nsteps);

    auto start = std::chrono::high_resolution_clock::now();

    diffusion_2d_irregular(X, Y, L, dt, nsteps, R, false);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Execution time (without applying CM) = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;

    start = std::chrono::high_resolution_clock::now();

    diffusion_2d_irregular(X, Y, L, dt, nsteps, R, true);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(stop - start);
    duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Execution time (applying CM) = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;
}
