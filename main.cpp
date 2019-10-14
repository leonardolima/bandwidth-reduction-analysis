/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include <iostream>
#include <chrono>
#include <vector>
#include <Eigen/Dense>
#include "bandwidth_minimization.h"
#include "gauss_jordan.h"

int main (void)
{
    // Eigen::MatrixXf A(10, 10);
    // A << 1, 1, 0, 1, 0, 0, 0, 0, 1, 0,
    //      1, 1, 1, 0, 0, 0, 0, 0, 1, 0,
    //      0, 1, 1, 0, 1, 0, 0, 0, 1, 0,
    //      1, 0, 0, 1, 1, 1, 0, 0, 1, 1,
    //      0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    //      0, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    //      0, 0, 0, 0, 0, 1, 1, 1, 0, 1,
    //      0, 0, 0, 0, 1, 0, 1, 1, 0, 1,
    //      1, 1, 1, 1, 1, 0, 0, 0, 1, 1,
    //      0, 0, 0, 1, 1, 1, 1, 1, 1, 1;

    // Eigen::MatrixXf A(5, 5);
    // A << 1, 1, 0, 1, 0,
    //      1, 1, 1, 0, 1,
    //      0, 1, 1, 0, 0,
    //      1, 0, 0, 1, 1,
    //      0, 1, 0, 1, 1;

    int dimension = 100;
    Eigen::MatrixXf A = Eigen::MatrixXf::Random(dimension, dimension);
    generate_binary_random_matrix(A, dimension);

    std::vector<int> starting_nodes = select_starting_nodes(A);

    Eigen::MatrixXf P = Eigen::MatrixXf::Zero(A.rows(), A.cols());
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(A.rows(), A.cols());
    Eigen::MatrixXf minR = Eigen::MatrixXf::Zero(A.rows(), A.cols());

    auto start = std::chrono::high_resolution_clock::now();

    compute_results (starting_nodes, A, P, R, minR);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Execution duration = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;    
}
