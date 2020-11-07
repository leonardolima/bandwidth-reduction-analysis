#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>
#include "test.h"
#include "adapted_cuthill_mckee.h"
#include "level.h"
#include "topological.h"
#include "io.h"
#include "basic.h"

/*******************************************************************************
 * Changes constraints format.
 *
 * Constraints are given in the form L_x > [L_yi], but in order to keep search
 * O(1) the format is changed to L_x < [L_yi].
 *
 *
 * @param O     Ordering constraints (in the form L_x > [L_yi]).
 * @param new_O Ordering constraints (in the form L_x < [L_yi]).
 ******************************************************************************/
void convert_vector(const std::vector<std::vector<int>>& O,
                    std::vector<std::vector<int>>& new_O)
{
    for (std::vector<std::vector<int>>::size_type i = 0; i < O.size(); ++i)
    {
        for (std::vector<int>::size_type j = 0; j < O[i].size(); ++j)
        {
            new_O[O[i][j]].push_back(i);
        }
    }
}

void adapted_cuthill_mckee(const Eigen::MatrixXf& A, Eigen::MatrixXf& P, Eigen::MatrixXf& R,
                           std::vector<int>& path1, std::vector<int>& path2)
{
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Applying adapted Cuthill-McKee approach: " << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    apply_adapted_cuthill_mckee(A, P, path1);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Execution time = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;

    R = (P*A*P.transpose());

    apply_topological(A, path2);

    print_bandwidth_comparison(A, R);
    print_peak_comparison(A, path1, path2);
    std::cout << "-------------------------------------------" << std::endl;
}

void levels(const Eigen::MatrixXf& A, Eigen::MatrixXf& P, Eigen::MatrixXf& R,
            std::vector<int>& path1, std::vector<int>& path2)
{
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Applying level approach: " << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    apply_levels(A, P, path1);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Execution time = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;

    R = (P*A*P.transpose());

    apply_topological(A, path2);

    print_bandwidth_comparison(A, R);
    print_peak_comparison(A, path1, path2);
    std::cout << "-------------------------------------------" << std::endl;
}

void topological(const Eigen::MatrixXf& A, Eigen::MatrixXf& P, Eigen::MatrixXf& R,
                 std::vector<int>& path)
{
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Applying topological approach: " << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    apply_topological_all(A, P, path);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::cout << "Execution time = " << duration_s.count() << "s" << std::endl;
    std::cout << std::endl;

    R = (P*A*P.transpose());

    print_bandwidth_comparison(A, R);
    print_peak(A, path);
    std::cout << "-------------------------------------------" << std::endl;
}

void test_adapted_cuthill_mckee(const std::string& file_name)
{
    int N = read_N_from_file(file_name);

    // Adjacency matrix representing the graph
    Eigen::MatrixXf A = Eigen::MatrixXf::Identity(N, N);

    // Permutation matrix
    Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N, N);

    // Resulting matrix
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N, N);

    // Paths
    std::vector<int> path1;
    std::vector<int> path2;

    read_file(file_name, A, N);

    // 1. Topological sorting
    adapted_cuthill_mckee(A, P, R, path1, path2);
}

void test_levels(const std::string& file_name)
{
    int N = read_N_from_file(file_name);

    // Adjacency matrix representing the graph
    Eigen::MatrixXf A = Eigen::MatrixXf::Identity(N, N);

    // Permutation matrix (in case of Cuthill-McKee algorithm)
    Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N, N);

    // Resulting matrix
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N, N);

    // Paths
    std::vector<int> path1;
    std::vector<int> path2;

    read_file(file_name, A, N);

    // 1. Topological sorting
    levels(A, P, R, path1, path2);
}

void test_topological(const std::string& file_name)
{
    int N = read_N_from_file(file_name);

    // Adjacency matrix representing the graph
    Eigen::MatrixXf A = Eigen::MatrixXf::Identity(N, N);

    // Permutation matrix
    Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N, N);

    // Resulting matrix
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N, N);

    // Paths
    std::vector<int> path;

    read_file(file_name, A, N);

    // 1. Topological sorting
    topological(A, P, R, path);
}

// void test_all(const std::string& file_name)
// {
//     int N = read_N_from_file(file_name);

//     // Adjacency matrix representing the graph
//     Eigen::MatrixXf A = Eigen::MatrixXf::Identity(N, N);

//     // Permutation matrix
//     Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N, N);

//     // Resulting matrix
//     Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N, N);

//     // Paths
//     std::vector<int> path;

//     read_file(file_name, A, N);

//     // 1. Topological sorting
//     topological(A, P, R, path);
//     path.clear();

//     // 2. Level approach
//     levels(A, P, R, path);
//     path.clear();

//     // 3. Adapted Cuthill-McKee approach
//     adapted_cuthill_mckee(A, P, R, path);
//     path.clear();
// }
