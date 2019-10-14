/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>
#include "gauss_jordan.h"

bool node_has_label (const Eigen::MatrixXf& P, unsigned int col)
{
    for (int i = 0; i < P.rows(); ++i)
    {
        if (P(i,col) != 0) return true;
    }

    return false;
}

std::vector<std::pair<int, int>> sort_row_deg (const Eigen::VectorXf& node_row, 
                                               const std::vector<int>& row_deg, 
                                               const Eigen::MatrixXf& P)
{
    std::vector<std::pair<int, int>> sorted_row_degs;

    for(int i = 0; i < node_row.size(); ++i)
    {
        // We are not interested if the element is on the diagonal
        // or if the node already has a label
        if (node_row[i] != 0 && !node_has_label(P, i))
        {
            // std::cout << "(" << row_deg[i] << ", " << i << ") ";
            std::pair<int, int> pair = std::make_pair(row_deg[i], i);
            sorted_row_degs.push_back(pair);
        }
    }

    std::sort(sorted_row_degs.begin(), sorted_row_degs.end());

    return sorted_row_degs;
}

// Returns index (column or row, it doesn't matter) in A
int node_index (const Eigen::MatrixXf& P, int label)
{
    // index in A corresponds to column (j)
    for(int j = 0; j < P.cols(); ++j)
    {
        if (P(label,j) != 0) return j;
    }

    return -1;
}

void nodal_numbering (const Eigen::MatrixXf& A, Eigen::MatrixXf& P, const std::vector<int>& row_deg)
{
    unsigned int new_label = 0;

    // 3. For every node we (i) list the adjacent nodes (ii) sort them in ascending order
    // of degree (iii) label them uniquely and accordingly (2, 3, ...)
    for(int label = 0; label < A.cols(); ++label)
    {
        int node_row_index = node_index(P, label);

        Eigen::VectorXf node_row = A.row(node_row_index);
        node_row[node_row_index] = 0; // Diagonal element should not be considered

        std::vector<std::pair<int, int>> sorted_row_degs = sort_row_deg(node_row, row_deg, P);

        // std::cout << "------==-------" << std::endl;
        // std::cout << "label = " << label << std::endl;
        // std::cout << "node_row_index = " << node_row_index << std::endl;
        // std::cout << "node_row = " << node_row << std::endl;
        // std::cout << "sorted_row_degs = " << sorted_row_degs.size() << std::endl;
        // for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_row_degs.size(); ++i)
        // {
        //     std::cout << "(" << sorted_row_degs[i].first << ", " << sorted_row_degs[i].second << ") ";
        // }
        // std::cout << std::endl;

        for (std::vector<std::pair<int, int>>::size_type j = 0; j < sorted_row_degs.size(); ++j)
        {
            // Update matrix P accordingly
            ++new_label;
            P(new_label, sorted_row_degs[j].second) = 1;
            // std::cout << "P = " << std::endl;
            // std::cout << P << std::endl;
        }
    }
}

int matrix_deg (Eigen::MatrixXf& R)
{
    std::vector<int> row_deg(R.rows(), 0); // Assuming R is square

    for(int j = 0; j < R.cols(); ++j)
    {
        for (int i = 0; i < R.rows(); ++i)
        {
            if (R(i,j) != 0 && i != j) row_deg[i] += 1;
        }
    }

    auto max_deg = *std::max_element(row_deg.begin(), row_deg.end());

    return max_deg;
}

void compareMatrices(const Eigen::MatrixXf& A, const Eigen::MatrixXf& R)
{
    Eigen::MatrixXf A_inv(A.rows(), A.cols());
    Eigen::MatrixXf R_inv(R.rows(), R.cols());
    A_inv = A;
    R_inv = R;
    
    auto start_A = std::chrono::high_resolution_clock::now();
    //A.inverse();
    compute_inverse(A_inv);
    auto stop_A = std::chrono::high_resolution_clock::now();
    auto duration_A = std::chrono::duration<float>(stop_A - start_A);
    auto duration_ns_A = std::chrono::duration_cast<std::chrono::nanoseconds>(duration_A);
    std::cout << "A.inverse() duration = " << duration_ns_A.count() << "ns" << std::endl;

    auto start_R = std::chrono::high_resolution_clock::now();
    compute_inverse(R_inv);
    auto stop_R = std::chrono::high_resolution_clock::now();
    auto duration_R = std::chrono::duration<float>(stop_R - start_R);
    auto duration_ns_R = std::chrono::duration_cast<std::chrono::nanoseconds>(duration_R);
    std::cout << "R.inverse() duration = " << duration_ns_R.count() << "ns" << std::endl;
    std::cout << std::endl;
}

std::vector<int> compute_row_deg (const Eigen::MatrixXf& A)
{
    std::vector<int> row_deg(A.rows(), 0); // Assuming A is square

    for(int j = 0; j < A.cols(); ++j)
    {
        for (int i = 0; i < A.rows(); ++i)
        {
            if (A(i,j) != 0 && i != j) row_deg[i] += 1;
        }
    }

    return row_deg;
}

std::vector<int> select_starting_nodes (Eigen::MatrixXf& A)
{
    std::vector<int> starting_nodes;

 
    // 1. Compute degree of every row, where deg(i) = \Sum_{j!=i} a_{ij} != 0   
    std::vector<int> row_deg = compute_row_deg(A);

    // for(std::vector<int>::size_type i = 0; i < row_deg.size(); ++i)
    // {
    //     std::cout << row_deg[i] << ' ';
    // }

    // 2. Select starting node (usually the one with minimum degree)
    // Permuting rows and columns of matrix A is done using the
    // permutation matrix P. P has one non-zero unit element in each
    // row and column. For row i: p_{ij} i is the new node label and
    // j is the original one. We should always keep in mind that
    // optimal starting nodes doesn't necessarely are the lowest
    // degree ones
    auto min_deg = *std::min_element(row_deg.begin(), row_deg.end());
    for(std::vector<int>::size_type i = 0; i < row_deg.size(); ++i)
    {
        if (row_deg[i] == min_deg) starting_nodes.push_back(i);
    }

    return starting_nodes;
}

void compute_results (std::vector<int>& starting_nodes, const Eigen::MatrixXf& A,
                      Eigen::MatrixXf& P, Eigen::MatrixXf& R, Eigen::MatrixXf& minR)
{   
    int max_mat_deg = std::numeric_limits<int>::max(); // max_mat_deg = \inf

    std::vector<int> row_deg = compute_row_deg(A);
    
    for(std::vector<int>::size_type i = 0; i < starting_nodes.size(); ++i)
    {
        P(0, starting_nodes[i]) = 1; // Starting node is labeled as 1
        // std::cout << std::endl;
        // std::cout << "Starting node: " << starting_nodes[i] << std::endl;

        nodal_numbering(A, P, row_deg);
        // std::cout << "P = " << std::endl;
        // std::cout << P << std::endl;
        // std::cout << std::endl;

        R = (P*A*P.transpose());
        int mat_deg = matrix_deg(R);
        // std::cout << "Matrix degree = " << mat_deg << std::endl;
        // std::cout << "P*A*P^T = " << std::endl;
        // std::cout << "A B C D E F G H I J" << std::endl;
        // std::cout << R << std::endl;
        // std::cout << std::endl;


        if (mat_deg < max_mat_deg) minR = R;

        P.setZero(); // Clear P matrix
    }

    // std::cout << "Lower degree matrix: " << std::endl;
    // std::cout << minR << std::endl;

    compareMatrices(A, R);
}

void generate_binary_random_matrix (Eigen::MatrixXf& A, int dimension)
{
    A = (A + Eigen::MatrixXf::Constant(dimension, dimension, 1.))*(1./2.);
    A = (A + Eigen::MatrixXf::Constant(dimension, dimension, 0.));
    A = A.array().round();
}
