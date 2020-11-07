#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <iostream>
#include "basic.h"

/*******************************************************************************
 * Calculates bandwidth of a particular matrix.
 *
 *
 * @param R Matrix.
 ******************************************************************************/
int matrix_bandwidth (const Eigen::MatrixXf& R)
{
    std::vector<int> rows_bandwidth(R.rows(), 0);

    for (int i = 0; i < R.rows(); ++i)
    {
        for(int j = 0; j < R.cols(); ++j)
        {
            if (R(i,j) != 0 && i != j && (fabs(j-i) > rows_bandwidth[i]))
            {
                rows_bandwidth[i] = fabs(j-i);
            }
        }
    }

    auto max_bandwidth = *std::max_element(rows_bandwidth.begin(), rows_bandwidth.end());

    return max_bandwidth;
}

void parents_from_matrix(const Eigen::MatrixXf& A, std::vector<std::vector<int>>& parents)
{
    int N = A.rows();

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (A(i, j) != 0 && i != j) parents[i].push_back(j);
        }
    }

    // for (std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    // {
    //     std::cout << i << " -> ";
    //     for (std::vector<int>::size_type j = 0; j < parents[i].size(); ++j)
    //     {
    //         std::cout << parents[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

void children_from_matrix(const Eigen::MatrixXf& A, std::vector<std::vector<int>>& children)
{
    int N = A.rows();

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (A(i, j) != 0 && i != j) children[j].push_back(i);
        }
    }
}

void mem_from_matrix(const Eigen::MatrixXf& A, std::vector<int>& mem)
{
    int N = A.rows();

    for (int i = 0; i < N; ++i)
    {
        mem[i] = A(i,i);
    }
}

void apply_symmetry (Eigen::MatrixXf& A)
{
    for (int i = 0; i < A.rows(); ++i)
        {
            for (int j = 0; j < A.cols(); ++j)
                {
                    if (A(i,j) != 0) A(j, i) = A(i, j);
                }
        }
}

/*******************************************************************************
 * Calculates the degree of each node (by definition it is the number of
 * non-zero off-diagonal elements).
 *
 *
 * @param A Adjacency matrix.
 ******************************************************************************/
std::vector<int> compute_nodes_deg(const Eigen::MatrixXf& A)
{
    std::vector<int> nodes_deg(A.rows(), 0); // A is square

    for(int j = 0; j < A.cols(); ++j)
    {
        for (int i = 0; i < A.rows(); ++i)
        {
            if (A(i,j) != 0 && i != j) nodes_deg[i] += 1;
        }
    }

    return nodes_deg;
}

void label_sorted_nodes(const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                        std::vector<std::pair<int, int>> sorted_pairs)
{

    // P_{ij}: where i is the new node label and j is the original one
    for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_pairs.size(); ++i)
        {
            P(sorted_pairs[i].first, sorted_pairs[i].second) = 1; // Update matrix P accordingly
        }
}

std::vector<std::pair<int, int>> make_sorted_pairs(const std::vector<int> sorted_nodes)
{
    std::vector<std::pair<int, int>> sorted_pairs;

    for (std::vector<int>::size_type i = 0; i < sorted_nodes.size(); ++i)
        {
            std::pair<int, int> pair = std::make_pair(i, sorted_nodes[i]);
            sorted_pairs.push_back(pair);
        }

    std::sort(sorted_pairs.begin(), sorted_pairs.end());

    return sorted_pairs;
}
