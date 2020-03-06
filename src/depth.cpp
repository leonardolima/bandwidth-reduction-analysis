#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <Eigen/Dense>
#include "basic.h"

void label_nodes(const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                 std::vector<std::pair<int, int>> sorted_depth)
{
    int new_label = 0;

    // P_{ij}: where i is the new node label and j is the original one
    for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_depth.size(); ++i)
    {
        P(new_label++, sorted_depth[i].second) = 1; // Update matrix P accordingly
    }
}

std::vector<std::pair<int, int>> sort_depth(const std::vector<int> depth)
{
    std::vector<std::pair<int, int>> sorted_depth;

    for (std::vector<int>::size_type i = 0; i < depth.size(); ++i)
    {
        std::pair<int, int> pair = std::make_pair(depth[i], i);
        // std::cout << "(depth[i], i) = (" << depth[i] << ", " << i << ")" << std::endl;
        sorted_depth.push_back(pair);
    }

    std::sort(sorted_depth.begin(), sorted_depth.end());

    return sorted_depth;
}

void calculate_depth(const std::vector<std::vector<int>>& parents, std::vector<int>& depth)
{
    int N = parents.size();

    std::vector<int> visited_nodes(N, 0);

    // 1. Find all root nodes
    for (std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if (parents[i].size() == 0) depth[i] = 0;
    }

    // 2. Assuming graph is acyclic (even if not it might work)
    for (std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if (parents[i].size() > 0)
        {
            auto min = std::numeric_limits<int>::max();

            for (std::vector<int>::size_type j = 0; j < parents[i].size(); ++j)
            {
                if (depth[parents[i][j]] + 1 < min) min = depth[parents[i][j]] + 1;
            }

            depth[i] = min;
        }
    }
}

void apply_depth(const Eigen::MatrixXf& A, Eigen::MatrixXf& P)
{
    // TODO: Declaring a variable N at each function doesn't seem to be
    // the best approach. Fix it.
    int N = A.rows(); // Number of nodes

    std::vector<std::vector<int>> parents(N), children(N);

    parents_from_matrix(A, parents);
    children_from_matrix(A, children);

    std::vector<int> depth(N, -1);

    calculate_depth(parents, depth);

    std::vector<std::pair<int, int>> sorted_depth = sort_depth(depth);

    label_nodes(A, P, sorted_depth);
}
