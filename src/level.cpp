#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "level.h"
#include "basic.h"
#include "io.h"

void label_nodes(const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                 std::vector<std::pair<int, int>> sorted_levels)
{
    int new_label = 0;

    // P_{ij}: where i is the new node label and j is the original one
    for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_levels.size(); ++i)
    {
        P(new_label++, sorted_levels[i].second) = 1; // Update matrix P accordingly
    }
}

std::vector<std::pair<int, int>> sort_levels(const std::vector<int> levels)
{
    std::vector<std::pair<int, int>> sorted_levels;

    for (std::vector<int>::size_type i = 0; i < levels.size(); ++i)
    {
        std::pair<int, int> pair = std::make_pair(levels[i], i);
        // std::cout << "(levels[i], i) = (" << levels[i] << ", " << i << ")" << std::endl;
        sorted_levels.push_back(pair);
    }

    std::sort(sorted_levels.begin(), sorted_levels.end());

    return sorted_levels;
}

void assign_level(const std::vector<std::vector<int>>& parents,
                  const std::vector<std::vector<int>>& children,
                  std::vector<int>& levels, int node)
{
    if (parents[node].empty()) levels[node] = 0;
    else {
        std::vector<int> parents_level;
        for (std::vector<int>::size_type j = 0; j < parents[node].size(); ++j)
        {
            parents_level.push_back(levels[parents[node][j]]);
        }

        levels[node] = (*max_element(parents_level.begin(), parents_level.end()))+1;
    }

    for (std::vector<int>::size_type j = 0; j < children[node].size(); ++j)
    {
        assign_level(parents, children, levels, children[node][j]);
    }
}

/*******************************************************************************
 * This approach is basically a depth-first search traversal starting from each
 * node that has no incoming links. When a node has more than one parent, we
 * consider the greatest level among them in order to compute current node's
 * level.
 *
 *
 * @param A     Adjacency matrix.
 * @param P     Permutation matrix.
 ******************************************************************************/
void compute_levels(const std::vector<std::vector<int>>& parents,
                    const std::vector<std::vector<int>>& children,
                    std::vector<int>& levels)
{
    for (std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if (parents[i].empty()) assign_level(parents, children, levels, i);
    }
}

/*******************************************************************************
 * Performs a level-based approach considering node's positions
 *
 *
 * @param A     Adjacency matrix.
 * @param P     Permutation matrix.
 ******************************************************************************/
void apply_levels(const Eigen::MatrixXf& A, Eigen::MatrixXf& P)
{
    size_t size = A.rows(); // Number of nodes

    std::vector<std::vector<int>> parents(size), children(size);

    parents_from_matrix(A, parents);
    children_from_matrix(A, children);

    std::vector<int> levels(size, -1);

    compute_levels(parents, children, levels);

    std::vector<std::pair<int, int>> sorted_levels = sort_levels(levels);

    label_nodes(A, P, sorted_levels);
}
