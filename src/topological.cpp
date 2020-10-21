#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "topological.h"
#include "basic.h"
#include "io.h"

void select_starting_nodes(const std::vector<std::vector<int>>& parents,
                           std::vector<int>& starting_nodes)
{
    for (std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if (parents[i].empty()) starting_nodes.push_back(i);
    }
}

void visit(const std::vector<std::vector<int>>& children,
           std::vector<int>& mark,
           std::vector<int>& sorted_nodes,
           int node)
{
    if(mark[node] == 1) return;

    for (std::vector<std::vector<int>>::size_type i = 0; i < children.size(); ++i)
    {
        visit(children, mark, sorted_nodes, i);
    }

    mark[node] == 1;
    sorted_nodes.push_back(node);
}

/*******************************************************************************
 * Implementation of a topological sorting algorithm, based on the depth-first
 * search algorithm. It goes through all nodes of the graph marking them until
 * it reaches a node that has already been visited. Since the graph is directed,
 * it can stop at that point.
 * mark[node] == 0 => unmarked node
 * mark[node] == 1 => marked node
 *
 *
 * @param children     List of children for each node
 * @param sorted_nodes Resulting vector with sorted nodes
 ******************************************************************************/
void topological_sorting(const std::vector<std::vector<int>>& children,
                         const std::vector<int>& starting_nodes,
                         std::vector<int>& sorted_nodes)
{
    std::vector<int> mark(size);

    // Check if std::vector<std::vector<int>>::size_type can be replaced with auto
    for (std::vector<std::vector<int>>::size_type i = 0; i < starting_nodes.size(); ++i)
    {
        visit(children, mark, sorted_nodes, i);
    }
}

/*******************************************************************************
 * Performs an adapted topological sorting to the graph
 *
 *
 * @param A     Adjacency matrix.
 * @param P     Permutation matrix.
 ******************************************************************************/
void apply_topological(const Eigen::MatrixXf& A, Eigen::MatrixXf& P)
{
    size_t size = A.rows(); // Number of nodes
    std::vector<std::vector<int>> children(size);
    std::vector<int> sorted_nodes(size), starting_nodes;

    children_from_matrix(A, children);

    select_starting_nodes(parents, starting_nodes);

    topological_sorting(children, starting_nodes, sorted_nodes);
}
