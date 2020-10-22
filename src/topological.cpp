#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "topological.h"
#include "basic.h"
#include "io.h"

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

    // for (std::vector<int>::size_type i = 0; i < sorted_pairs.size(); ++i)
    // {
    //     std::cout << "sorted_pairs[" << i << "] = (" << sorted_pairs[i].first << ", " << sorted_pairs[i].second << ")" << std::endl;
    // }

    return sorted_pairs;
}

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

    for (std::vector<std::vector<int>>::size_type i = 0; i < children[node].size(); ++i)
    {
        visit(children, mark, sorted_nodes, children[node][i]);
    }

    mark[node] = 1;
    sorted_nodes.insert(sorted_nodes.begin(), node);
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
                         std::vector<int>& sorted_nodes,
                         int starting_node)
{
    std::vector<int> mark(children.size());

    visit(children, mark, sorted_nodes, starting_node); // If graph is connected this is enough
    std::vector<int>::iterator iter = std::find(mark.begin(), mark.end(), 0);

    while(iter != mark.end())
    {
        visit(children, mark, sorted_nodes, iter-mark.begin());

        iter = std::find(mark.begin(), mark.end(), 0);
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
    size_t A_dim = A.rows(); // Matrix dimension

    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(A_dim, A_dim);
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(A_dim, A_dim);

    std::vector<std::vector<int>> children(A_dim), parents(A_dim);

    std::vector<int> sorted_nodes, starting_nodes;

    int min_bandwidth = std::numeric_limits<int>::max();

    children_from_matrix(A, children);
    parents_from_matrix(A, parents);
    select_starting_nodes(parents, starting_nodes);

    // Check if std::vector<std::vector<int>>::size_type can be replaced with auto
    for (std::vector<std::vector<int>>::size_type i = 0; i < starting_nodes.size(); ++i)
    {
        topological_sorting(children, sorted_nodes, starting_nodes[i]);
        label_sorted_nodes(A, P, make_sorted_pairs(sorted_nodes));
        M = (P*A*P.transpose());
        int M_bandwidth = matrix_bandwidth(M);
        if (M_bandwidth < min_bandwidth)
        {
            R = P;
            min_bandwidth = M_bandwidth;
        }
        sorted_nodes.clear();
        P.setZero();
    }
    P = R;
}
