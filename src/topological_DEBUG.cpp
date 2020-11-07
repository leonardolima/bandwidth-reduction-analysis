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
                         std::vector<std::vector<int>>& possible_paths,
                         std::vector<int>& path,
                         std::vector<int>& indegree,
                         std::vector<bool>& marked,
                         const size_t num_nodes)
{
    // std::cout << "NEW LOOP" << std::endl;
    for(size_t i = 0; i < num_nodes; ++i)
    {
        // std::cout << "@ node = " << i << std::endl;
        if(indegree[i] == 0 && !marked[i])
        {
            // std::cout << "> Status: ";
            // std::cout << "indegree = ";
            // for(std::vector<int>::size_type j = 0; j < indegree.size(); ++j)
            // {
            //     std::cout << indegree[j] << " ";
            // }
            // std::cout << std::endl;

            // std::cout << "marked = ";
            // for(std::vector<int>::size_type j = 0; j < marked.size(); ++j)
            // {
            //     std::cout << marked[j] << " ";
            // }
            // std::cout << std::endl;

            // std::cout << "path = ";
            // for(std::vector<int>::size_type j = 0; j < path.size(); ++j)
            // {
            //     std::cout << path[j] << " ";
            // }
            // std::cout << std::endl;

            for(std::vector<int>::size_type j = 0; j < children[i].size(); ++j)
            {
                indegree[children[i][j]]--;
            }

            path.push_back(i);
            marked[i] = true;
            topological_sorting(children, possible_paths, path, indegree, marked, num_nodes);

            // Backtracking
            for(std::vector<int>::size_type j = 0; j < children[i].size(); ++j)
            {
                indegree[children[i][j]]++;
            }

            path.pop_back();
            marked[i] = false;
        }
    }

    if(path.size() == num_nodes) possible_paths.push_back(path);
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
    const size_t num_nodes = A.rows(); // Matrix dimension

    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(num_nodes, num_nodes);
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(num_nodes, num_nodes);

    std::vector<std::vector<int>> children(num_nodes), parents(num_nodes), possible_paths;
    std::vector<int> indegree(num_nodes), path;
    std::vector<bool> marked(num_nodes);

    int min_bandwidth = std::numeric_limits<int>::max();

    children_from_matrix(A, children);
    parents_from_matrix(A, parents);

    for(std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        indegree[i] = parents[i].size();
    }

    topological_sorting(children, possible_paths, path, indegree, marked, num_nodes);

    for(std::vector<std::vector<int>>::size_type i = 0; i < possible_paths.size(); ++i)
    {
        label_sorted_nodes(A, P, make_sorted_pairs(possible_paths[i]));
        M = (P*A*P.transpose());
        int M_bandwidth = matrix_bandwidth(M);
        if(M_bandwidth < min_bandwidth)
        {
            R = P;
            min_bandwidth = M_bandwidth;
        }
        path.clear();
        P.setZero();
    }
    P = R;

    std::cout << "Sorted nodes: ";
    for(std::vector<std::vector<int>>::size_type i = 0; i < possible_paths.size(); ++i)
    {
        for(size_t j = 0; j < num_nodes; ++j)
        {
            std::cout << possible_paths[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
