#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "topological.h"
#include "peak_mem.h"
#include "basic.h"
#include "io.h"

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
    for(size_t i = 0; i < num_nodes; ++i)
    {
        if(indegree[i] == 0 && !marked[i])
        {
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
void apply_topological(const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                       std::vector<int>& worst_path, std::vector<int>& chosen_path)
{
    const size_t num_nodes = A.rows(); // Matrix dimension

    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(num_nodes, num_nodes);

    std::vector<std::vector<int>> children(num_nodes), parents(num_nodes), possible_paths;
    std::vector<int> indegree(num_nodes), path;
    std::vector<bool> marked(num_nodes);

    int min_peak = std::numeric_limits<int>::max();
    int max_peak = 0;

    children_from_matrix(A, children);
    parents_from_matrix(A, parents);

    for(std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        indegree[i] = parents[i].size();
    }

    topological_sorting(children, possible_paths, path, indegree, marked, num_nodes);

    for(std::vector<std::vector<int>>::size_type i = 0; i < possible_paths.size(); ++i)
    {
        int peak = peak_mem(A, possible_paths[i]);
        label_sorted_nodes(A, R, make_sorted_pairs(possible_paths[i]));

        // Maximum peak memory usage
        if(peak > max_peak)
        {
            max_peak = peak;
            worst_path = possible_paths[i];
        }

        // Minimum peak memory usage
        if(peak < min_peak)
        {
            P = R;
            min_peak = peak;
            chosen_path = possible_paths[i];
        }

        R.setZero();
    }
}
