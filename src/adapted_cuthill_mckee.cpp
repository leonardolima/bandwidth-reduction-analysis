#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <Eigen/Dense>
#include "adapted_cuthill_mckee.h"
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
void nodal_numbering(const std::vector<std::vector<int>>& children,
                     std::vector<int>& path,
                     std::vector<int>& indegree,
                     std::vector<int>& outdegree,
                     std::vector<bool>& marked,
                     int num_nodes,
                     int starting_node)
{
    std::deque<int> d;
    int seen = 0;

    // Disconnected nodes
    for(int i = 0; i < num_nodes; ++i)
    {
        if(indegree[i] == 0 && outdegree[i] == 0 && !marked[i])
        {
            path.push_back(i);
            marked[i] = true;
            seen++;
        }
    }

    path.push_back(starting_node);
    marked[starting_node] = true;
    seen++;

    // Connected nodes
    while(seen < num_nodes)
    {
        for(int i = 0; i < num_nodes; ++i)
        {
            if(indegree[i] == 0 && !marked[i])
            {
                for(std::vector<int>::size_type j = 0; j < children[i].size(); ++j)
                {
                    path.push_back(children[i][j]);
                    indegree[children[i][j]]--;
                    marked[children[i][j]] = true;
                    seen++;
                }

                path.push_back(i);
                marked[i] = true;
            }

        }
    }
}

std::vector<int> select_starting_nodes(const Eigen::MatrixXf& A,
                                       const std::vector<std::vector<int>>& parents)
{
    std::vector<int> starting_nodes;

    for(std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if(parents[i].size() == 0) starting_nodes.push_back(i);
    }

    return starting_nodes;
}

void sort_children(const Eigen::MatrixXf& A, std::vector<std::vector<int>>& children)
{
    std::vector<std::pair<int, int>> pairs;

    std::vector<int> nodes_deg = compute_nodes_deg(A);

    for(std::vector<std::vector<int>>::size_type i = 0; i < children.size(); ++i)
    {
        for(std::vector<int>::size_type j = 0; j < children[i].size(); ++j)
        {
            std::pair<int, int> pair = std::make_pair(nodes_deg[children[i][j]], children[i][j]);
            pairs.push_back(pair);
        }

        std::sort(pairs.begin(), pairs.end());
        children[i].clear();

        for(std::vector<std::pair<int, int>>::size_type j = pairs.size(); j >= 0; --j)
        {
            children[i].push_back(pairs[j].second);
        }

        pairs.clear();
    }
}

/*******************************************************************************
 * Performs an adapted topological sorting to the graph
 *
 *
 * @param A     Adjacency matrix.
 * @param P     Permutation matrix.
 ******************************************************************************/
void apply_adapted_cuthill_mckee(const Eigen::MatrixXf& A, Eigen::MatrixXf& P)
{
    int num_nodes = A.rows(); // Matrix dimension

    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(num_nodes, num_nodes);
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(num_nodes, num_nodes);

    std::vector<std::vector<int>> children(num_nodes), parents(num_nodes), possible_paths;
    std::vector<int> indegree(num_nodes), outdegree(num_nodes), starting_nodes;
    std::vector<bool> marked(num_nodes);
    std::vector<int> path;

    // int min_bandwidth = std::numeric_limits<int>::max();
    int min_peak = std::numeric_limits<int>::max();

    children_from_matrix(A, children);
    parents_from_matrix(A, parents);

    for(std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        indegree[i] = parents[i].size();
        outdegree[i] = children[i].size();
    }

    sort_children(A, children);

    starting_nodes = select_starting_nodes(A, parents);

    for(std::vector<int>::size_type i = 0; i < starting_nodes.size(); ++i)
    {
        nodal_numbering(children, path, indegree, outdegree, marked, num_nodes, starting_nodes[i]);
        possible_paths.push_back(path);

        int peak = peak_mem(A, path);

        for(std::vector<int>::size_type j = 0; j < path.size(); ++j)
        {
            std::cout << path[j] << " ";
        }
        std::cout << std::endl;

        std::cout << "peak = " << peak << std::endl;
        // label_sorted_nodes(A, P, make_sorted_pairs(possible_paths[i]));
        // M = (P*A*P.transpose());
        // if(peak < min_peak)
        // {
        //     R = P;
        //     min_peak = peak;
        // }
        // P.setZero();
    }
    // std::cout << "min_peak = " << min_peak << std::endl;
    // P = R;
}
