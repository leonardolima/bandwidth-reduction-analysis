#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <Eigen/Dense>
#include "adapted_cuthill_mckee.h"
#include "peak_mem.h"
#include "basic.h"
#include "io.h"

std::vector<int> list_disconnected_nodes(const std::vector<std::vector<int>>& parents,
                                         const std::vector<std::vector<int>>& children)
{
    std::vector<int> disconnected_nodes;

    // Disconnected nodes
    for(std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if(parents[i].empty() && children[i].empty()) disconnected_nodes.push_back(i);
    }

    return disconnected_nodes;
}

void process_node(const std::vector<std::vector<int>>& children,
                  std::vector<int>& path,
                  std::vector<int>& indegree,
                  std::vector<bool>& marked,
                  int& seen,
                  int num_nodes,
                  int node)
{
    if(indegree[node] == 0 && !marked[node])
    {
        // std::cout << "processing: " << node << std::endl;
        for(std::vector<int>::size_type j = 0; j < children[node].size(); ++j)
        {
            indegree[children[node][j]]--;
        }
        path.push_back(node);
        marked[node] = true;
        seen++;
    }
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
void nodal_numbering(const std::vector<std::vector<int>>& parents,
                     const std::vector<std::vector<int>>& children,
                     const std::vector<int>& starting_nodes,
                     std::vector<int>& path,
                     std::vector<bool>& marked,
                     std::vector<int> indegree, // Copy
                     int num_nodes)
{
    int seen = 0, counter = 0;
    std::vector<int> disconnected_nodes = list_disconnected_nodes(parents, children);

    // Process disconnected nodes
    for(std::vector<int>::size_type i = 0; i < disconnected_nodes.size(); ++i)
    {
        process_node(children, path, indegree, marked, seen, num_nodes, i);
    }

    // Connected nodes
    while(seen < num_nodes)
    {
        // Starting node
        // std::cout << "starting_node = " << starting_nodes[counter] << std::endl;
        process_node(children, path, indegree, marked, seen, num_nodes, starting_nodes[counter]);
        counter++;

        // std::cout << "seen = " << seen << std::endl;
        for(int i = 0; i < num_nodes; ++i)
        {
            // std::cout << "i = " << i << std::endl;
            process_node(children, path, indegree, marked, seen, num_nodes, i);
        }
    }
}

std::vector<int> select_starting_nodes(const std::vector<std::vector<int>>& parents,
                                       const std::vector<std::vector<int>>& children)
{
    std::vector<int> starting_nodes;

    // Starting nodes in each component of the graph
    for(std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if(parents[i].empty() && !children[i].empty()) starting_nodes.push_back(i);
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

        ////////////////////////////////////////////////////////
        // std::cout << "Pairs: ";
        // for(std::vector<int>::size_type j = 0; j < pairs.size(); ++j)
        // {
        //     std::cout << "(" << pairs[j].first << ", " << pairs[j].second << ") | ";
        // }
        // std::cout << std::endl;
        ////////////////////////////////////////////////////////

        std::sort(pairs.begin(), pairs.end(), std::greater<>());
        children[i].clear();

        for(std::vector<std::pair<int, int>>::size_type j = 0; j < pairs.size(); ++j)
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

    std::vector<std::vector<int>> children(num_nodes), parents(num_nodes);
    std::vector<int> indegree(num_nodes), outdegree(num_nodes), starting_nodes;
    std::vector<bool> marked(num_nodes, false);
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
    starting_nodes = select_starting_nodes(parents, children);

    // for(std::vector<int>::size_type j = 0; j < starting_nodes.size(); ++j)
    // {
    //     std::cout << starting_nodes[j] << " ";
    // }
    // std::cout << std::endl;

    do {
        nodal_numbering(parents, children, starting_nodes, path, marked, indegree, num_nodes);
        int peak = peak_mem(A, path);

        // for(std::vector<int>::size_type j = 0; j < path.size(); ++j)
        // {
        //     std::cout << path[j] << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "peak = " << peak << std::endl;

        label_sorted_nodes(A, P, make_sorted_pairs(path));
        M = (P*A*P.transpose());
        if(peak < min_peak)
        {
            R = P;
            min_peak = peak;
        }
        P.setZero();

        // Clearing/resetting auxiliary vectors
        std::fill(marked.begin(), marked.end(), false);
        path.clear();
    } while(std::next_permutation(starting_nodes.begin(), starting_nodes.end()));

    std::cout << "min_peak = " << min_peak << std::endl;
    P = R;
}
