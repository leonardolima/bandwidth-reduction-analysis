#include <iostream>
#include <vector>
#include <deque>
#include <Eigen/Dense>
#include "peak_mem.h"
#include "basic.h"

void update_peak(int cur_mem, int& peak)
{
    if (cur_mem > peak) peak = cur_mem;
}

void free_parents(const std::vector<std::vector<int>>& parents,
                  const std::vector<int>& outdegree,
                  const std::vector<int>& mem,
                  int& cur_mem, int node)
{
    // Free nodes that don't have children yet to be processed
    for(std::vector<int>::size_type j = 0; j < parents[node].size(); ++j)
    {
        if (outdegree[parents[node][j]] == 0)
        {
            // std::cout << "Freeing: " << parents[node][j] << std::endl;
            cur_mem -= mem[parents[node][j]];
        }
    }
}

bool process_node(const std::vector<std::vector<int>>& children,
                  const std::vector<std::vector<int>>& parents,
                  const std::vector<int>& mem, std::vector<int>& indegree,
                  std::vector<int>& outdegree, int& peak, int& cur_mem, int node)
{
    // 1. Process a node if it can be processed, i.e., if all parents have already
    //    been processed
    if (indegree[node] == 0)
    {
        cur_mem += mem[node];
        update_peak(cur_mem, peak);

        // std::cout << "Processing: " << node << std::endl;
        // std::cout << "cur_mem = " << cur_mem << std::endl;
        // std::cout << "peak = " << peak << std::endl;

        for(std::vector<int>::size_type j = 0; j < children[node].size(); ++j)
        {
            indegree[children[node][j]]--;
        }

        for(std::vector<int>::size_type j = 0; j < parents[node].size(); ++j)
        {
            outdegree[parents[node][j]]--;
        }

        free_parents(parents, outdegree, mem, cur_mem, node);
        return true;
    } else {
        free_parents(parents, outdegree, mem, cur_mem, node);
        return false;
    }
}

int compute_peak_mem(const std::vector<std::vector<int>>& children,
                     const std::vector<std::vector<int>>& parents,
                     const std::vector<int>& mem, const std::vector<int>& path,
                     std::vector<int>& indegree, std::vector<int>& outdegree)
{
    int cur_mem = 0, peak = 0;
    std::deque<int> to_free;

    for(std::vector<int>::size_type i = 0; i < path.size(); ++i)
    {
        for (std::deque<int>::size_type j = 0; j < to_free.size(); ++j)
        {
            int node = to_free.front();
            if (process_node(children, parents, mem, indegree, outdegree, peak, cur_mem, node))
            {
                to_free.pop_front();
            }
        }
        if (!process_node(children, parents, mem, indegree, outdegree, peak, cur_mem, path[i]))
        {
            to_free.push_back(path[i]);
        }
    }
    return peak;
}

/*******************************************************************************
 * Considering a directed acyclic graph corresponding to a computer program and
 * a particular path, computes its peak memory usage
 *
 *
 * @param A matrix representation of the graph
 * @param path Contains a possible path in the graph
 ******************************************************************************/
int peak_mem(const Eigen::MatrixXf& A, const std::vector<int>& path)
{
    const size_t num_nodes = A.rows();

    std::vector<std::vector<int>> parents(num_nodes), children(num_nodes);
    std::vector<int> mem(num_nodes), indegree(num_nodes), outdegree(num_nodes);

    parents_from_matrix(A, parents);
    children_from_matrix(A, children);
    mem_from_matrix(A, mem);

    for(std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        indegree[i] = parents[i].size();
        outdegree[i] = children[i].size();
    }

    return compute_peak_mem(children, parents, mem, path, indegree, outdegree);
}
