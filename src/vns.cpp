#include <iostream>
#include <random>
#include <vector>
#include <limits>
#include <Eigen/Dense>
#include "basic.h"

// TODO: Change comments to Doxygen format

/* Function: max_adj_label
 *
 * Inputs:
 *         v        - vertex
 *         adj_list - adjacency list representation of the graph G(M)
 *         label    - labeling function
 *
 */
int max_adj_label (int v, const std::vector<std::vector<int>>& adj_list, const std::vector<int> label)
{
    int max = 0;

    for (int i = 0; i < adj_list[v].size(); ++i)
    {
        if (label[adj_list[v][i]] > max) max = label[adj_list[v][i]];
    }

    return max;
}

/* Function: min_adj_label
 *
 * Inputs:
 *         v        - vertex
 *         adj_list - adjacency list representation of the graph G(M)
 *         label    - labeling function
 *
 */
int min_adj_label (int v, const std::vector<std::vector<int>>& adj_list, const std::vector<int> label)
{
    int min = std::numeric_limits<int>::max();

    for (int i = 0; i < adj_list[v].size(); ++i)
    {
        if (label[adj_list[v][i]] < min) min = label[adj_list[v][i]];
    }

    return min;
}

/* Function: critical_v
 *
 * Inputs:
 *         v        - vertex
 *         adj_list - adjacency list representation of the graph G(M)
 *         label    - labeling function
 *
 * Returns critical vertex index, i.e., considering a vertex u, it is
 * the one which satisfies |f(u) - f(v)| = B(f, u)
 */
int critical_v (int v, const std::vector<std::vector<int>>& adj_list, const std::vector<int> label)
{
    int max = 0, crit_v;

    for (int i = 0; i < adj_list[v].size(); ++i)
    {
        if (label[adj_list[v][i]] - label[v] > max)
        {
            max = label[adj_list[v][i]] - label[v];
            crit_v = adj_list[v][i];
        }
    }

    return crit_v;
}

/* Function: v_bandwidth
 *
 * Inputs:
 *         v        - vertex
 *         adj_list - adjacency list representation of the graph G(M)
 *         label    - labeling function
 *
 * Computes vertex bandwidth considering a particular labeling function
 */
int v_bandwidth (int v, const std::vector<std::vector<int>>& adj_list, const std::vector<int> label)
{
    int max = 0;

    for (int i = 0; i < adj_list[v].size(); ++i)
    {
        if (label[adj_list[v][i]] - label[v] > max) max = label[adj_list[v][i]] - label[v];
    }

    return max;
}

/* Function: shake_1
 *
 * Inputs:
 *         adj_list - adjacency list representation of the graph G(M)
 *         k        - decimal value representing kth neighborhood of labeling
 *         label    - labeling function
 *         rng      - random number generator (mt19937)
 *
 */
void shake_1 (const std::vector<std::vector<int>>& adj_list, int k, std::vector<int> label,
              std::mt19937 rng)
{
    int N = adj_list.size();

    std::vector<std::pair<int, int>> v_b;
    std::vector<std::pair<int, int>> K(k);

    int min = std::numeric_limits<int>::max();

    std::uniform_int_distribution<std::mt19937::result_type> uni_dist_k(0, k);

    // Compute v_bandwidth for each vertex
    for (int i = 0; i < N; ++i)
    {
        std::pair<int, int> pair = std::make_pair(v_bandwidth(i, adj_list, label), i);
        v_b.push_back(pair);
    }

    // Sort and select k vertices that have the greatest bandwidths
    std::sort(v_b.begin(), v_b.end());
    for (int i = N-k-1; i < N; ++i) K[i] = v_b[i];

    for (int i = 0; i < k; i++)
    {
        int rnd = uni_dist_k(rng);
        int u = K[rnd].second, v;

        // Find critical vertex of u, i.e., |f(u) - f(v)| = B(f, u)
        v = critical_v(u, adj_list, label);

        // Find w s.t. max{f(v) - f_min(w), f_max(w) - f(v)} is minimum
        // where f_min(u) <= f(w) <= f_max(u)
        int w = 0;
        for (int j = 0; j != u && j < N; ++j)
        {
            if (label[j] <= max_adj_label(u) && label[j] >= min_adj_label(u))
            {
                int max = 0;

                if (label[v] - min_adj_label(w) > max_adj_label(w) - label[v])
                {
                    max = label[v] - min_adj_label(w);
                } else {
                    max = max_adj_label(w) - label[v];
                }

                if (min > max)
                {
                    min = max;
                    w = j;
                }
            }
        }

        std::swap(label[v], label[w]);
    }
}

/* Function: initial_solution
 *
 * Inputs:
 *         adj_list - adjacency list representation of the graph G(M)
 *         label    - labeling function
 *         rng      - random number generator (mt19937)
 *
 */
void initial_solution (const std::vector<std::vector<int>>& adj_list, std::vector<int> label,
                       std::mt19937 rng)
{
    int N = adj_list.size(); // N = |V|

    std::uniform_int_distribution<std::mt19937::result_type> uni_dist_N(0, N);

    std::vector<bool> mark(N);
    std::vector<int> q, s;
    int ql; // cardinality of q vector
    int l;  // loop variable

    // 1. Initially, assign false to mark[i], i \in (1, ..., N)
    for (int i = 0; i < N; i++)
    {
        mark[i] = false;
    }

    // 2. Pick randomly starting vertex
    int k = uni_dist_N(rng);
    mark[k] = true;
    q.push_back(k);
    ql = 1;
    label[k] = 0;
    l = 1;

    // 3.
    while (l < N)
    {
        std::uniform_int_distribution<std::mt19937::result_type> uni_dist_ql(0, ql);

        s.clear();

        for (int i1 = 0; i1 < ql; ++i1)
        {
            int i = q[i1];
            for (int j1 = 0; j1 < adj_list[i].size(); ++j1)
            {
                int j = adj_list[i][j1];

                if (!mark[j])
                {
                    s.push_back(j);
                    mark[j] = true;
                }
            }
        }

        int j_star = uni_dist_ql(rng);

        for (int j1 = 0; j1 < ql; ++j1)
        {
            int j = q[j_star];
            l++;
            label[j] = l;
            j_star++;

            if (j_star >= ql-1) j_star = 0;
        }

        for (int i = 0; i < s.size(); i++) q.push_back(s[i]);
        ql = s.size();
    }
}

/* Function: construct_adj_list
 *
 * Inputs:
 *         M              - original matrix representing the graph
 *
 * When considering just the upper triangular part of the matrix M,
 * ordering is preserved, hence constructing adj_list becomes straightforward
 */
void construct_adj_list (const Eigen::MatrixXf& M, std::vector<std::vector<int>>& adj_list)
{
    for (int i = 0; i < M.rows(); ++i)
    {
        for(int j = 0; j < M.cols(); ++j)
        {
            if (M(i,j) != 0 && i != j) adj_list[i].push_back(j);
        }
    }
}

/* Function: compute_matrices
 *
 * Inputs:
 *         M              - original matrix representing the graph
 *         P              - permutation matrix
 *         R              - lowest degree resulting matrix (i.e. P*A*P^T)
 *         O              - vector corresponding to the constraint ordering L_x < [L_yi]
 *
 * Perform VNS metaheuristic considering ordering constraint
 */
void apply_vns (const Eigen::MatrixXf& M, Eigen::MatrixXf& P, Eigen::MatrixXf& R,
                       const std::vector<std::vector<int>>& O)
{
    int N = M.rows();

    // Initializing random number generator
    std::random_device dev;
    std::mt19937 rng(dev());

    std::vector<std::vector<int>> adj_list(N);

    // P(0, 0) = 1;                          // Starting node = 0 is labeled as 0
    // nodal_numbering(A, P, rows_deg, O);   // Execute algorithm on A
    // R = (P*A*P.transpose());              // Compute resulting matrix R

    // std::cout << "dim(A) = " << A.rows() << "x" << A.cols() << std::endl;
    // std::cout << "bandwidth(A) = " << matrix_bandwidth(A) << std::endl;
    // std::cout << "bandwidth(R) = " << matrix_bandwidth(R) << std::endl;
}
