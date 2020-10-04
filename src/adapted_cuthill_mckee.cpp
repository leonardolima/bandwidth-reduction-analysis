#include <iostream>
#include <vector>
#include <deque>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>
#include "bandwidth_analysis/cuthill_mckee.h"
#include "adapted_cuthill_mckee.h"
#include "basic.h"
#include "io.h"

/*******************************************************************************
 * Check if constraints are satisfied, i.e., in order to label the current
 * node, it is necessary to check if L_x < [L_yi] for every i.
 *
 *
 * @param node  Node's index.
 * @param label Possible label at this step for this particular node.
 * @param O     Ordering constraints (in the form L_x < [L_yi]).
 * @param P     Permutation matrix.
 ******************************************************************************/
bool check_constraints (int node, int label, Eigen::MatrixXf& P,
                        const std::vector<std::vector<int>>& prec_O,
                        const std::vector<std::vector<int>>& succ_O)
{
    bool satisfy_ordering = true, have_labelled_parents = true;

    for (std::vector<int>::size_type j = 0; j < succ_O[node].size(); ++j)
    {
        int adj_node_label = node_label(P, succ_O[node][j]);
        if (adj_node_label != -1 && adj_node_label < label) satisfy_ordering = false;
    }

    for (std::vector<int>::size_type j = 0; j < prec_O[node].size(); ++j)
    {
        int adj_node_label = node_label(P, prec_O[node][j]);
        if (adj_node_label == -1) have_labelled_parents = false;
    }

    return satisfy_ordering && have_labelled_parents;
}

void apply_label_to_node (Eigen::MatrixXf& P,
                          const std::vector<std::vector<int>>& prec_O,
                          const std::vector<std::vector<int>>& succ_O,
                          std::deque<int>& d, int node, int& new_label,
                          std::vector<int>& visited_nodes,
                          std::vector<bool>& node_in_d)
{
    if (!node_has_label(P, node))
    {
        if (check_constraints(node, new_label, P, prec_O, succ_O))
        {
            P(new_label++, node) = 1;
            visited_nodes[node] = 1;
        } else {
            d.push_back(node);
            node_in_d[node] = true;
        }
    }
}

/*******************************************************************************
 * Adapted version of the original algorithm in order to consider the constraints
 * imposed. When labeling the adjacent nodes, it is necessary to check if their
 * label respect the constraints or not. If so, we label them accordingly.
 * Otherwise we keep a queue with all the remaining nodes that haven't been
 * labeled already, checking again at each step (i.e., when considering
 * other nodes if labeling is possible).
 *
 *
 * @param A             Adjacency matrix.
 * @param P             Permutation matrix.
 * @param rows_deg      Vector with the degree of each node's rows.
 * @param O             Ordering constraints (in the form L_x < [L_yi]).
 * @param initial_nodes List of nodes of each subgraph of the graph.
 ******************************************************************************/
void adapted_nodal_numbering(const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                             const std::vector<int>& nodes_deg,
                             const std::vector<std::vector<int>>& prec_O,
                             const std::vector<std::vector<int>>& succ_O,
                             const std::vector<int>& initial_nodes)
{
    int N = A.cols();
    int new_label = 0;

    std::deque<int> d;
    std::vector<bool> node_in_d(N, false); // Existence of node in d
    std::vector<int> visited_nodes(N, 0);

    for (std::vector<int>::size_type i = 0; i < initial_nodes.size(); ++i)
    {
        // Mark initial_nodes[i] as visited
        visited_nodes[initial_nodes[i]] = 1;

        // Push initial_nodes[i] to d
        d.push_back(initial_nodes[i]);
        node_in_d[initial_nodes[i]] = true;

        while (!d.empty())
        {
            int u = d.front();
            d.pop_front();
            node_in_d[u] = false;

            apply_label_to_node(P, prec_O, succ_O, d, u, new_label, visited_nodes, node_in_d);

            for (std::deque<int>::size_type i = 0; i < d.size(); ++i)
            {
                int n = d.front();
                d.pop_front();
                node_in_d[n] = false;

                apply_label_to_node(P, prec_O, succ_O, d, n, new_label, visited_nodes, node_in_d);
            }

            if (nodes_deg[u] > 0)
            {
                Eigen::VectorXf adj_nodes = A.row(u);
                std::vector<std::pair<int, int>> sorted_adj_nodes = sort_adj_nodes(adj_nodes, nodes_deg, P);

                for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_adj_nodes.size(); ++i)
                {
                    int cur_node = sorted_adj_nodes[i].second;
                    if (visited_nodes[cur_node] == 0 && check_constraints(cur_node, new_label, P, prec_O, succ_O))
                    {
                        P(new_label++, cur_node) = 1;
                        visited_nodes[cur_node] = 1;
                    }

                    if (node_in_d[cur_node] == false)
                    {
                        d.push_back(cur_node);
                        node_in_d[cur_node] = true;
                    }
                }

                if (node_in_d[u] == false) visited_nodes[u] = 2;
            }
        }
    }
}

std::vector<int> select_initial_nodes (const Eigen::MatrixXf& A)
{
    int N = A.rows();

    std::vector<int> initial_nodes;
    std::vector<std::vector<int>> parents(N);

    parents_from_matrix(A, parents);

    // Find all source nodes
    for (std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    {
        if (parents[i].size() == 0) initial_nodes.push_back(i);
    }

    return initial_nodes;
}

/*******************************************************************************
 * Performs Cuthill-McKee algorithm considering ordering constraint.
 *
 *
 * @param A     Adjacency matrix.
 * @param A_sym Symmetric adjacency matrix.
 * @param P     Permutation matrix.
 * @param R     Resulting matrix (i.e., P*A*P^T).
 * @param O     Ordering constraints (in the form L_x < [L_yi]).
 ******************************************************************************/
void apply_adapted_cuthill_mckee(const Eigen::MatrixXf& A, const Eigen::MatrixXf& A_sym,
                                 Eigen::MatrixXf& P, Eigen::MatrixXf& R,
                                 const std::vector<std::vector<int>>& prec_O,
                                 const std::vector<std::vector<int>>& succ_O)
{
    // std::cout << A << std::endl;

    std::vector<int> nodes_deg = compute_nodes_deg(A_sym);

    // Starting node of each subgraph of A
    std::vector<int> initial_nodes = select_initial_nodes(A);

    adapted_nodal_numbering(A_sym, P, nodes_deg, prec_O, succ_O, initial_nodes);
}
