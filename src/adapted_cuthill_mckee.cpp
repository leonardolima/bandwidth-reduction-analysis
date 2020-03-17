#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>
#include "bandwidth_analysis/cuthill_mckee.h"
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

void apply_label_to_node (Eigen::MatrixXf& P, const std::vector<std::vector<int>>& prec_O,
                          const std::vector<std::vector<int>>& succ_O,
                          std::queue<int>& q, int node, int& new_label)
{
    if (check_constraints(node, new_label, P, prec_O, succ_O)) {
        P(new_label++, node) = 1;
        q.pop();
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
void adapted_nodal_numbering (const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                              const std::vector<int>& rows_deg,
                              const std::vector<std::vector<int>>& prec_O,
                              const std::vector<std::vector<int>>& succ_O,
                              const std::vector<int>& initial_nodes)
{
    int N = A.cols();
    // std::cout << "N = " << N << std::endl;
    int new_label = 0;

    std::queue<int> q;
    std::vector<bool> node_in_q(N, false); // Existence of node in q
    std::vector<int> visited_nodes(N, 0);

    for(std::vector<int>::size_type i = 0; i < initial_nodes.size(); ++i)
    {
        // Mark initial_nodes[i] as visited
        visited_nodes[initial_nodes[i]] = 1;

        // Push initial_nodes[i] to q
        q.push(initial_nodes[i]);
        node_in_q[initial_nodes[i]] = true;

        while (!q.empty())
        {
            // std::cout << "--------------------------" << std::endl;
            // print_queue(q);
            // std::cout << std::endl;

            int u = q.front();
            q.pop();
            node_in_q[u] = false;

            // std::cout << "select:" << u << std::endl;

            if (!node_has_label(P, u))
            {
                if (check_constraints(u, new_label, P, prec_O, succ_O))
                {
                    // std::cout << "new_label = " << new_label << "| cur_node = " << u << std::endl;
                    P(new_label++, u) = 1;
                    visited_nodes[u] = 1;
                } else {
                    q.push(u);
                    node_in_q[u] = true;
                }
            }

            // print_queue(q);
            // std::cout << std::endl;

            if (rows_deg[u] > 0)
            {
                Eigen::VectorXf node_row = A.row(u);
                std::vector<std::pair<int, int>> sorted_rows_deg = sort_rows_deg(node_row, rows_deg, P);

                for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_rows_deg.size(); ++i)
                {
                    int cur_node = sorted_rows_deg[i].second;
                    if (visited_nodes[cur_node] == 0 && check_constraints(cur_node, new_label, P, prec_O, succ_O))
                    {
                        // std::cout << "new_label = " << new_label << "| cur_node = " << cur_node << std::endl;
                        P(new_label++, cur_node) = 1;
                        visited_nodes[cur_node] = 1;
                    }

                    if (node_in_q[cur_node] == false)
                    {
                        q.push(cur_node);
                        node_in_q[cur_node] = true;
                    }
                }

                if (node_in_q[u] == false) visited_nodes[u] = 2;
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
void apply_adapted_cuthill_mckee (const Eigen::MatrixXf& A, const Eigen::MatrixXf& A_sym,
                                  Eigen::MatrixXf& P, Eigen::MatrixXf& R,
                                  const std::vector<std::vector<int>>& prec_O,
                                  const std::vector<std::vector<int>>& succ_O)
{
    // std::cout << A << std::endl;

    std::vector<int> rows_deg = compute_rows_deg(A_sym);

    // Starting node of each subgraph of A
    std::vector<int> initial_nodes = select_initial_nodes(A);

    adapted_nodal_numbering(A_sym, P, rows_deg, prec_O, succ_O, initial_nodes);
}
