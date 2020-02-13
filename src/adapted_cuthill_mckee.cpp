#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>
#include "basic.h"
#include "bandwidth_analysis/cuthill_mckee.h"

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
bool check_constraints (int node, int label, const std::vector<std::vector<int>>& O,
                        Eigen::MatrixXf& P)
{
    // std::cout << "Checking node = " << node << std::endl;

    // for (std::vector<std::vector<int>>::size_type i; i < O.size(); ++i)
    // {
    //     std::cout << "O[" << i << "] = ";
    //     for (std::vector<int>::size_type j; j < O[i].size(); ++j) std::cout << O[i][j] << " ";
    //     std::cout << std::endl;
    // }

    for (std::vector<int>::size_type j = 0; j < O[node].size(); ++j)
    {
        int adj_node_label = node_label(P, O[node][j]);
        if (adj_node_label != -1 && adj_node_label < label) return false;
    }

    return true;
}

void apply_label_to_node (Eigen::MatrixXf& P, const std::vector<std::vector<int>>& O, std::queue<int>& q,
                          int node, int& new_label)
{
    if (check_constraints(node, new_label, O, P)) {
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
 * @param A        Adjacency matrix.
 * @param P        Permutation matrix.
 * @param rows_deg Vector with the degree of each node's rows.
 * @param O        Ordering constraints (in the form L_x < [L_yi]).
 ******************************************************************************/
void adapted_nodal_numbering (const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                              const std::vector<int>& rows_deg,
                              const std::vector<std::vector<int>>& O)
{
    int new_label = 1;


    std::queue<int> q;

    // 3. First we need to label unconnected nodes
    for(int j = 0; j < A.cols(); ++j)
    {
        if (rows_deg[j] == 0) P(new_label++, j) = 1; // Update matrix P accordingly
    }

    // 4. For every node we (i) list the adjacent nodes (ii) sort them in ascending order
    // of degree (iii) label them uniquely and accordingly (2, 3, ...) considering the
    // ordering O
    for(int j = 0; j < A.cols(); ++j)
    {
        // Apply label to first element on the queue
        if (!q.empty()) apply_label_to_node(P, O, q, q.front(), new_label);

        if (rows_deg[j] > 0)
        {
            int node_row_index = node_index(P, j);
            Eigen::VectorXf node_row = A.row(node_row_index);
            std::vector<std::pair<int, int>> sorted_rows_deg = sort_rows_deg(node_row, rows_deg, P);

            for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_rows_deg.size(); ++i)
            {
                int cur_node = sorted_rows_deg[i].second;
                if (check_constraints(cur_node, new_label, O, P)) P(new_label++, cur_node) = 1;
                else q.push(cur_node);
            }
        }
    }
}

/*******************************************************************************
 * Performs Cuthill-McKee algorithm considering ordering constraint.
 *
 *
 * @param A Adjacency matrix.
 * @param P Permutation matrix.
 * @param R Resulting matrix (i.e., P*A*P^T).
 * @param O Ordering constraints (in the form L_x < [L_yi]).
 ******************************************************************************/
void apply_adapted_cuthill_mckee (const Eigen::MatrixXf& A, Eigen::MatrixXf& P,
                                  Eigen::MatrixXf& R, const std::vector<std::vector<int>>& O)
{
    std::vector<int> rows_deg = compute_rows_deg(A);

    P(0, 0) = 1;                                  // Starting node = 0 is labeled as 0
    adapted_nodal_numbering(A, P, rows_deg, O);   // Execute algorithm on A
    R = (P*A*P.transpose());                      // Compute resulting matrix R

    std::cout << "dim(A) = " << A.rows() << "x" << A.cols() << std::endl;
    std::cout << "bandwidth(A) = " << matrix_bandwidth(A) << std::endl;
    std::cout << "bandwidth(R) = " << matrix_bandwidth(R) << std::endl;
}
