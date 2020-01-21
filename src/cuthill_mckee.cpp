/******************************
 *                            *
 *    Leonardo Lima, 2019     *
 *                            *
/******************************/

#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>
#include "basic.h"

/* Function: node_label
 *
 * Inputs:
 *         P   - permutation matrix
 *         col - column's index
 *
 */
int node_label (const Eigen::MatrixXf& P, unsigned int col)
{
    for (int i = 0; i < P.rows(); ++i)
    {
        if (P(i,col) != 0) return P(i, col);
    }

    return -1; // Node doesn't have a label
}

/* Function: sort_rows_deg
 *
 * Inputs:
 *         node_row - row in A of node
 *         rows_deg - nodes' degree
 *         P        - permutation matrix
 *
 */
std::vector<std::pair<int, int>> sort_rows_deg (const Eigen::VectorXf& node_row,
                                                const std::vector<int>& rows_deg,
                                                const Eigen::MatrixXf& P)
{
    std::vector<std::pair<int, int>> sorted_rows_deg;

    for(int i = 0; i < node_row.size(); ++i)
    {
        // We are not interested if the element is on the diagonal
        // or if the node already has a label
        if (node_row[i] != 0 && node_label(P, i) == -1)
        {
            std::pair<int, int> pair = std::make_pair(rows_deg[i], i);
            sorted_rows_deg.push_back(pair);
        }
    }

    std::sort(sorted_rows_deg.begin(), sorted_rows_deg.end());

    return sorted_rows_deg;
}

/* Function: node_index
 *
 * Inputs:
 *         P     - permutation matrix
 *         label - corresponding row of permutation matrix
 *
 * Returns index (column or row, it doesn't matter) in A (original matrix)
 *
 */
int node_index (const Eigen::MatrixXf& P, int label)
{
    // index in A corresponds to column (j)
    for(int j = 0; j < P.cols(); ++j)
    {
        if (P(label, j) != 0) return j;
    }

    return -1;
}

/* Function: check_constraints
 *
 * Inputs:
 *         node  - node index \in [0, ..., N]
 *         label - possible labeling at this step for the node
 *         O     - dependency vector
 *         P     - permutation matrix
 *
 * Checks if the label respect the constraints or not
 *
 */
bool check_constraints (int node, int label, const std::vector<std::vector<int>>& O, Eigen::MatrixXf& P)
{
    // std::cout << "Checking node = " << node << std::endl;

    // for (std::vector<std::vector<int>>::size_type i; i < O.size(); ++i)
    // {
    //     std::cout << "O[" << i << "] = ";
    //     for (std::vector<int>::size_type j; j < O[i].size(); ++j) std::cout << O[i][j] << " ";
    //     std::cout << std::endl;
    // }

    for (std::vector<int>::size_type j; j < O[node].size(); ++j)
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

/* Function: nodal_numbering
 *
 * Inputs:
 *         A        - original matrix
 *         P        - permutation matrix
 *         rows_deg - vector of nodes' degree
 *
 * Adapted version of the original algorithm in order to consider the constraints imposed.
 * When labeling the adjacent nodes, we check if their labeling respect the constraints
 * or not. If so, we label them accordingly. Otherwise we keep a queue with all the remaining
 * nodes that haven't been labeled already, checking again at each step (i.e., when considering
 * other nodes if labeling is possible)
 *
 */
void nodal_numbering (const Eigen::MatrixXf& A, Eigen::MatrixXf& P, const std::vector<int>& rows_deg,
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

/* Function: matrix_deg
 *
 * Inputs:
 *         R - resulting matrix after applying algorithm
 *
 * Outputs:
 *         max_deg - matrix degree (i.e. maximum{rows_deg})
 */
int matrix_deg (const Eigen::MatrixXf& R)
{
    std::vector<int> rows_deg(R.rows(), 0);

    for (int i = 0; i < R.rows(); ++i)
    {
        for(int j = 0; j < R.cols(); ++j)
        {
            if (R(i,j) != 0 && i != j) rows_deg[i] += 1;
        }
    }

    auto max_deg = *std::max_element(rows_deg.begin(), rows_deg.end());

    return max_deg;
}

/* Function: compute_rows_deg
 *
 * Inputs:
 *         A - original matrix representing the graph
 *
 * Outputs:
 *         rows_deg - for each node (i.e., column or row index of matrix A)
 *                    its degree (number of non-zero off-diagonal elements
 *                    in column or row) is computed
 *
 */
std::vector<int> compute_rows_deg (const Eigen::MatrixXf& A)
{
    std::vector<int> rows_deg(A.rows(), 0); // Assuming A is square

    for(int j = 0; j < A.cols(); ++j)
    {
        for (int i = 0; i < A.rows(); ++i)
        {
            if (A(i,j) != 0 && i != j) rows_deg[i] += 1;
        }
    }

    return rows_deg;
}

/* Function: compute_matrices
 *
 * Inputs: starting nodes - possible starting nodes considering the lowest
 *                          degree strategy (keep in mind it's not always optimal)
 *         A              - original matrix representing the graph
 *         P              - permutation matrix
 *         R              - lowest degree resulting matrix (i.e. P*A*P^T)
 *
 * Perform the nodal numbering algorithm in each one of the lowest degree nodes
 * and check which one generates the lowest degree resulting matrix
 */
void compute_matrices (const Eigen::MatrixXf& A, Eigen::MatrixXf& P, Eigen::MatrixXf& R,
                       const std::vector<std::vector<int>>& O)
{
    std::vector<int> rows_deg = compute_rows_deg(A);

    P(0, 0) = 1;                          // Starting node = 0 is labeled as 0
    nodal_numbering(A, P, rows_deg, O);   // Execute algorithm on A
    R = (P*A*P.transpose());              // Compute resulting matrix R

    std::cout << "dim(A) = " << A.rows() << "x" << A.cols() << std::endl;
    std::cout << "bandwidth(A) = " << matrix_bandwidth(A) << std::endl;
    std::cout << "bandwidth(R) = " << matrix_bandwidth(R) << std::endl;
}
