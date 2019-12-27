/******************************
 *                            *
 *    Leonardo Lima, 2019     *
 *                            *
/******************************/

#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>

/* Function: node_has_label
 *
 * Inputs:
 *         P   - permutation matrix
 *         col - column's index
 *
 */
bool node_has_label (const Eigen::MatrixXf& P, unsigned int col)
{
    for (int i = 0; i < P.rows(); ++i)
    {
        if (P(i,col) != 0) return true;
    }

    return false;
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
        if (node_row[i] != 0 && !node_has_label(P, i))
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

/* Function: nodal_numbering
 *
 * Inputs:
 *         A        - original matrix
 *         P        - permutation matrix
 *         rows_deg - vector of nodes' degree
 *
 */
void nodal_numbering (const Eigen::MatrixXf& A, Eigen::MatrixXf& P, const std::vector<int>& rows_deg,
                      std::vector<std::vector<int>>& O)
{
    unsigned int new_label = 1;

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
        if (rows_deg[j] > 0)
        {
            int node_row_index = node_index(P, j);
            Eigen::VectorXf node_row = A.row(node_row_index);
            std::vector<std::pair<int, int>> sorted_rows_deg = sort_rows_deg(node_row, rows_deg, P);

            for (std::vector<std::pair<int, int>>::size_type j = 0; j < sorted_rows_deg.size(); ++j)
            {
                P(new_label++, sorted_rows_deg[j].second) = 1; // Update matrix P accordingly
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

/* Function: matrix_bandwidth
 *
 * Inputs:
 *         R - resulting matrix after applying algorithm
 *
 * Outputs:
 *         max_bandwidth - matrix bandwidth
 */
int matrix_bandwidth (const Eigen::MatrixXf& R)
{
    std::vector<int> rows_bandwidth(R.rows(), 0);

    for (int i = 0; i < R.rows(); ++i)
    {
        for(int j = 0; j < R.cols(); ++j)
        {
            if (R(i,j) != 0 && i != j && (j-i) > rows_bandwidth[i]) rows_bandwidth[i] = j-i;
        }
    }

    auto max_bandwidth = *std::max_element(rows_bandwidth.begin(), rows_bandwidth.end());

    return max_bandwidth;
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

/* Function: select_starting_nodes
 *
 * Inputs:
 *         A - original matrix representing the graph
 *
 * Outputs:
 *         starting_nodes - indexes of starting nodes where we
 *                          consider the lowest degree ones
 *                          (this approach is not always optimal)
 *
 */
std::vector<int> select_starting_nodes (const Eigen::MatrixXf& A)
{
    std::vector<int> starting_nodes;

    // 1. Compute degree of every row, where deg(i) = \Sum_{j!=i} a_{ij} != 0
    std::vector<int> rows_deg = compute_rows_deg(A);

    // 2. Select starting node (usually the one with minimum degree)
    // Permuting rows and columns of matrix A is done using the
    // permutation matrix P. P has one non-zero unit element in each
    // row and column. For row i: p_{ij} i is the new node label and
    // j is the original one. We should always keep in mind that
    // optimal starting nodes doesn't necessarely are the lowest
    // degree ones
    auto min_deg = std::numeric_limits<int>::max();

    // We should not consider rows that only have a single element on
    // the main diagonal when computing starting_nodes
    for(std::vector<int>::size_type i = 0; i < rows_deg.size(); ++i)
    {
        if (rows_deg[i] > 0 && rows_deg[i] < min_deg) min_deg = rows_deg[i];
    }

    for(std::vector<int>::size_type i = 0; i < rows_deg.size(); ++i)
    {
        if (rows_deg[i] == min_deg) starting_nodes.push_back(i);
    }

    return starting_nodes;
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
void compute_matrices (const std::vector<int>& starting_nodes, const Eigen::MatrixXf& A,
                       Eigen::MatrixXf& P, Eigen::MatrixXf& R, std::vector<std::vector<int>>& O)
{
    int max_mat_deg = std::numeric_limits<int>::max(); // max_mat_deg = \infty

    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(A.rows(), A.cols());

    std::vector<int> rows_deg = compute_rows_deg(A);

    for(std::vector<int>::size_type i = 0; i < starting_nodes.size(); ++i)
    {
        P(0, starting_nodes[i]) = 1;          // Starting node is labeled as 1
        nodal_numbering(A, P, rows_deg, O);   // Execute algorithm on A
        M = (P*A*P.transpose());              // Compute resulting matrix M
        int mat_deg = matrix_deg(M);          // Compute resulting matrix degree
        if (mat_deg < max_mat_deg) R = M;     // Keep matrix M with lowest degree
        P.setZero();                          // Clear P matrix
    }

    // std::cout << "dim(A) = " << A.rows() << "x" << A.cols() << std::endl;
    // std::cout << "bandwidth(A) = " << matrix_bandwidth(A) << std::endl;
    // std::cout << "bandwidth(R) = " << matrix_bandwidth(R) << std::endl;
}

void run_algorithm (const Eigen::MatrixXf& A, Eigen::MatrixXf& P, Eigen::MatrixXf& R, std::vector<std::vector<int>>& O)
{
    std::vector<int> starting_nodes = select_starting_nodes(A);
    compute_matrices(starting_nodes, A, P, R, O);
}
