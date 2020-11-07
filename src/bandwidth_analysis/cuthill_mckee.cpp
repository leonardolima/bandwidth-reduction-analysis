#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>
#include "cuthill_mckee.h"
#include "gauss_jordan.h"
#include "tridiagonal.h"
#include "../basic.h"

/*******************************************************************************
 * Checks if a node has already been labeled by the algorithm. This is done by
 * checking if there's a non-zero entry at the node's column on the permutation
 * matrix.
 *
 *
 * @param P Permutation matrix.
 * @param col Node's column index on the permutation matrix.
 ******************************************************************************/
bool node_has_label(const Eigen::MatrixXf& P, unsigned int col)
{
    for (int i = 0; i < P.rows(); ++i) if (P(i,col) != 0) return true;
    return false;
}

int node_label(const Eigen::MatrixXf& P, unsigned int col)
{
    for (int i = 0; i < P.rows(); ++i) if (P(i, col) != 0) return P(i, col);
    return -1;
}

/*******************************************************************************
 * Returns node's index on the adjacency matrix (the matrix representing the graph)
 * given its label.
 *
 *
 * @param P Permutation matrix.
 * @param label Node's label.
 ******************************************************************************/
int node_index_adj_matrix(const Eigen::MatrixXf& P, int label)
{
    // Index on the adjacency matrix corresponds to column's index
    for(int j = 0; j < P.cols(); ++j) if (P(label, j) != 0) return j;
    return -1;
}

/*******************************************************************************
 * TODO
 * There might be something wrong here
 *
 *
 * @param node_row Node's row on adjacency matrix.
 * @param nodes_deg Vector with the degree of each node's rows.
 * @param P Permutation matrix.
 ******************************************************************************/
std::vector<std::pair<int, int>> sort_adj_nodes(const Eigen::VectorXf& adj_nodes,
                                                const std::vector<int>& nodes_deg,
                                                const Eigen::MatrixXf& P)
{
    std::vector<std::pair<int, int>> sorted_adj_nodes;

    for(int i = 0; i < adj_nodes.size(); ++i)
    {
        // We are not interested if the element is on the diagonal
        // or if the node already has a label
        if (adj_nodes[i] != 0 && !node_has_label(P, i))
        {
            std::pair<int, int> pair = std::make_pair(nodes_deg[i], i);
            sorted_adj_nodes.push_back(pair);
        }
    }

    std::sort(sorted_adj_nodes.begin(), sorted_adj_nodes.end());

    return sorted_adj_nodes;
}

/*******************************************************************************
 * Nodal numbering of the Cuthill-McKee algorithm. First, unconnected nodes are
 * labeled. Next, for every node, (1) list the adjacent nodes (2) sort them in
 * ascending order of degree, then (3) label them uniquely and accordingly.
 *
 *
 * @param A Adjacency matrix.
 * @param P Permutation matrix.
 * @param nodes_deg Vector with the degree of each node's rows.
 ******************************************************************************/
void nodal_numbering(const Eigen::MatrixXf& A, Eigen::MatrixXf& P, const std::vector<int>& nodes_deg)
{
    unsigned int new_label = 1;

    // 3. Label unconnected nodes
    for(int j = 0; j < A.cols(); ++j)
    {
        if (nodes_deg[j] == 0) P(new_label++, j) = 1; // Update matrix P accordingly
    }

    // 4. For every node we (i) list the adjacent nodes (ii) sort them in ascending order
    // of degree (iii) label them uniquely and accordingly (2, 3, ...)
    for(int j = 0; j < A.cols(); ++j)
    {
        if (nodes_deg[j] > 0)
        {
            int node_index = node_index_adj_matrix(P, j);
            Eigen::VectorXf adj_nodes = A.row(node_index);
            std::vector<std::pair<int, int>> sorted_adj_nodes = sort_adj_nodes(adj_nodes, nodes_deg, P);

            for (std::vector<std::pair<int, int>>::size_type i = 0; i < sorted_adj_nodes.size(); ++i)
            {
                P(new_label++, sorted_adj_nodes[i].second) = 1; // Update matrix P accordingly
            }
        }
    }
}

/*******************************************************************************
 * Computes the matrix degree, i.e., the maximum of each node's degree.
 *
 *
 * @param R Resulting matrix after applying Cuthill-McKee algorithm.
 ******************************************************************************/
int matrix_deg(const Eigen::MatrixXf& R)
{
    std::vector<int> nodes_deg(R.rows(), 0);

    for (int i = 0; i < R.rows(); ++i)
    {
        for(int j = 0; j < R.cols(); ++j)
        {
            if (R(i,j) != 0 && i != j) nodes_deg[i] += 1;
        }
    }

    auto max_deg = *std::max_element(nodes_deg.begin(), nodes_deg.end());
    return max_deg;
}

/*******************************************************************************
 * Select starting nodes, uses the heuristic described on the paper of picking
 * all the lowest degree nodes.
 *
 *
 * @param A Adjacency matrix.
 ******************************************************************************/
std::vector<int> select_starting_nodes(const Eigen::MatrixXf& A)
{
    std::vector<int> starting_nodes;

    // 1. Compute degree of every row, where deg(i) = \Sum_{j!=i} a_{ij} != 0
    std::vector<int> nodes_deg = compute_nodes_deg(A);

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
    for(std::vector<int>::size_type i = 0; i < nodes_deg.size(); ++i)
    {
        if (nodes_deg[i] > 0 && nodes_deg[i] < min_deg) min_deg = nodes_deg[i];
    }

    for(std::vector<int>::size_type i = 0; i < nodes_deg.size(); ++i)
    {
        if (nodes_deg[i] == min_deg) starting_nodes.push_back(i);
    }

    return starting_nodes;
}

/*******************************************************************************
 * Performs the nodal numbering process starting with each one of the lowest
 * degree nodes and checks which one generates the lowest degree resulting
 * matrix.
 *
 *
 * @param starting_nodes Vector of possible starting nodes, considering the
 *                       lowest degree strategy.
 * @param A Adjacency matrix.
 * @param P Permutation matrix.
 * @param R Lowest degree resulting matrix (R = P*A*P^T).
 ******************************************************************************/
void compute_matrices(const std::vector<int>& starting_nodes, const Eigen::MatrixXf& A,
                      Eigen::MatrixXf& P, Eigen::MatrixXf& R)
{
    int max_bandwidth = std::numeric_limits<int>::max(); // max_bandwidth = \infty

    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(A.rows(), A.cols());
    std::vector<int> nodes_deg = compute_nodes_deg(A);

    for(std::vector<int>::size_type i = 0; i < starting_nodes.size(); ++i)
    {
        P(0, starting_nodes[i]) = 1;           // Starting node is labeled as 0
        nodal_numbering(A, P, nodes_deg);      // Execute algorithm on A
        M = (P*A*P.transpose());               // Compute resulting matrix M
        int M_bandwidth = matrix_bandwidth(M); // Compute resulting matrix bandwidth
        if (M_bandwidth < max_bandwidth)
        {
            R = M;                             // Keep matrix M with lowest bandwidth
            max_bandwidth = M_bandwidth;       // Update max_bandwidth
        }
        P.setZero();                           // Clear P matrix
    }

    std::cout << "dim(A) = " << A.rows() << "x" << A.cols() << std::endl;
    std::cout << "bandwidth(A) = " << matrix_bandwidth(A) << std::endl;
    std::cout << "bandwidth(R) = " << matrix_bandwidth(R) << std::endl;
}

/*******************************************************************************
 * Generate a matrix with 1's on the main diagonal and 1's and 0's randomly
 * distributed outside the main diagonal.
 *
 *
 * @param A Zeroed adjacency matrix.
 ******************************************************************************/
void generate_binary_random_matrix(Eigen::MatrixXf& A)
{
    int dimension = A.rows();

    A = (A + Eigen::MatrixXf::Constant(dimension, dimension, 1.))*(1./2.);
    A = (A + Eigen::MatrixXf::Constant(dimension, dimension, 0.));
    A = A.array().round();

    // Fill main diagonal with 1's
    for(int i = 0; i < A.rows(); ++i)
    {
        for (int j = 0; j < A.cols(); ++j) if (i == j) A(i,j) = 1;
    }
}

/*******************************************************************************
 * Performs the Cuthill-McKee on a particular symmetric matrix.
 *
 *
 * @param A Zeroed adjacency matrix.
 * @param P Permutation matrix.
 * @param R Lowest degree resulting matrix (R = P*A*P^T).
 ******************************************************************************/
void run_algorithm(const Eigen::MatrixXf& A, Eigen::MatrixXf& P, Eigen::MatrixXf& R)
{
    std::vector<int> starting_nodes = select_starting_nodes(A);
    compute_matrices(starting_nodes, A, P, R);
}

void run_tests(int dimension)
{
    // Eigen::MatrixXf A(10, 10);
    // A << 1, 1, 0, 1, 0, 0, 0, 0, 1, 0,
    //      1, 1, 1, 0, 0, 0, 0, 0, 1, 0,
    //      0, 1, 1, 0, 1, 0, 0, 0, 1, 0,
    //      1, 0, 0, 1, 1, 1, 0, 0, 1, 1,
    //      0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    //      0, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    //      0, 0, 0, 0, 0, 1, 1, 1, 0, 1,
    //      0, 0, 0, 0, 1, 0, 1, 1, 0, 1,
    //      1, 1, 1, 1, 1, 0, 0, 0, 1, 1,
    //      0, 0, 0, 1, 1, 1, 1, 1, 1, 1;

    // Eigen::MatrixXf A(5, 5);
    // A << 1, 1, 0, 1, 0,
    //      1, 1, 1, 0, 1,
    //      0, 1, 1, 0, 0,
    //      1, 0, 0, 1, 1,
    //      0, 1, 0, 1, 1;

    // Eigen::MatrixXf A = Eigen::MatrixXf::Random(dimension, dimension);
    // generate_binary_random_matrix(A);

    // std::vector<int> starting_nodes = select_starting_nodes(A);

    // Eigen::MatrixXf P = Eigen::MatrixXf::Zero(A.rows(), A.cols());
    // Eigen::MatrixXf R = Eigen::MatrixXf::Zero(A.rows(), A.cols());

    // auto start = std::chrono::high_resolution_clock::now();

    // compute_matrices(starting_nodes, A, P, R);
    // std::cout << "R = " << std::endl;
    // std::cout << R << std::endl;

    // R = Eigen::MatrixXf::Zero(A.rows(), A.cols());

    // auto stop = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration<float>(stop - start);
    // auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    // std::cout << "Execution duration = " << duration_s.count() << "s" << std::endl;
    // std::cout << std::endl;
}
