#include <Eigen/Dense>
#include <algorithm>
#include <vector>

/*******************************************************************************
 * Calculates bandwidth of a particular matrix.
 *
 *
 * @param R Matrix.
 ******************************************************************************/
int matrix_bandwidth (const Eigen::MatrixXf& R)
{
    std::vector<int> rows_bandwidth(R.rows(), 0);

    for (int i = 0; i < R.rows(); ++i)
    {
        for(int j = 0; j < R.cols(); ++j)
        {
            if (R(i,j) != 0 && i != j && (fabs(j-i) > rows_bandwidth[i]))
            {
                rows_bandwidth[i] = fabs(j-i);
            }
        }
    }

    auto max_bandwidth = *std::max_element(rows_bandwidth.begin(), rows_bandwidth.end());

    return max_bandwidth;
}

void parents_from_matrix(const Eigen::MatrixXf& A, std::vector<std::vector<int>>& parents)
{
    int N = A.rows();

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (A(i, j) != 0 && i != j) parents[i].push_back(j);
        }
    }

    // for (std::vector<std::vector<int>>::size_type i = 0; i < parents.size(); ++i)
    // {
    //     std::cout << i << " -> ";
    //     for (std::vector<int>::size_type j = 0; j < parents[i].size(); ++j)
    //     {
    //         std::cout << parents[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

void children_from_matrix(const Eigen::MatrixXf& A, std::vector<std::vector<int>>& children)
{
    int N = A.rows();

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (A(i, j) != 0 && i != j) children[j].push_back(i);
        }
    }
}

void apply_symmetry (Eigen::MatrixXf& A)
{
    for (int i = 0; i < A.rows(); ++i)
        {
            for (int j = 0; j < A.cols(); ++j)
                {
                    if (A(i,j) != 0) A(j, i) = A(i, j);
                }
        }
}
