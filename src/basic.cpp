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
