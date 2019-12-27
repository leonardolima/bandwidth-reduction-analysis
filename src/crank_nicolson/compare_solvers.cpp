/******************************
 *                            *
 *    Leonardo Lima, 2019     *
 *                            *
/******************************/

#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include "lu_decomposition.h"
#include "tridiagonal.h"
#include "../io.h"

/* Function: compare_lubksb_tridag
 *
 * Inputs: N - maximum matrix dimension
 *
 */
void compare_lubksb_tridag(int N)
{
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N/10, 3);

    int j = 0;

    for (int i = 900; i < N+900; i=i+10)
    {
        Eigen::RowVectorXf row_vec(3);

        row_vec[0] = i;

        auto start = std::chrono::high_resolution_clock::now();

        Eigen::MatrixXf a = Eigen::MatrixXf::Random(i, i);
        Eigen::VectorXf indx = Eigen::VectorXf::Zero(i);
        double c = 0.0;
        ludcmp(a, indx, c);
        Eigen::VectorXf b = Eigen::VectorXf::Zero(i);
        lubksb(a, indx, b);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        auto duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);

        row_vec[1] = duration_ms.count();

        start = std::chrono::high_resolution_clock::now();

        Eigen::VectorXf d = Eigen::VectorXf::Random(i);
        Eigen::VectorXf e = Eigen::VectorXf::Random(i);
        Eigen::VectorXf f = Eigen::VectorXf::Random(i);
        Eigen::VectorXf r = Eigen::VectorXf::Random(i);
        Eigen::VectorXf u = Eigen::VectorXf::Zero(i);
        tridag(d, e, f, r, u);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration<float>(stop - start);
        duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);

        row_vec[2] = duration_ms.count();

        R.row(j) << row_vec[0], row_vec[1], row_vec[2];
        j++;
    }

    matrix_to_csv(R);
}
