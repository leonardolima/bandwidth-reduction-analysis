/******************************
 *                            *
 *    Leonardo Lima, 2019     *
 *                            *
/******************************/

#include <cmath>
#include <iostream>
#include <vector>
#include "gauss_jordan.h"

/* Function: check_if_inversible
 *
 * Inputs: M - original matrix
 *
 * Check if M*M^T is a diagonal matrix of 1's
 */
bool check_if_inversible (const Eigen::MatrixXf& M)
{
    Eigen::MatrixXf R = M*M.transpose();

    for(int j = 0; j < R.cols(); ++j)
    {
        for (int i = 0; i < R.rows(); ++i)
        {
            if (i != j && R(i,j) != 0) return false;
            if (i == j && R(i,j) != 1) return false;
        }
    }

    return true;
}

/* Function: gauss_jordan
 *
 * Inputs: M - original matrix
 *
 * Implementation of Gauss-Jordan's method.
 * The inverse of M is gradually built up in M itself.
 */
void compute_inverse (Eigen::MatrixXf& M)
{
    int icol = 0, irow = 0;
    double big = 0.0, dum = 0.0, pivinv = 0.0;

    int n = M.rows();
    int m = M.cols();

    Eigen::MatrixXf B = Eigen::MatrixXf::Zero(n, m);

    std::vector<int> indxc(n, 0), indxr(n, 0), ipiv(n, 0);

    for (int i = 0; i < n; ++i)
    {
        big = 0.0;
        for (int j = 0; j < n; ++j)
        {
            if (ipiv[j] != 1)
            {
                for (int k = 0; k < n; ++k)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(M(j, k)) >= big)
                        {
                            big = std::fabs(M(j, k));
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (int l = 0; l < n; ++l) std::swap(M(irow, l), M(icol, l));
            for (int l = 0; l < m; ++l) std::swap(B(irow, l), B(icol, l));
        }

        indxr[i] = irow;
        indxc[i] = icol;

        if (M(icol, icol) == 0.0) std::cout << "GaussJ: Singular Matrix" << std::endl;

        pivinv = 1.0/M(icol, icol);
        M(icol, icol) = 1.0;

        for (int l = 0; l < n; ++l) M(icol, l) *= pivinv;
        for (int l = 0; l < m; ++l) B(icol, l) *= pivinv;

        for (int ll = 0; ll < n; ++ll)
        {
            if (ll != icol)
            {
                dum = M(ll, icol);
                M(ll, icol) = 0.0;
                for (int l = 0; l < n; ++l) M(ll, l) -= M(icol, l)*dum;
                for (int l = 0; l < m; ++l) B(ll, l) -= B(icol, l)*dum;
            }
        }
    }

    for (int l = n-1; l >= 0; --l)
    {
        if (indxr[l] != indxc[l])
        {
            for (int k = 0; k < n; ++k) std::swap(M(k, indxr[l]), M(k, indxc[l]));
        }
    }
}

void gaussj_elim (const Eigen::MatrixXf& H, const Eigen::VectorXf& U, Eigen::VectorXf& b)
{
    Eigen::MatrixXf M = H;
    M.conservativeResize(Eigen::NoChange, M.cols()+1);
    M.col(H.cols()-1) = U;

    int c, N = M.rows();
    float p = 0;

    for (int i = 0; i < N; ++i)
    {
        if (M(i,i) == 0)
        {
            c = 1;
            while ((i + c) < N && M(i+c,i) == 0) c++;
            if ((i + c) == N) break;

            for (int j = i, k = 0; k <= N; ++k) std::swap(M(j, k), M(j+c, k));
        }

        for (int j = 0; j < N; ++j)
        {
            if (i != j)
            {
                p = M(j, i) / M(i, i);

                for (int k = 0; k <= N; ++k) M(j, k) = M(j, k) - (M(i, k)) * p;
            }
        }
    }

    for (int i = 0; i < N; ++i) b[i] = M(i, N) / M(i, i);
}
