#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "lu_decomposition.h"
#include "../basic.h"

/*******************************************************************************
 * Numerical Recipes' implementation of a method for performing LU
 * decomposition.
 *
 *
 * @param a    Matrix.
 * @param indx Output vector that records the row permutation effected
 *             by the partial pivoting
 * @param d    Output is +/-1 depending on whether the number of row
 *             interchanges was even or odd.
 ******************************************************************************/
void ludcmp(Eigen::MatrixXf& a, Eigen::VectorXf& indx, double& d)
{
    const double TINY = 1.0e-20;
    int i, imax, j, k;
    double big, dum, sum, temp;

    int n = a.rows();
    Eigen::VectorXf vv(n);

    d = 1.0;

    for (i = 0; i < n; i++)
    {
        big = 0.0;
        for (j = 0; j < n; j++) if ((temp=fabs(a(i, j))) > big) big = temp;
        if (big == 0.0) std::cout << "LUdcmp is traversing a singular matrix" << std::endl;
        vv[i] = 1.0/big;
    }

    for (j = 0; j < n; j++) {
        for (i = 0;i < j; i++) {
            sum = a(i, j);
            for (k = 0; k < i; k++) sum -= a(i, k)*a(k, j);
            a(i, j) = sum;
        }
        big = 0.0;
        for (i = j; i < n; i++) {
            sum = a(i, j);
            for (k = 0; k < j; k++) sum -= a(i, k)*a(k, j);
            a(i, j) = sum;
            if ((dum = vv[i]*fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 0; k < n; k++) {
                dum = a(imax, k);
                a(imax, k) = a(j, k);
                a(j, k) = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a(j, j) == 0.0) a(j, j) = TINY;
        if (j != n-1) {
            dum = 1.0/(a(j, j));
            for (i = j+1; i < n; i++) a(i, j) *= dum;
        }
    }
}

/*******************************************************************************
 * Numerical Recipes' implementation for solving a linear system of equations
 * of the form A*X = B.
 *
 *
 * @param a    Matrix after LU decomposition is applied.
 * @param indx Input vector that comes from LU decomposition.
 ******************************************************************************/
void lubksb(const Eigen::MatrixXf& a, const Eigen::VectorXf& indx, Eigen::VectorXf& b)
{
    int i, ii = 0, ip, j;
    double sum;

    int n = a.rows();

    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii != 0)
            for (j = ii-1; j < i; j++) sum -= a(i, j)*b[j];
        else if (sum != 0.0)
            ii = i+1;
        b[i]=sum;
    }

    for (i = n-1; i >= 0; i--) {
        sum = b[i];
        for (j = i+1; j < n; j++) sum -= a(i, j)*b[j];
        b[i] = sum/a(i, i);
    }
}

/*******************************************************************************
 * Numerical Recipes' implementation of a band decomposer.
 *
 *
 * @param A Band-diagonal matrix.
 * @param Au Upper triangular matrix.
 * @param Al Lower triangular matrix.
 * @param indx Input vector related to LU decomposition.
 ******************************************************************************/
void band_decomposer (const Eigen::MatrixXf& A, Eigen::MatrixXf& Au, Eigen::MatrixXf& Al,
                      Eigen::VectorXf& indx)
{
    const double TINY = 1.0e-40;
    int m1 = matrix_bandwidth(A), n = A.rows();
    int m2 = m1, mm = m1+m2+1, i, j, k, l;
    double dum, d;

    l = m1;

    for (i = 0; i < m1; ++i)
    {
        for (int j = m1-i; j < mm; ++j) Au(i, j-1) = Au(i, j);
        l--;
        for (int j = mm-l-1; j < mm; j++) Au(i, j) = 0.0;
    }

    d = 1.0;
    l = m1;

    for (k = 0; k < n; ++k)
    {
        dum = Au(k, 0);
        i = k;
        if (l < n) l++;
        for (j = k+1; j < l; ++j)
        {
            if (std::abs(Au(j, 0)) > std::abs(dum))
            {
                dum = Au(j, 0);
                i = j;
            }
        }
        indx[k] = i+1;
        if (dum == 0.0) Au(k, 0) = TINY;
        if (i != k)
        {
            d = -d;
            for (j = 0; j < mm; ++j) std::swap(Au(k, j), Au(i, j));
        }
        for (i = k+1; i < l; ++i)
        {
            dum = Au(i, 0)/Au(k, 0);
            Al(k, i-k-1) = dum;
            for (j = 1; j < mm; ++j)
            {
                Au(i, j-1) = Au(i, j)-dum*Au(k, j);
                Au(i, mm-1) = 0.0;
            }
        }
    }
}

/*******************************************************************************
 * LU decomposition method to solve a system of linear equations of the form
 * A*x = b, considering that A is band-diagonal.
 *
 *
 * @param A Band-diagonal matrix.
 * @param x Vector.
 * @param b Vector.
 ******************************************************************************/
void band_solver (const Eigen::MatrixXf& A, Eigen::VectorXf& x, Eigen::VectorXf& b)
{
    int m1 = matrix_bandwidth(A), n = A.rows();
    int m2 = m1, mm = m1+m2+1, i, j, k, l;
    double dum;

    Eigen::MatrixXf Au = Eigen::MatrixXf::Zero(n, n); // Upper triangular matrix
    Eigen::MatrixXf Al = Eigen::MatrixXf::Zero(n, m1); // Lower triangular matrix
    Eigen::VectorXf indx = Eigen::VectorXf::Zero(n);

    band_decomposer(A, Au, Al, indx);

    l = m1;

    for (k = 0; k < n; ++k) x[k] = b[k];
    for (k = 0; k < n; ++k)
    {
        j = indx[k]-1;
        if (j != k) std::swap(x[k], x[j]);
        if (l < n) l++;
        for (j = k+1; j < l; ++j) x[j] -= Al(k, j-k-1)*x[k];
    }

    l = 1;

    for (i = n-1; i >= 0; --i)
    {
        dum = x[i];
        for (k = 1; k < l; ++k) dum -= Au(i, k)*x[k+i];
        x[i] = dum/Au(i, 0);
        if (l < mm) l++;
    }
}
