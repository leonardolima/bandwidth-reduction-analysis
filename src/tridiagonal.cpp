/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include <cmath>
#include <iostream>
#include <vector>
#include "tridiagonal.h"

/* Function: compute_tridiagonal
 *
 * Inputs: A - original matrix
 *         R - tridiagonal matrix
 *
 * Numerical Recipes' implementation of Householder's method.
 */

void compute_tridiagonal (Eigen::MatrixXf& A, Eigen::MatrixXf& R)
{
    int n = A.rows(), l = 0;
    double scale, hh, h, g, f;

    std::vector<int> d(n, 0), e(n, 0);

    for (int i = n-1; i > 0; --i)
    {
        l = i-1;
        h = scale = 0.0;
        if (l > 0)
        {
            for (int k = 0; k < l+1; ++k) scale += std::fabs(A(i, k));
            if (scale == 0.0) e[i] = A(i, l);
            else {
                for (int k = 0; k < l+1; ++k)
                {
                    A(i, k) /= scale;
                    h += A(i, k) * A(i, k);                    
                }
                f = A(i, l);
                g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i] = scale*g;
                h -= f*g;
                A(i, l) = f-g;
                f = 0.0;
                for (int j = 0; j < l+1; ++j)
                {
                    A(j, i) = A(i, j)/h;
                    g = 0.0;
                    for (int k = 0; k < j+1; ++k) g += A(j, k)*A(i, k);
                    for (int k = j+1; k < l+1; ++k) g += A(k, j)*A(i, k);
                    e[j] = g/h;
                    f += e[j]*A(i, j);
                }
                hh = f/(h+h);
                for (int j = 0; j < l+1; ++j)
                {
                    f = A(i, j);
                    e[j] = g = e[j]-hh*f;
                    for (int k = 0; k < j+1; ++k) A(j, k) -= (f*e[k]+g*A(i, k));
                }
            }
        } 
        else e[i] = A(i, l);
        d[i] = h;
    }

    d[0] = 0.0;
    e[0] = 0.0;

    for (int i = 0; i < n; ++i)
    {
        l = i;
        if (d[i] != 0.0)
        {
            for (int j = 0; j < l; ++j)
            {
                g = 0.0;
                for (int k = 0; k < l; ++k) g += A(i, k)*A(k, j);
                for (int k = 0; k < l; ++k) A(k, j) -= g*A(k, i);
            }
        }
        d[i] = A(i, i);
        A(i, i) = 1.0;
        for (int j = 0; j < l; j++) A(j, i) = A(i, j) = 0.0;
    }

    std::cout << "d = " << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << d[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "e = " << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << e[i] << " ";
    }
    std::cout << std::endl;

    // Construct R
    for (int i = 0; i < n; ++i)
    {
        R(i, i) = d[i];
        if (i < n-1)
        {
            // We ignore the first element of e
            R(i, i+1) = e[i+1];
            R(i+1, i) = e[i+1];
        }
    }
}
