/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include "gauss_jordan.h"

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

void compute_inverse (Eigen::MatrixXf& M)
{
    M = M.inverse();
    // int k = 0, c = 0, flag = 0;

    // for (int i = 0; i < M.rows(); ++i)
    // {
    //     if (M(i,i) == 0)
    //     {
    //         c = 1;

    //         while (M(i + c, i) == 0 && (i + c) < M.rows()) ++c;
            
    //         if ((i + c) == M.rows())
    //         {
    //             flag = 1;
    //             break;
    //         }

    //         for (int j = i, k = 0; k <= M.rows(); ++k)
    //         {
    //             int temp = M(j, k);
    //             M(j, k) = M(j + c, k);
    //             M(j + c, k) = temp;
    //         }
    //     }

    //     for (int j = 0; j < M.cols(); ++j)
    //     {
    //         if (i != j)
    //         {
    //             float p = M(j, i) / M(i, i);

    //             for (k = 0; k <= M.cols(); k++)
    //             {
    //                 M(j, k) = M(j, k) - (M(i, k)) * p;
    //             }
    //         }
    //     }
    // }
    
    // return flag;
}
