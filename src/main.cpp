/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include "bandwidth_minimization.h"
#include "crank_nicolson.h"
#include "compare_solvers.h"

int main (void)
{
    // Cuthill-Mckee algorithm
    // run_tests(20);

    // 1D diffusion
    // int N = 146;
    // float L = 1;
    // float dt = 5.e-4;
    // int nsteps = 620;
    // diffusion_1d_results_to_csv(N, L, dt, nsteps);

    // 2D diffusion
    // int N = 18;
    // float L = 1;
    // float dt = 5.e-4;
    // int nsteps = 50;
    // diffusion_2d_results_to_csv(N, L, dt, nsteps);
    // diffusion_2d_compare(N, L, dt, nsteps);

    // compare solvers LU decomposition and tridiagonal
    compare_lubksb_tridag(600);
}
