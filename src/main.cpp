/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include "bandwidth_minimization.h"
#include "crank_nicolson.h"

int main (void)
{
    //execute_algorithm(20);

    int N = 51;
    float L = 1;
    float dt = 5.e-4;
    int nsteps = 620;
    diffusion_1d(N, L, dt, nsteps);
}
