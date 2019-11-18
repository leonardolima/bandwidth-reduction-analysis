#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include <fstream>

void diffusion_1d (int, float, float, int, std::ofstream&);

void diffusion_2d (int, float, float, int);

void diffusion_results_to_csv (int, float, float, int);

#endif
