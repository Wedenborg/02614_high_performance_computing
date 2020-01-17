/* gauss_omp.h - Poisson problem
 *
 */
#ifndef _GAUSS_OPENMP_H
#define _GAUSS_OPENMP_H

int Gauss_openMP(int N, double ***u, double ***v, double ***f, int iter_max, double tolerance);

#endif
