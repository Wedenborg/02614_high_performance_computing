/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi_gpu2.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBIGPU2_H
#define _JACOBIGPU2_H
__global__ 
void jacobi_per_elem(int N, double ***u, double ***v, double ***f, int iter_max);
#endif
