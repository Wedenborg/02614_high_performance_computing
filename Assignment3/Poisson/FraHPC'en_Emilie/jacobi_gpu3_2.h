/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi_gpu2.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBIGPU3_H
#define _JACOBIGPU3_H

__global__ 
void jacobi_dual2(int N, double ***u, double ***v, double ***f,double ***v1, int iter_max);
#endif
