/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi_gpu2.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBIGPU4_H
#define _JACOBIGPU4_H
__global__ 
void jacobi_stopTest(int N, double ***u, double ***v, double ***f, int iter_max, double *res);
#endif
