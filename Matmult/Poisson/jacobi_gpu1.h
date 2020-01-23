/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi_gpu1.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBIGPU1_H
#define _JACOBIGPU1_H

__global__
void jacobi_serial(int N, double ***u, double ***v, double ***f, int iter_max);

#endif
