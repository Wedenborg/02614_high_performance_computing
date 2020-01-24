/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

__device__ 
void jacobi_gpu(int N, double ***u, double ***v, double ***f, int iter_max) {
    //int counter = 0;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    //v[i][j][k] = u[i][j][k];

    //} while (counter <iter_max);
    if(i > 0 && j > 0 && k > 0 && i<N-1 && j<N-1 && k<N-1){
            //v[i][j][k] = u[i][j][k];
            u[i][j][k] = 1./6.*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + 1./((N)*(N)) * f[i][j][k]);    
            
            printf("i=%i j=%i k=%i | u=%f v=%f f=%f\n", i, j, k, u[i][j][k], v[i][j][k], f[i][j][k]);
        }
}

// Kernel to be launched on a single thread
__global__
void jacobi_per_elem(int N, double ***u, double ***v, double ***f, int iter_max)
{
    jacobi_gpu(N, u, v, f, iter_max);
} 