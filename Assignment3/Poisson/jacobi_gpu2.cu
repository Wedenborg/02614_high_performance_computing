/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

__device__ 
void jacobi_gpu(int N, double ***u, double ***v, double ***f, int iter_max) {
    int counter =0;
    int i = blockIdx.x * blockDim.x + threadIdx.x+1;
    int j = blockIdx.y * blockDim.y + threadIdx.y+1;
    int k = blockIdx.z * blockDim.z + threadIdx.z+1;
    if(i<N && j<N && k<N){
        //do{
            v[i][j][k] = u[i][j][k];
            u[i][j][k] = 1./6.*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + 1./((N-2)*(N-2)) * f[i][j][k]);    
            
            //counter++;

        //} while (counter <iter_max);
}
}

// Kernel to be launched on a single thread
__global__
void jacobi_per_elem(int N, double ***u, double ***v, double ***f, int iter_max)
{
    jacobi_gpu(N, u, v, f, iter_max);
} 