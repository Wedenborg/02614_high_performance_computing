/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// en for gpu1 og en for gpu2
__device__ 
void jacobi_gpuDual1(int N, double ***u, double ***v, double ***f, double ***v2, int iter_max) {
    //int counter = 0;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    //v[i][j][k] = u[i][j][k];

    //} while (counter <iter_max);
    // if for halvdelen af gpu (kun halvdelen af z)
    if(i > 0 && j > 0 && k > 0 && i<((N/2)-1) && j<(N-1) && k<(N-1)){
        //v[i][j][k] = u[i][j][k];
        u[i][j][k] = 1./6.*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + 1./((N)*(N)) * f[i][j][k]);    
        }
    if(i > 0 && j > 0 && k > 0 && i==(N/2-1) && j<(N-1) && k<(N-1)){
        u[i][j][k] = 1./6.*(v[i-1][j][k]+v2[0][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1]+ 1./((N)*(N))* f[i][j][k]);    

    }
}
__device__
void jacobi_gpuDual2(int N, double ***u, double ***v, double ***f, double ***v1, int iter_max) {
    //int counter = 0;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    //v[i][j][k] = u[i][j][k];

    //} while (counter <iter_max);
    // if for halvdelen af gpu (kun halvdelen af z)
    // Den skal være større end k/2
    if(i > 0 && j > 0 && k > 0 && i<((N/2)-1) && j<(N-1) && k<(N-1)){
        //v[i][j][k] = u[i][j][k];
        u[i][j][k] = 1./6.*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + 1./((N)*(N)) * f[i][j][k]);    
        
        //printf("i=%i j=%i k=%i | u=%f v=%f f=%f\n", i, j, k, u[i][j][k], v[i][j][k], f[i][j][k]);
    }
    if(j > 0 && k>0 && i==0 &&j<N-1 && k<N-1){
        u[i][j][k] = 1./6.*(v1[(N/2)-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + 1./((N)*(N)) * f[i][j][k]);    
        
        //printf("i=%i j=%i k=%i | u=%f v=%f f=%f\n", i, j, k, u[i][j][k], v[i][j][k], f[i][j][k]);
    }
}
// Kernel to be launched on a single thread
__global__
void jacobi_dual1(int N, double ***u, double ***v, double ***f, double ***v2, int iter_max)
{
    jacobi_gpuDual1(N, u, v, f, v2, iter_max);
} 

__global__
void jacobi_dual2(int N, double ***u, double ***v, double ***f, double ***v1, int iter_max)
{
    jacobi_gpuDual2(N, u, v, f, v1,iter_max);
} 
