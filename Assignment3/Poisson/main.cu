/* main.c - Poisson problem in 3D
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <helper_cuda.h>
#include <stdio.h>
#include <stdlib.h>

#include "alloc3d_gpu.h"
#include "transfer3d_gpu.h"
#include "alloc3d.h"
#include "print.h"

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _JACOBIGPU1
#include "jacobi_gpu1.h"
#endif

#ifdef _JACOBIGPU2
#include "jacobi_gpu2.h"
#endif

//#define N 2880
#define N_DEFAULT 1000

int main(int argc, char *argv[]) {
    printf("%d \n", 1000000);
    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int     output_type =1;
    #ifdef _JACOBIGPU1
    char    *output_prefix ="Jacobi_gpu1.res";
    #endif

    #ifdef _JACOBIGPU2
    char    *output_prefix ="Jacobi_gpu1.res";
    #endif

    char    *output_ext = "";
    char    output_filename[FILENAME_MAX];


    double *** h_u = NULL;
    double *** h_v = NULL;
    double *** h_f = NULL;
    double *** d_u = NULL;
    double *** d_v = NULL;
    double *** d_f = NULL;


    int i,j,k; // indicies

    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    N = N+2;
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc==6){
        output_type = atoi(argv[5]);
    }
    printf("%d \n",N);
    // speed up test 
    double ts, te;
    
    // allocate memory in CPU
    if ( (h_u = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array h_u: allocation failed");
        exit(-1);
    }
    // Initialize first guess as zero. 
    if ( (h_v = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array h_v: allocation failed");
        exit(-1);
    }
    // Allocating a f matrix
    if ( (h_f = d_malloc_3d(N, N, N)) == NULL ) {
    perror("array h_f: allocation failed");
    exit(-1);
    }
    

    // allocate memory in GPU
    if ( (d_u = d_malloc_3d_gpu(N, N, N)) == NULL ) {
            perror("array d_u: allocation failed");
            exit(-1);
    }
    
        // Initialize first guess as zero. 
    if ( (d_v = d_malloc_3d_gpu(N, N, N)) == NULL ) {
            perror("array d_v: allocation failed");
            exit(-1);
    }
    
        // Allocating a f matrix
    if ( (d_f = d_malloc_3d_gpu(N, N, N)) == NULL ) {
        perror("array d_f: allocation failed");
        exit(-1);
    }
    
    // init u, f
    #pragma omp parallel for default(none) shared(h_u,h_f, N, start_T) private( i, j, k)
    for( i =0; i < N; i++){
        for( j = 0; j < N; j++){
            for( k = 0; k < N; k++){
                if( i==0 || i==N-1 || j == N-1|| k == 0 || k == N-1  ){
                    h_u[i][j][k] = 20;
                } else if (j == 0){
                    h_u[i][j][k] = 0;
                } else {
                    h_u[i][j][k] = start_T ;
                }
                if( i >= 0 && i <= 5./8.*N*0.5 - 1 && j >= 0 && j <= .5*N*0.5-1 && k >= 1./3.*0.5*N-1  && k <= N*0.5-1 ){
                    h_f[i][j][k] = 200;
                } else{
                    h_f[i][j][k] = 0;
                }
            }
        }
    }
    // slut omp

    for( int i =0; i < N; i++){
         for( int j = 0; j < N; j++){
             for( int k = 0; k < N; k++){
                 if( i >= 0 && i <= 5./8.*N*0.5 - 1 && j >= 0 && j <= .5*N*0.5-1 && k >= 1./3.*0.5*N-1  && k <= N*0.5-1 ){
                     h_f[i][j][k] = 200;
                 } else{
                     h_f[i][j][k] = 0;
                }
             }
         }
     }
    
    // Transfer matrix to device
    transfer_3d(d_u, h_u, N, N, N , cudaMemcpyHostToDevice);
    transfer_3d(d_v, h_v, N, N, N , cudaMemcpyHostToDevice);
    transfer_3d(d_f, h_f, N, N, N , cudaMemcpyHostToDevice);   

    #ifdef _JACOBIGPU1

    ts = omp_get_wtime();
    jacobi_serial<<<1, 1>>>(N, d_u, d_v, d_f, iter_max);
    te = omp_get_wtime() - ts;
    printf("%lf \n",te);

    #endif


    #ifdef _JACOBIGPU2

    #define K 10
    
    dim3 dimGrid(ceil((double)N/K),ceil((double)N/K),ceil((double)N/K));
    dim3 dimBlock(K,K,K);
    ts = omp_get_wtime();
    int counter =0;
    do{
        
    jacobi_per_elem<<<dimGrid,dimBlock>>>(N-2, d_u, d_v, d_f, iter_max);
    //jacobi_per_elem<<<1,1>>>(N+2, d_u, d_v, d_f, iter_max);
    checkCudaErrors(cudaDeviceSynchronize());
    counter++;

    } while (counter <iter_max);
    te = omp_get_wtime() - ts;
    printf("%lf \n",te);

    #endif 

    // checkCudaErrors(cudaDeviceSynchronize());
    
    // Transfer result to host
    transfer_3d(h_u, d_u, N, N, N , cudaMemcpyDeviceToHost);
    transfer_3d(h_v, d_v, N, N, N , cudaMemcpyDeviceToHost);
    transfer_3d(h_f, d_f, N, N, N , cudaMemcpyDeviceToHost); 

    // print h_U til .vyk burde virke... 



    // de-allocate memory
    //free_gpu(d_u);
    //free_gpu(d_v);
    //free_gpu(d_f);

    cudaFree(d_u);
    cudaFree(d_v);
    cudaFree(d_f);

    //cudaFreeHost(h_u);
    //cudaFreeHost(h_v);
    //cudaFreeHost(h_f);



    switch(output_type){
        case 0:
        // no output 
        break;

        case 1:
        output_ext=".vtk";
        sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N, h_u);
	    break;
    }
    free(h_u);
    free(h_v);
    free(h_f);

    return(0);
}
