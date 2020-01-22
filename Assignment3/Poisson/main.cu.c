/* main.c - Poisson problem in 3D
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _JACOBIGPU1
#include "jacobi_gpu1.h"
#endif


#define N_DEFAULT 100

int
main(int argc, char *argv[]) {

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	*c**h_u = NULL;
    double *** h_v = NULL;
    double *** h_f = NULL;
    int i,j,k;


    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    N = N+2;
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

// speed up test 
    double ts, te;
    //int Threads;
    //Threads = N-2;
    //N = 200 + 2;

    // allocate memory
    if ( (h_u = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }

    // Initialize first guess as zero. 
    if ( (h_v = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array v: allocation failed");
        exit(-1);
    }

    // Allocating a f matrix
    if ( (h_f = d_malloc_3d(N, N, N)) == NULL ) {
    perror("array f: allocation failed");
    exit(-1);
    }

    
    // init u, f
    #pragma omp parallel for default(none) shared(u,f, N, start_T) private( i, j, k)
    //for( i =0; i < N; i++){
        for( j = 0; j < N; j++){
            for( k = 0; k < N; k++){
                if( i==0 || i==N-1 || j == N-1|| k == 0 || k == N-1  ){
                    h_u[i][j][k] = 20;
                }
                else if (j == 0){
                    h_u[i][j][k] = 0;
                }
                else
                {
                    h_u[i][j][k] = start_T ;
                }
                if( i >= 0 && i <= 5./8.*N*0.5 - 1 && j >= 0 && j <= .5*N*0.5-1 && k >= 1./3.*0.5*N-1  && k <= N*0.5-1 ){
                    h_f[i][j][k] = 200;
                }
                else{
                    h_f[i][j][k] = 0;
                }
            }
        }
    //}
    // slut omp

    // for( int i =0; i < N; i++){
    //     for( int j = 0; j < N; j++){
    //         for( int k = 0; k < N; k++){
    //             if( i >= 0 && i <= 5./8.*N*0.5 - 1 && j >= 0 && j <= .5*N*0.5-1 && k >= 1./3.*0.5*N-1  && k <= N*0.5-1 ){
    //                 f[i][j][k] = 200;
    //             }
    //             else{
    //                 f[i][j][k] = 0;
    //             }
    //         }
    //     }
    // }





    #define N 2880

    // CPU reference transpose for checking result
    jacobi(N, h_u, h_v, h_f, iter_max,tolerance);

    ts = omp_get_wtime();
    // Transfer matrix to device
    cudaMemcpy(h_u, d_u, u_size, cudaMemcpyHostToDevice);
    cudaMemcpy(h_v, d_v, v_size, cudaMemcpyHostToDevice);
    cudaMemcpy(h_f, d_f, f_size, cudaMemcpyHostToDevice);
    jacobi_serial<<<1, 1>>>(N, d_u, d_v, d_f, iter_max,tolerance);
    cudaDeviceSynchronize();
    // Transfer result to host

    cudaMemcpy(h_u, d_u, u_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_v, d_v, v_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_f, d_f, f_size, cudaMemcpyDeviceToHost);   
    // Allocating a uu if the Gauss-Seidel


    te = omp_get_wtime() - ts;
    printf("%d ", Threads);
    printf("%lf \n",te);



    // dump  results if wanted 
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    cudaFree(d_u);
    cudaFree(d_v);
    cudaFree(d_f);

    return(0);
}
