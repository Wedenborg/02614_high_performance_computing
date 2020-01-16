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

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#define N_DEFAULT 100

int
main(int argc, char *argv[]) {
    double ts, te;


    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;
    double *** v = NULL;
    double *** f = NULL;


    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    N = N+2;
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }
    int Threads;
    Threads = N-2;

    for(int n =0; n<4; n++){
        N = n*25 + 75+2;
        // allocate memory
        if ( (u = d_malloc_3d(N, N, N)) == NULL ) {
            perror("array u: allocation failed");
            exit(-1);
        }
        // init u
        for( int i =0; i < N; i++){
            for( int j = 0; j < N; j++){
                for( int k = 0; k < N; k++){
                    if( i==0 || i==N-1 || j == N-1|| k == 0 || k == N-1  ){
                        u[i][j][k] = 20;
                    }
                    else if (j == 0){
                        u[i][j][k] = 0;
                    }
                    else
                    {
                        u[i][j][k] = start_T ;
                    }
                }
            }
        }

        // Initialize first guess as zero. 
        if ( (v = d_malloc_3d(N, N, N)) == NULL ) {
            perror("array v: allocation failed");
            exit(-1);
        }

        // Allocating a f matrix
        if ( (f = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array f: allocation failed");
        exit(-1);
        }

        for( int i =0; i < N; i++){
            for( int j = 0; j < N; j++){
                for( int k = 0; k < N; k++){
                    if( i >= 0 && i <= 5./8.*N*0.5 - 1 && j >= 0 && j <= .5*N*0.5-1 && k >= 1./3.*0.5*N-1  && k <= N*0.5-1 ){
                        f[i][j][k] = 200;
                    }
                    else{
                        f[i][j][k] = 0;
                    }
                }
            }
        }





        // Allocating a uu if the Gauss-Seidel
        ts = omp_get_wtime();
        #ifdef _JACOBI
        jacobi(N, u, v, f, iter_max,tolerance);
        #endif
        te = omp_get_wtime() - ts;
        printf("%d ", Threads);
        printf("%d ",N-2);
        printf("%lf \n",te);

        #ifdef _GAUSS_SEIDEL
        gauss_seidel(N, u, v, f, iter_max,tolerance);
        #endif

    }   

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
    free(u);

    return(0);
}
