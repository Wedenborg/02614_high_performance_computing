/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>



void jacobi(int N, double ***u, int max_iter) {

    double cond =0.005;
    double stopTest = 100000;
    double *** v = NULL;

    if ( (v = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);



    while( sqrt(stopTest)<max_iter){

    for( int i =0; i < N; i++){
        for( int j = 0; j < N; j++){
            for( int k = 0; k < N; k++){
                v[i][j][k] = u[i][j][k];
            }
        }
    }
        for( int i =0; i < N; i++){
            for( int j = 0; j < N; j++){
                for( int k = 0; k < N; k++){
	                u[i][j][k]= 1/6*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1]); //formula and matrix
	                stopTest +=(u[i][j][k]-v[i][j][k])*(u[i][j][k]-v[i][j][k]);

	                  }

                  }
            }
      }
    }

