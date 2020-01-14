/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
jacobi(int N, double ***u) {

    double ***u_old = NULL;
    double ***u_new = NULL;
    double cond =0.005;
    double stopTest = 100000;
    double *** v = 0;

    if ( (u_old = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);

    for( int i =0; i < N; i++){
        for( int j = 0; j < N; j++){
            for( int k = 0; k < N; k++){
              	  while( sqrt(stopTest)<cond){
	                  stopTest +=(u_new[i][j][k]-u_old[i][j][k])*(u_new[i][j][k]-u_old[i][j][k]);
	                  u_old[i][j][k] = u_new[i][j][k];
	                  v = 1/6*(v[i-1][j][k]+u[i+1][j][k]+u[i][j-1][k]+u[i][j+1][k]+u[i][j][k-1]+u[i][j][k+1]); //formula and matrix


	                  }

                  }
            }
      }
    }
    
