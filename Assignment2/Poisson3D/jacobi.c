/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
jacobi(int N, double ***u) {

    double ***u_old = NULL;
    double cond =0.005;
    double stopTest = 100000;

    if ( (u_old = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);

    for( int i =0; i < N; i++){
        for( int j = 0; j < N; j++){
            for( int k = 0; k < N; k++){
              	  while( sqrt(stopTest)<cond){
	                  stopTest +=(u_new[i][j][k]-u_old[i][j][k])*(u_new[i][j][k]-u_old[i][j][k])
	                  u_old[i][j][k] = u_new[i][j][k];
	                  u_new = 4 //formula and matrix


	                  }

                  }
            }
      }
    }
    
