/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
jacobi(int N, double ***u, int max_iter) {

    double ***u_new = NULL;
    double cond =0.005;
    double stopTest = 100000;
    double *** v = v;

    if ( (uu = d_malloc_3d(N, N, N)) == NULL ) {
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
	                  stopTest +=(u_new[i][j][k]-u_old[i][j][k])*(u_new[i][j][k]-u_old[i][j][k]);
	                   u[i][j][k]= 1/6*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1]); //formula and matrix


	                  }

                  }
            }
      }
    }

