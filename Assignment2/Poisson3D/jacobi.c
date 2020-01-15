/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>

void jacobi(int N, double ***u, double ***v, double ***f, int iter_max) {
    double cond =0.005;
    double stopTest = 100000;

    while( sqrt(stopTest)<iter_max){

    for( int i =0; i < N; i++){
        for( int j = 0; j < N; j++){
            for( int k = 0; k < N; k++){
                v[i][j][k] = u[i][j][k];
            }
        }
    }
        for( int i =1; i < N-1; i++){
            for( int j = 1; j < N-1; j++){
                for( int k = 1; k < N-1; k++){
	                u[i][j][k]= 0.2*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + N*N * f[i][j][k]); //formula and matrix      
                    stopTest +=(u[i][j][k]-v[i][j][k])*(u[i][j][k]-v[i][j][k]);
	                  }
                  }
            }
      }
    free(v);
    }

