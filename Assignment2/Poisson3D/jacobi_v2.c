/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void jacobi(int N, double ***u, double ***v, double ***f, int iter_max, double tolerance) {
    double stopTest = 100000;
    int counter =0;
    int i,j,k;

    while(stopTest>tolerance && counter < iter_max){
        stopTest =0.0;

        #pragma omp parallel default(none) shared(u, v, f, N, stopTest) private(i, j, k) 
        {
        #pragma omp for
        for( i =0; i < N; i++){
            for( j = 0; j < N; j++){
                for( k = 0; k < N; k++){
                    v[i][j][k] = u[i][j][k];
                }
            }
        }

        #pragma omp for reduction(+: stopTest)
        for( i =1; i < N-1; i++){
            for( j = 1; j < N-1; j++){
                for( k = 1; k < N-1; k++){
	                u[i][j][k] = 1./6.*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + 1./((N-2)*(N-2)) * f[i][j][k]); //formula and matrix      

                    stopTest +=(u[i][j][k]-v[i][j][k])*(u[i][j][k]-v[i][j][k]);
                    
	            }
            }
        }
        } //End Parallel
        counter++;
        
    }
}

