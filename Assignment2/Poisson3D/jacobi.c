/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void jacobi(int N, double ***u, double ***v, double ***f, int iter_max, double tolerance) {
    double stopTest = 100000;
    int counter =0 ;

    while( sqrt(stopTest)>tolerance && counter < iter_max){
        stopTest =0.0;
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
	                u[i][j][k] = 1./6.*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + (N-2)*(N-2) * f[i][j][k]); //formula and matrix      
                    //printf("u(%d,%d,%d) = %lf \n",i,j,k, u[i][j][k]);

                    //printf("v(%d,%d,%d) = %lf \n",i,j,k, v[i][j][k]);

                    stopTest +=sqrt((u[i][j][k]-v[i][j][k])*(u[i][j][k]-v[i][j][k]));
                    
	            }
            }
            
        }
        //printf("stopTest: %lf \n",stopTest);
        //printf("count: %d \n",counter); 
        counter++;
    }
    free(v);
}

