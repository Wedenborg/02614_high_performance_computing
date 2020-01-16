/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void gauss_seidel(int N, double ***u, double ***v, double ***f, int iter_max, double tolerance) {
    double stopTest = 100000;
    int counter =0 ;

    while(sqrt(stopTest)>tolerance && counter < iter_max){
        stopTest =0.0;
        for( int i =1; i < N-1; i++){
            for( int j = 1; j < N-1; j++){
                for( int k = 1; k < N-1; k++){
	                u[i][j][k] = 1./6.*(u[i-1][j][k]+u[i+1][j][k]+u[i][j-1][k]+u[i][j+1][k]+u[i][j][k-1]+u[i][j][k+1] + 1./((N-2)*(N-2)) * f[i][j][k]); //formula and matrix      
                    //printf("u(%d,%d,%d) = %lf \n",i,j,k, u[i][j][k]);

                    //printf("v(%d,%d,%d) = %lf \n",i,j,k, v[i][j][k]);

                    stopTest +=(u[i][j][k]-v[i][j][k])*(u[i][j][k]-v[i][j][k]);
                    
	            }
            }
            
        }
	//printf("%d ",counter);
        //printf("%lf \n ",stopTest);
         
        counter++;
    }
    printf("Total iterations: %d ",counter);
}


