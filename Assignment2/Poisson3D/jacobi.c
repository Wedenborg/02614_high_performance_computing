/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>

int f(int i,int j,int k){
    if( i >= -1 && i <= 0.375 ){
        if( j >= -1 && j <= -0.5){
            if( k >= -0.66666666667  && k <= 0){
                return 200;
            }
        }
    }
    else{
        return 0;
    }
}

void jacobi(int N, double ***u, int iter_max) {

    double cond =0.005;
    double stopTest = 100000;
    double *** v = NULL;


    // Initialize first guess as zero. 
    if ( (v = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }

    // init u
    for( int i =0; i < N; i++){
        for( int j = 0; j < N; j++){
            for( int k = 0; k < N; k++){
                if( i==1 || i==-1 || j == 1|| k == -1 || k == 1  ){
                    u[i][j][k] = 20;
                }
                else
                {
                    u[i][j][k] = 0;
                }
            }
        }
    }



    while( sqrt(stopTest)<iter_max){

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
	                u[i][j][k]= 1/6*(v[i-1][j][k]+v[i+1][j][k]+v[i][j-1][k]+v[i][j+1][k]+v[i][j][k-1]+v[i][j][k+1] + N*N * f(i,j,k)); //formula and matrix      
                    stopTest +=(u[i][j][k]-v[i][j][k])*(u[i][j][k]-v[i][j][k]);
	                  }
                  }
            }
      }
    }

