#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void matmult_mnk(int m, int n, int k, double **A,double **B,double **C){
/*
    int m = 3;
    int n = 2;
    int k = 2;

    double A[3][2] = {{2,2},{2,2},{2,2}}; // m x k
    double B[2][2] = {{1,1},{1,1}}; // k x n
    double C[3][2] = {{0}}; // m x n
*/

    
    double C[m][n] = {{0}};
    for(int i=0; i<m;i++){
        for (int j = 0; j<n;j++){
            for (int h= 0;h<k;h++){
                C[i][j] +=  A[i][h]*B[h][j];
            }    
        }
    }
    
/*
    for(int i=0; i<m;i++){
        printf("\n");
        for (int j = 0; j<n;j++){
            printf("%lf ",C[i][j]);
        }
    }
*/
}

