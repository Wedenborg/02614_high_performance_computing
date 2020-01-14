#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

//void matmult_mnk(int m, int n, int k, double **A,double **B,double **C, int bs){
int main(void){

    int m = 4;
    int n = 4;
    int k = 4;
    int bs = 2;

    double A[4][4] = {{2,2,2,2},{2,2,2,2},{2,2,2,2},{2,2,2,2}}; // m x k
    double B[4][4] = {{2,2,2,2},{2,2,2,2},{2,2,2,2},{2,2,2,2}}; // k x n
    double C[4][4] = {{0}}; // m x n


    
    for (int i = 0; i<m;i++){
        for (int j= 0;j<n;j++){
            C[i][j] =0;
        }    
    }
    for(int i=0; i<m;i+=bs){
        for (int j = 0; j<n;j+=bs){
            for (int h= 0;h<k;h+=bs){
                for (int ii=0;ii<MIN(m-i, bs);ii++){
                    for (int jj=0;jj<MIN(n-j, bs);jj++){
                        for (int hh=0;hh<MIN(n-h, bs);hh++){
                            C[i+ii][j+jj] +=  A[i+ii][h+hh]*B[h+hh][j+jj];
                        }

                    }
                }
            }    
        }
    }
    

    for(int i=0; i<m;i++){
        printf("\n");
        for (int j = 0; j<n;j++){
            printf("%lf ",C[i][j]);
        }
    }

}

