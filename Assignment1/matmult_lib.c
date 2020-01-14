 #include "cblas.h"
 #include <stdio.h>
 /*
 cblas_dgemm(	
    char	TRANSA,
    char	TRANSB,
    int 	M,
    int	    N,
    int 	K,
    double  ALPHA,
    double  A,
    int 	lda,
    double 	B,
    int 	LDB,
    double  BETA,
    double  C,
    int 	ldc
);
*/

void matmult_lib (int m, int n, int k, double **A, double **B, double **C){
  /*
    int m = 2, n=2, k=2;
    double A[2][2] = {{2,2},{2,2}};
    double B[2][2] = {{1,1},{1,1}};
    double C[2][2] = {{0}};
  */
    int TRANS = 111;
    int order = 101;

    double alpha = 1.0;
    int beta = 0;

    int lda = m;
    int ldb = n;
    int ldc = m;
    
    

    cblas_dgemm(order,TRANS,TRANS, m,n, k ,alpha, *A, lda, *B, ldb, beta, *C, ldc);
/*
    for(int i=0; i<m;i++){
        printf("\n");
        for (int j = 0; j<n;j++){
            printf("%lf ",C[i][j]);
        }
    }
*/
}

