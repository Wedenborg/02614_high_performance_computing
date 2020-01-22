#include <stdio.h>
extern "C"{

#include "cblas.h"

__global__ void matcal_1(int m,int n,int k, double *A,double *B, double *C){


    for (int i = 0; i<m;i++){
        for (int j= 0;j<n;j++){
            C[i*n + j] =0;
        }    
    }
    for(int i=0; i<m;i++){
        for (int h = 0; h<k;h++){
            for (int j= 0;j<n;j++){
                C[i*n + j] +=  A[i*k + h]*B[h*n + j];
            }    
        }
    }
};

void matmult_lib (int m, int n, int k, double *A, double *B, double *C){



    double alpha = 1.0;
    int beta = 0;

    int lda = k;
    int ldb = n;
    int ldc = n;
    
    

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, m,n, k ,alpha, A, lda, B, ldb, beta, C, ldc);
}

void matmult_gpu1(int m, int n, int k, double *h_A,double *h_B,double *h_C){
    /*
        int m = 3;
        int n = 2;
        int k = 2;
    
        double A[3][2] = {{2,2},{2,2},{2,2}}; // m x k
        double B[2][2] = {{1,1},{1,1}}; // k x n
        double C[3][2] = {{0}}; // m x n
    */


    double *d_A;
    double *d_B;
    double *d_C;


    cudaMalloc((void**)&d_A, m*k * sizeof(double));
    cudaMalloc((void**)&d_B, n*k * sizeof(double));
    cudaMalloc((void**)&d_C, m*n * sizeof(double));

    cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);

    matcal_1<<<1,1>>>(m,n,k,d_A,d_B,d_C);
    cudaDeviceSynchronize();

    cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}


void matmult_gpu2(int m, int n, int k, double *h_A,double *h_B,double *h_C){
/*
    int m = 3;
    int n = 2;
    int k = 2;

    double A[3][2] = {{2,2},{2,2},{2,2}}; // m x k
    double B[2][2] = {{1,1},{1,1}}; // k x n
    double C[3][2] = {{0}}; // m x n
*/


double *d_A;
double *d_B;
double *d_C;


cudaMalloc((void**)&d_A, m*k * sizeof(double));
cudaMalloc((void**)&d_B, n*k * sizeof(double));
cudaMalloc((void**)&d_C, m*n * sizeof(double));

cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);

dim3 dimBlock(32,32,1); // threads per block
dim3 dimGrid(m/32,n/32,1);// blocks in total

matcal_1<<<dimGrid, dimBlock>>>(m,n,k,d_A,d_B,d_C);
cudaDeviceSynchronize();

cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

cudaFree(d_A);
cudaFree(d_B);
cudaFree(d_C);
}

}
