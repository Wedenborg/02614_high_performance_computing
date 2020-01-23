#include <stdio.h>
#include <stdlib.h>
extern "C"{
#ifdef _OPENMP
#include <omp.h>
#endif
#include "cblas.h"

#define size 8

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

__global__ void matcal_2(int m,int n,int k, double *A,double *B, double *C){

    // 2D thread indices defining row and col of element
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if ( i<m && j<n ){

        double c1 = 0;
        for (int h= 0;h<k;h++){
            c1 +=  A[i*k + h]*B[h*n + j];
        } 
        C[i*n + j] = c1;
    }
};

// below version
// __global__ void matcal_3(int m,int n,int k, double *A,double *B, double *C){

//     // 2D thread indices defining row and col of element
//     int i = 2*(blockIdx.x * blockDim.x + threadIdx.x);
//     int j = blockIdx.y * blockDim.y + threadIdx.y;

//     double c1 = 0;
//     double c2 = 0;

//         if ( i<m && j<n  && (i+1)<m ){
//             for (int h= 0;h<k;h++){
//                 c1 +=  A[i*k + h]*B[h*n + j];

//                 c2 +=  A[(i+1)*k + h]*B[h*n + j];
//             } 
//             C[i*n + j] =c1;
//             C[(i+1)*n + j]=c2;
//         }
//         else if( i<m && j<n){
//             for (int h= 0;h<k;h++){
//                 c1+=  A[i*k + h]*B[h*n + j];
//             }
//             C[i*n + j] =c1;
//         }

    
// };

// right version 
__global__ void matcal_3(int m,int n,int k, double *A,double *B, double *C){

    // 2D thread indices defining row and col of element
    int i = (blockIdx.x * blockDim.x + threadIdx.x);
    int j = 2*(blockIdx.y * blockDim.y + threadIdx.y);

    double c1 = 0;
    double c2 = 0;

        if ( i<m && j<n  && (j+1)<m ){
            for (int h= 0;h<k;h++){
                c1 +=  A[i*k + h]*B[h*n + j];

                c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
            } 
            C[i*n + j] =c1;
            C[(i)*n + j+1]=c2;
        }
        else if( i<m && j<n){
            for (int h= 0;h<k;h++){
                c1+=  A[i*k + h]*B[h*n + j];
            }
            C[i*n + j] =c1;
        }

    
};

__global__ void matcal_4(int m,int n,int k, double *A,double *B, double *C){

    // 2D thread indices defining row and col of element
    int i = (blockIdx.x * blockDim.x + threadIdx.x);
    int j = size*(blockIdx.y * blockDim.y + threadIdx.y);

    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    double c4 = 0;
    double c5 = 0;
    double c6 = 0;
    double c7 = 0;
    double c8 = 0;

        if (size ==2){
            if ( i<m && j<n  && (j+1)<m ){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];

                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
            }
            else if( i<m && j<n){
                for (int h= 0;h<k;h++){
                    c1+=  A[i*k + h]*B[h*n + j];
                }
                C[i*n + j] =c1;
            }
        }
        else if(size ==4){
            if ( i<m && j<n  && (j+1)<m && (j+2)<m && (j+3)<m ){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];
                    c4 +=  A[(i)*k + h]*B[h*n + (j+3)];
                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
                C[(i)*n + j+3]=c4;
            }
            else if( i<m && j<n  && (j+1)<m && (j+2)<m){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];

                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
            } 
            else if( i<m && j<n  && (j+1)<m){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];

                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
            } 
            else if( i<m && j<n){
                for (int h= 0;h<k;h++){
                    c1+=  A[i*k + h]*B[h*n + j];
                }
                C[i*n + j] =c1;
            }
        }
        else if(size ==8){
            if (i<m && j<n  && (j+1)<m && (j+2)<m && (j+3)<m && (j+4) && (j+5) && (j+6) && (j+7)){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];
                    c4 +=  A[(i)*k + h]*B[h*n + (j+3)];
                    c5 +=  A[(i)*k + h]*B[h*n + (j+4)];
                    c6 +=  A[(i)*k + h]*B[h*n + (j+5)];
                    c7 +=  A[(i)*k + h]*B[h*n + (j+6)];
                    c8 +=  A[(i)*k + h]*B[h*n + (j+7)];
                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
                C[(i)*n + j+3]=c4;
                C[(i)*n + j+4]=c5;
                C[(i)*n + j+5]=c6;
                C[(i)*n + j+6]=c7;
                C[(i)*n + j+7]=c8;
            }
            else if ( i<m && j<n  && (j+1)<m && (j+2)<m && (j+3)<m && (j+4) && (j+5) && (j+6) ){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];
                    c4 +=  A[(i)*k + h]*B[h*n + (j+3)];
                    c5 +=  A[(i)*k + h]*B[h*n + (j+4)];
                    c6 +=  A[(i)*k + h]*B[h*n + (j+5)];
                    c7 +=  A[(i)*k + h]*B[h*n + (j+6)];
                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
                C[(i)*n + j+3]=c4;
                C[(i)*n + j+4]=c5;
                C[(i)*n + j+5]=c6;
                C[(i)*n + j+6]=c7;
            }
            else if ( i<m && j<n  && (j+1)<m && (j+2)<m && (j+3)<m && (j+4) && (j+5)){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];
                    c4 +=  A[(i)*k + h]*B[h*n + (j+3)];
                    c5 +=  A[(i)*k + h]*B[h*n + (j+4)];
                    c6 +=  A[(i)*k + h]*B[h*n + (j+5)];
                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
                C[(i)*n + j+3]=c4;
                C[(i)*n + j+4]=c5;
                C[(i)*n + j+5]=c6;
            }
            else if ( i<m && j<n  && (j+1)<m && (j+2)<m && (j+3)<m && (j+4)){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];
                    c4 +=  A[(i)*k + h]*B[h*n + (j+3)];
                    c5 +=  A[(i)*k + h]*B[h*n + (j+4)];
                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
                C[(i)*n + j+3]=c4;
                C[(i)*n + j+4]=c5;
            }
            else if ( i<m && j<n  && (j+1)<m && (j+2)<m && (j+3)<m ){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];
                    c4 +=  A[(i)*k + h]*B[h*n + (j+3)];
                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
                C[(i)*n + j+3]=c4;
            }
            else if( i<m && j<n  && (j+1)<m && (j+2)<m){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];
                    c3 +=  A[(i)*k + h]*B[h*n + (j+2)];

                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
                C[(i)*n + j+2]=c3;
            } 
            else if( i<m && j<n  && (j+1)<m){
                for (int h= 0;h<k;h++){
                    c1 +=  A[i*k + h]*B[h*n + j];
                    c2 +=  A[(i)*k + h]*B[h*n + (j+1)];

                } 
                C[i*n + j] =c1;
                C[(i)*n + j+1]=c2;
            } 
            else if( i<m && j<n){
                for (int h= 0;h<k;h++){
                    c1+=  A[i*k + h]*B[h*n + j];
                }
                C[i*n + j] =c1;
            }
        }

    
};

void matmult_lib (int m, int n, int k, double *A, double *B, double *C){



    double alpha = 1.0;
    int beta = 0;
    double ts, te;
    int lda = k;
    int ldb = n;
    int ldc = n;

    ts = omp_get_wtime();
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, m,n, k ,alpha, A, lda, B, ldb, beta, C, ldc);
    te = omp_get_wtime() - ts;
    printf("%lf \n",te);
}

void matmult_gpu1(int m, int n, int k, double *h_A,double *h_B,double *h_C){

    double *d_A;
    double *d_B;
    double *d_C;

    // Allocating memory on the device  
    cudaMalloc((void**)&d_A, m*k * sizeof(double));
    cudaMalloc((void**)&d_B, n*k * sizeof(double));
    cudaMalloc((void**)&d_C, m*n * sizeof(double));

    // Copying input matrix A and B to the device 
    cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);

    // Lounching the kernel 
    matcal_1<<<1,1>>>(m,n,k,d_A,d_B,d_C);
    cudaDeviceSynchronize();

    // Copying the calculated C matrix to host 
    cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

    // Free all allocated memory on the hos 
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}

void matmult_gpu2(int m, int n, int k, double *h_A,double *h_B,double *h_C){

    double *d_A;
    double *d_B;
    double *d_C;


    cudaMalloc((void**)&d_A, m*k * sizeof(double));
    cudaMalloc((void**)&d_B, n*k * sizeof(double));
    cudaMalloc((void**)&d_C, m*n * sizeof(double));

    cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);

    dim3 dimBlock(32,32,1); // threads per block
    dim3 dimGrid(m/32 + 1,n/32 + 1,1);// blocks in total

    matcal_2<<<dimGrid, dimBlock>>>(m,n,k,d_A,d_B,d_C);
    cudaDeviceSynchronize();

    cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}

// Below version
// void matmult_gpu3(int m, int n, int k, double *h_A,double *h_B,double *h_C){

//     double *d_A;
//     double *d_B;
//     double *d_C;


//     cudaMalloc((void**)&d_A, m*k * sizeof(double));
//     cudaMalloc((void**)&d_B, n*k * sizeof(double));
//     cudaMalloc((void**)&d_C, m*n * sizeof(double));

//     cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);

//     dim3 dimBlock(32,32,1); // threads per block
//     dim3 dimGrid((m/(32*2)+ 1),n/32 + 1,1);// blocks in total

//     matcal_3<<<dimGrid, dimBlock>>>(m,n,k,d_A,d_B,d_C);
//     cudaDeviceSynchronize();

//     cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

//     cudaFree(d_A);
//     cudaFree(d_B);
//     cudaFree(d_C);
// }

// right version
void matmult_gpu3(int m, int n, int k, double *h_A,double *h_B,double *h_C){

    double *d_A;
    double *d_B;
    double *d_C;


    cudaMalloc((void**)&d_A, m*k * sizeof(double));
    cudaMalloc((void**)&d_B, n*k * sizeof(double));
    cudaMalloc((void**)&d_C, m*n * sizeof(double));

    cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);

    dim3 dimBlock(32,32,1); // threads per block
    dim3 dimGrid(((m/32)+ 1),n/(32*2) + 1,1);// blocks in total

    matcal_3<<<dimGrid, dimBlock>>>(m,n,k,d_A,d_B,d_C);
    cudaDeviceSynchronize();

    cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}

void matmult_gpu4(int m, int n, int k, double *h_A,double *h_B,double *h_C){

    double *d_A;
    double *d_B;
    double *d_C;


    cudaMalloc((void**)&d_A, m*k * sizeof(double));
    cudaMalloc((void**)&d_B, n*k * sizeof(double));
    cudaMalloc((void**)&d_C, m*n * sizeof(double));

    cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);

    dim3 dimBlock(32,32,1); // threads per block
    dim3 dimGrid(((m/32)+ 1),n/(32*size) + 1,1);// blocks in total

    matcal_4<<<dimGrid, dimBlock>>>(m,n,k,d_A,d_B,d_C);
    cudaDeviceSynchronize();

    cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}


}
