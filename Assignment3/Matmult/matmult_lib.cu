#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
extern "C"{
#ifdef _OPENMP
#include <omp.h>
#endif
#include "cblas.h"

#define size 8 // gpu4

// Thread block size
#define BLOCK_SIZE 16 // gpu5

// Matrices are stored in row-major order:
// M(row, col) = *(M.elements + row * M.stride + col)
typedef struct {
    int width;
    int height;
    int stride; 
    double* elements;
} Matrix;

// Get a matrix element
__device__ double GetElement(const Matrix A, int row, int col)
{
    return A.elements[row * A.stride + col];
}

// Set a matrix element
__device__ void SetElement(Matrix A, int row, int col,double value)
{
    A.elements[row * A.stride + col] = value;
}

// Get the BLOCK_SIZExBLOCK_SIZE sub-matrix Asub of A that is
// located col sub-matrices to the right and row sub-matrices down
// from the upper-left corner of A
 __device__ Matrix GetSubMatrix(Matrix A, int row, int col) 
{
    Matrix Asub;
    Asub.width    = BLOCK_SIZE;
    Asub.height   = BLOCK_SIZE;
    Asub.stride   = A.stride;
    Asub.elements = &A.elements[A.stride * BLOCK_SIZE * row
                                         + BLOCK_SIZE * col];
    return Asub;
}

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
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Col index in C
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Row index in C

    if ( i<m && j<n ){

        double c1 = 0;
        for (int h= 0;h<k;h++){
            c1 +=  A[i*k + h]*B[h*n + j];
        } 
        C[i*n + j] = c1;
    }
};

// below version
__global__ void matcal_3(int m,int n,int k, double *A,double *B, double *C){

    // 2D thread indices defining row and col of element
    int j = (blockIdx.x * blockDim.x + threadIdx.x);
    int i = 2*(blockIdx.y * blockDim.y + threadIdx.y);

    double c1 = 0;
    double c2 = 0;

        if ( i<m && j<n  && (i+1)<m ){
            for (int h= 0;h<k;h++){
                c1 +=  A[i*k + h]*B[h*n + j];

                c2 +=  A[(i+1)*k + h]*B[h*n + j];
            } 
            C[i*n + j] =c1;
            C[(i+1)*n + j]=c2;
        }
        else if( i<m && j<n){
            for (int h= 0;h<k;h++){
                c1+=  A[i*k + h]*B[h*n + j];
            }
            C[i*n + j] =c1;
        }

    
// };

// right version 
// __global__ void matcal_3(int m,int n,int k, double *A,double *B, double *C){

//     // 2D thread indices defining row and col of element
//     int j = 2*(blockIdx.x * blockDim.x + threadIdx.x);
//     int i = (blockIdx.y * blockDim.y + threadIdx.y);

//     double c1 = 0;
//     double c2 = 0;

//         if ( i<m && j<n  && (j+1)<n ){
//             for (int h= 0;h<k;h++){
//                 c1 +=  A[i*k + h]*B[h*n + j];
//                 c2 +=  A[i*k + h]*B[h*n + (j+1)];
//             } 
//             C[i*n + j] =c1;
//             C[i*n + j+1]=c2;
//         }
//         else if( i<m && j<n){
//             for (int h= 0;h<k;h++){
//                 c1+=  A[i*k + h]*B[h*n + j];
//             }
//             C[i*n + j] =c1;
//         }

    
// };

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

// Matrix multiplication kernel called by MatMul()
__global__ void MatMulKernel(Matrix A, Matrix B, Matrix C)
{
    // Block row and column
    int blockRow = blockIdx.y;
    int blockCol = blockIdx.x;

    // Each thread block computes one sub-matrix Csub of C
    Matrix Csub = GetSubMatrix(C, blockRow, blockCol);

    // Each thread computes one element of Csub
    // by accumulating results into Cvalue
    double Cvalue = 0;

    // Thread row and column within Csub
    int row = threadIdx.y;
    int col = threadIdx.x;

    // Loop over all the sub-matrices of A and B that are
    // required to compute Csub
    // Multiply each pair of sub-matrices together
    // and accumulate the results
    for (int m = 0; m < (A.width / BLOCK_SIZE); ++m) {

        // Get sub-matrix Asub of A
        Matrix Asub = GetSubMatrix(A, blockRow, m);

        // Get sub-matrix Bsub of B
        Matrix Bsub = GetSubMatrix(B, m, blockCol);

        // Shared memory used to store Asub and Bsub respectively
        __shared__ double As[BLOCK_SIZE][BLOCK_SIZE];
        __shared__ double Bs[BLOCK_SIZE][BLOCK_SIZE];

        // Load Asub and Bsub from device memory to shared memory
        // Each thread loads one element of each sub-matrix
        As[row][col] = GetElement(Asub, row, col);
        Bs[row][col] = GetElement(Bsub, row, col);

        // Synchronize to make sure the sub-matrices are loaded
        // before starting the computation
        __syncthreads();
        // Multiply Asub and Bsub together
        for (int e = 0; e < BLOCK_SIZE; ++e)
            Cvalue += As[row][e] * Bs[e][col];

        // Synchronize to make sure that the preceding
        // computation is done before loading two new
        // sub-matrices of A and B in the next iteration
        __syncthreads();
    }

    // Write Csub to device memory
    // Each thread writes one element
    SetElement(Csub, row, col, Cvalue);
}
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

void matmult_gpulib (int m, int n, int k, double *h_A, double *h_B, double *h_C){

    
    const double alf = 1;
    const double bet = 0;
    const double *alpha = &alf;
    const double *beta = &bet;
    
    double ts, te;
    int lda = k;
    int ldb = n; // .
    int ldc = n;
    double *d_A;
    double *d_B;
    double *d_C;


    // cublasDgemm(cublasHandle_t handle,
    //     cublasOperation_t transa, cublasOperation_t transb,
    //     int m, int n, int k,
    //     const double          *alpha,
    //     const double          *A, int lda,
    //     const double          *B, int ldb,
    //     const double          *beta,
    //     double          *C, int ldc)

    // Allocating memory on the device  
    cudaMalloc((void**)&d_A, m*k * sizeof(double));
    cudaMalloc((void**)&d_B, n*k * sizeof(double));
    cudaMalloc((void**)&d_C, m*n * sizeof(double));

    // Copying input matrix A and B to the device 
    cudaMemcpy(d_A, h_A,  m*k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B,  n*k * sizeof(double), cudaMemcpyHostToDevice);
    
    // Create a handle for CUBLAS
    cublasHandle_t handle;
    cublasCreate(&handle);
    
    // ts = omp_get_wtime();
    cublasDgemm(handle,CUBLAS_OP_N,CUBLAS_OP_N, n,m, k ,alpha, d_B, ldb, d_A, lda, beta, d_C, ldc);
    // te = omp_get_wtime() - ts;
    // printf("%lf \n",te);    // Copying the calculated C matrix to host 
    
    // Destroy the handle
    cublasDestroy(handle);
    cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

    // Free all allocated memory on the hos 
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
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
    dim3 dimGrid(n/32 + 1,m/32 + 1,1);// blocks in total

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
    dim3 dimGrid((n/(32*2)+ 1),m/(32) + 1,1);// blocks in total

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
    dim3 dimGrid(((n/32)+ 1),m/(32*size) + 1,1);// blocks in total

    matcal_4<<<dimGrid, dimBlock>>>(m,n,k,d_A,d_B,d_C);
    cudaDeviceSynchronize();

    cudaMemcpy(h_C, d_C, m*n * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}

void matmult_gpu5(int mm, int nn, int kk, double *h_A,double *h_B,double *h_C){

    // Load A and B to device memory
    Matrix d_A;
    d_A.width = kk; 
    d_A.stride = kk; 
    d_A.height = mm;
    // size = kk * mm * sizeof(double);
    cudaMalloc(&d_A.elements, kk * mm * sizeof(double));
    cudaMemcpy(d_A.elements, h_A, kk * mm * sizeof(double),cudaMemcpyHostToDevice);
    Matrix d_B;
    d_B.width = nn;
    d_B.stride = nn; 
    d_B.height = kk;
    // size = nn * kk * sizeof(double);
    cudaMalloc(&d_B.elements,  nn * kk * sizeof(double));
    cudaMemcpy(d_B.elements, h_B, nn * kk * sizeof(double),cudaMemcpyHostToDevice);

    // Allocate C in device memory
    Matrix d_C;
    d_C.width = nn; 
    d_C.stride = nn; 
    d_C.height = mm;
    // size = nn * mm * sizeof(double);
    cudaMalloc(&d_C.elements,  nn * mm * sizeof(double));

    // Invoke kernel
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(nn / dimBlock.x, mm / dimBlock.y);
    MatMulKernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_C);

    // Read C from device memory
    cudaMemcpy(h_C, d_C.elements,  nn * mm * sizeof(double),cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_A.elements);
    cudaFree(d_B.elements);
    cudaFree(d_C.elements);
}

}
