#include <cuda_runtime_api.h>
#include <helper_cuda.h>

void
transfer_3d(double ***dst, double ***src, int m, int n, int k, int flag)
{
    long nPtr = m + m * n;
    long nBlk = m * n * k;

    // we only transfer the value block
    checkCudaErrors( cudaMemcpy((double *) dst + nPtr,
                                (double *) src + nPtr,
                                nBlk * sizeof(double),
                                (cudaMemcpyKind) flag) );
}
