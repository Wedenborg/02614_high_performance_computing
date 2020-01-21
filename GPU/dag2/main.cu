#include <stdio.h>
#include <stdlib.h>
#include "mandelgpu.h"
#include "writepng.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <helper_cuda.h>




int
main(int argc, char *argv[]) {

    int   width, height;
    int	  max_iter;
    int   *image;
    double ts, te;

    width    = 2601;
    height   = 2601;
    max_iter = 400;

    // command line argument sets the dimensions of the image
    if ( argc == 2 ) width = height = atoi(argv[1]);

    // image = (int *)malloc( width * height * sizeof(int));
    // if ( image == NULL ) {
    //    fprintf(stderr, "memory allocation failed!\n");
    //    return(1);
    // }
    cudaMallocManaged((void **)&image, width * height * sizeof(int));
    dim3 dimBlock(17,17,1); // threads per block
    dim3 dimGrid(153,153,1);// blocks in total
    ts = omp_get_wtime();
    // mandelgpu<<<dim3(5,5),dim3(width/5,height/5)>>>(width, height, image, max_iter);
    mandelgpu<<< dimGrid, dimBlock >>>(width, height, image, max_iter);
    checkCudaErrors(cudaDeviceSynchronize());
    te = omp_get_wtime() - ts;
    printf("Time :%lf \n",te);
    writepng("mandelbrot.png", image, width, height);

    cudaFree(image);

    return(0);
}
