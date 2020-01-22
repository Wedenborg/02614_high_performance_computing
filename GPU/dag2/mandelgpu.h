#ifndef __MANDEL_H
#define __MANDEL_H

__global__ void mandelgpu(int width, int height, int *image, int max_iter);

#endif
