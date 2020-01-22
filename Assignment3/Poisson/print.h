#ifndef _PRINT_H
#define _PRINT_H

__global__
void print_binary(char *fname, int num, double ***u);

__global__
void print_vtk(const char *fname, int n, double ***u);

#endif
