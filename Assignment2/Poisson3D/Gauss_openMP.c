/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef _OMP
#include<omp.h>
#endif

void Gauss_openMP(int N, double ***u, double ***v, double ***f, int iter_max, double tolerance) {
    int i,j,k;
    int static t;
    #pragma omp threadprivate(t)
    t=0;

    #pragma omp parallel default(none) shared(u, N, tolerance,f,iter_max) private(i, j, k) copyin(t)
    {
    do {
        #pragma omp for ordered(2) // shared(u, N, tolerance,f) private(i, j, k) copyin(t)
        for(i =1; i < N-1; i++){
            for(j = 1; j < N-1; j++){
                #pragma omp ordered depend(sink:i-1,j) depend(sink: i,j-1)
                for( k = 1; k < N-1; k++){
	                u[i][j][k] = 1./6.*(u[i-1][j][k]+u[i+1][j][k]+u[i][j-1][k]+u[i][j+1][k]+u[i][j][k-1]+u[i][j][k+1] + 1./((N-2)*(N-2)) * f[i][j][k]); //formula and matrix      
                }
                #pragma omp ordered depend(source)
                 
            }
            
        }
        t++;
    } while (t< iter_max);
    
}
}
