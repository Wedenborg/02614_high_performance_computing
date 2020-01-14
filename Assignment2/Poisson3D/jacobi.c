/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>

void
jacobi(int N, double ***u) {

    double ***uu = NULL;

    if ( (uu = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);

    
}
