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

    for( int i =0; i < N; i++){
        for( int j = 0; j < N; j++){
            for( int k = 0; k < N; k++){
                uu[i][j][k] = u[i][j][k];
            }
        }
    }
    
}
