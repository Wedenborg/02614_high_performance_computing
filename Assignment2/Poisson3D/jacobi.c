/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>

void
jacobi() {
    // fill in your code here

  double u_new = 0;
  double u_old = 0;
  double cond =0.005;
  double stopTest = 100000;
  double stopTest =0;
  for(int i=0;i<N;i++){
    for(int j=0;j<M;j++){
      for(int k=0;k<L;k++){
	  while( sqrt(stopTest)<cond){
	    stopTest +=(u_new[i][j][k]-u_old[i][j][k])*(u_new[i][j][k]-u_old[i][j][k])
	    u_old = u_new
	    u_new = 4 //formula and matrix


	     }
      }

    }
  }

  
  
}
