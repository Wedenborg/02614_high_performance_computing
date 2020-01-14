#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif



void main(void){

    #pragma omp parallel default(none) \shared(N,pi) private(i,inte,fx,prod)
    {
    
    int N =1000000000;
    double inte =0;
    double prod =0;
    double pi=0;
    double fx=0;
    #pragma omp for
    for (int i=0; i<N;i++){
        prod = ((i-0.5)/N)*((i-0.5)/N);
        inte =4/(1+prod);
        fx = inte;
    #pragma omp critical
        pi += fx/N;
    }
    

    printf("%lf /N", pi);

    }

} // end 



