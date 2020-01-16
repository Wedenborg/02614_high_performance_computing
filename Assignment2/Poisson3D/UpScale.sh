#!/bin/bash
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -W 15
#BSUB -J bench_poisson
#BSUB -o bench_poisson_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here (uncomment and adjust compiler and version!)

module load gcc

THREADS="1 2 4 8 12 16 20 24"
CMD=poisson_j
N=1
IT=2000
TOL=0.0005
TS=15

for t in $THREADS
do
    OMP_NUM_THREADS=$t ./$CMD $t $IT $TOL $TS 
done