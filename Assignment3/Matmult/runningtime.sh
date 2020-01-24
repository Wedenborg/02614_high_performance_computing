#!/bin/bash
# 02614 - High-Performance Computing, January 2018
# 
# batch script to run matmult on a decidated server in the hpcintro
# queue
#
# Author: Bernd Dammann <bd@cc.dtu.dk>
#
#BSUB -J mm_batch
#BSUB -o mm_batch_%J.out
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process:mps=yes"
#BSUB -R "rusage[mem=7GB]"
#BSUB -W 5

module load cuda/10.2
module load gcc/6.3.0

# define the driver name to use
# valid values: matmult_c.studio, matmult_f.studio, matmult_c.gcc or
# matmult_f.gcc
#
EXECUTABLE=matmult_f.nvcc

# define the mkn values in the MKN variable
#
SIZES="500 1000 2000 4000 8000"

# define the permutation type in PERM
#
PERM="gpu4"

# uncomment and set a reasonable BLKSIZE for the blk version
#
# BLKSIZE=1

# enable(1)/disable(0) result checking
export MATMULT_COMPARE=0

# start the collect command with the above settings
for S in $SIZES
do
    MATMULT_COMPARE=0 nvprof --print-gpu-summary ./$EXECUTABLE $PERM $S $S $S
    #ATMULT_COMPARE=0 ./$EXECUTABLE $PERM $S $S $S
done
