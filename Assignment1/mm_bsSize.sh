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
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 00:59

# define the driver name to use
# valid values: matmult_c.studio, matmult_f.studio, matmult_c.gcc or
# matmult_f.gcc
#
EXECUTABLE=matmult_c.gcc

# define the mkn values in the MKN variable
#
SIZES="2000"

# define the permutation type in PERM
#
PERM="blk"

# uncomment and set a reasonable BLKSIZE for the blk version
#
BLKSIZE="10 15 20 30 40 50 65 80 90 120 150 180 250 300 400 500 600 700 800 900 1000 1200 1500 1800 2000"

# enable(1)/disable(0) result checking
export MATMULT_COMPARE=0

# start the collect command with the above settings
for B in $BLKSIZE
do
    ./$EXECUTABLE $PERM $SIZES $SIZES $SIZES $B
done
