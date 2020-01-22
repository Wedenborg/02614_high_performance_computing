#!/bin/bash
# 02614 - High-Performance Computing, January 2018
# 
# batch script to run collect on a decidated server in the hpcintro
# queue
#
# Author: Bernd Dammann <bd@cc.dtu.dk>
#
#BSUB -J Jacobi_gpu1
#BSUB -o jacobi_gpu1_%J
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process:mps=yes"
#BSUB -R "rusage[mem=7GB]"
#BSUB -W 00:59


module load gcc/8.3.0
module load cuda/10.2

# define the driver name to use
# valid values: matmult_c.studio, matmult_f.studio, matmult_c.gcc or
# matmult_f.gcc
#
EXECUTABLE=jacobi_gpu1

# define the mkn values in the MKN variable
#
Size="20 30 40"

# define the permutation type in PERM
#
max_iter=1000

# uncomment and set a reasonable BLKSIZE for the blk version
#

T_time=15

# define the max no. of iterations the driver should use - adjust to
# get a reasonable run time.  You can get an estimate by trying this
# on the command line, i.e. "MFLOPS_MAX_IT=10 ./matmult_...." for the
# problem size you want to analyze.
#
export MFLOPS_MAX_IT=50000

# experiment name 
#
JID=${LSB_JOBID}
EXPOUT="$LSB_JOBNAME.${JID}.er"

# uncomment the HWCOUNT line, if you want to use hardware counters
# define an option string for the harwdware counters (see output of
# 'collect -h' for valid values.  The format is:
# -h cnt1,on,cnt2,on,...  (up to four counters at a time)
#
# the example below is for L1 hits, L1 misses, L2 hits, L2 misses
#
HWCOUNT="-h dch,on,dcm,on,l2h,on,l2m,on,l3h,on,l3m,on"
tol=10

# start the collect command with the above settings
for S in $Size
do
    nvprof --print-gpu-summary ./$EXECUTABLE $S $max_iter $tol $T_time 
done


##collect -o $EXPOUT $HWCOUNT ./$EXECUTABLE $S $ $BLKSIZE
