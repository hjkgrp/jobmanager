#!/bin/bash

#SBATCH -o log.o%j
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=96:00:00

# user environment/variables
export ORCA_IN=$1
export ORCA_OUT=$2

# set up calculation environment/variables
module load orca
export SUBDIR=$PWD

# ignore no CUDA and no internet warnings from MPI
export PARAMS="--mca opal_warn_on_missing_libcuda 0 --mca btl ^openib"

echo "Temporary directory: "$TMPDIR
cd $TMPDIR
cp $SUBDIR/$ORCA_IN ./
cp $SUBDIR/*.xyz ./
cp $SUBDIR/*.gbw ./

$ORCA_PATH $ORCA_IN "$PARAMS" > $SUBDIR/$ORCA_OUT

cp -r $TMPDIR/* $SUBDIR/
