#!/bin/bash
for (( nt=1; nt<5; nt++ ))
do
export OMP_NUM_THREADS=$nt
mkdir $nt
cd $nt
../../../src/cli.o ../1d.cfg
cd ..
done
