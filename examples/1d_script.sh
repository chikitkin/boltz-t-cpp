#!/bin/bash
# export OMP_NUM_THREADS=$nt
mkdir RESULTS_1D
cd RESULTS_1D

mkdir 1d_40
cd 1d_40
../../../src/cli.o ../../1d_40/1d.cfg
python3 ../../plot.py T.txt
python3 ../../plot.py res.txt
cd ..

mkdir 1d_80
cd 1d_80
../../../src/cli.o ../../1d_80/1d.cfg
python3 ../../plot.py T.txt
python3 ../../plot.py res.txt
cd ..

mkdir 1d_160
cd 1d_160
../../../src/cli.o ../../1d_160/1d.cfg
python3 ../../plot.py T.txt
python3 ../../plot.py res.txt
cd ..


# gnuplot -persist <<-EOFMarker
#     set output 'n.png'
#     plot "T.txt" using 1:2 with linespoints
#     set output 'ux.png'
#     plot "T.txt" using 1:3 with linespoints
#     set output 'T.png'
#     plot "T.txt" using 1:4 with linespoints
#     set output 'res.png'
#     plot "res.txt" with linespoints
# EOFMarker
