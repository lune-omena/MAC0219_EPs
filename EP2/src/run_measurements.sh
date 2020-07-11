#! /bin/bash

set -o xtrace

MEASUREMENTS=15
ITERATIONS=5
INITIAL_SIZE=4096

SIZE=$INITIAL_SIZE

NAMES=('mandelbrot_cuda')

GRID=2;


make
mkdir results

for NAME in ${NAMES[@]}; do
    mkdir results/$NAME

    for ((i=1; i<=$ITERATIONS; i++)); do
            perf stat -r $MEASUREMENTS ./$NAME -0.188 -0.012 0.554 0.754 $SIZE 1 i >> triple_spiral.log 2>&1
            GRID=2*$GRID
    done

    SIZE=$INITIAL_SIZE

    mv *.log results/$NAME
    rm output.ppm
done
