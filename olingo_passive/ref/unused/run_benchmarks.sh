#!/bin/bash

mkdir -p benchmarks

for i in $(seq 1 1000); do
    echo "Running iteration $i..."
    ./test/ghkss256 > benchmarks/${i}.txt
done

echo "Done. Output stored in ./benchmarks/"
