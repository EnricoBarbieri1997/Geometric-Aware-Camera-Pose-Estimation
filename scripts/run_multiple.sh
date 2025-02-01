#!/bin/bash

noises=($(seq 0 0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0))
for i in {1..1}; do
    for j in $noises; do
        nohup ./scripts/run.sh "$i" "$j" true &> ./tmp/output_seed_${i}_noise_${j}.log &
    done
done