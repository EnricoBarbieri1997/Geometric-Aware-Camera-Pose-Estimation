#!/bin/bash

for i in {1..5}; do
    for j in {0..10}; do
        nohup ./scripts/run.sh "$i" "$j" &> ./tmp/output_seed_${i}_noise_${j}.log &
    done
done