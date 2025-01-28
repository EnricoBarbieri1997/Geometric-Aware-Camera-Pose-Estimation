#!/bin/bash

for i in {1..5}; do
    for j in {1..10}; do
        pwd
        nohup ./scripts/run.sh "$i" "$j" &> ./tmp/output_seed_${i}_noise_${j}.log &
    done
done