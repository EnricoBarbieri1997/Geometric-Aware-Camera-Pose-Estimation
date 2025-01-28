#!/bin/bash

for i in {1..5}; do
    for j in {1..10}; do
        nohup ./run.sh &> ./tmp/output_seed_$i_noise_$j.log &
    done
done