#!/bin/bash

for i in {1..5}
    for j in {1..10}
        nohup ./run.sh &> ./tmp/output_seed_$i_noise_$j.log &
    end
end