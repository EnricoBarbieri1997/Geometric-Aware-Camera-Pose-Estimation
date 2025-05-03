#!/bin/bash

seed=$1
noise=$2
output=$3

command_common='using Pkg; Pkg.activate("./"); using CylindersBasedCameraResectioning; CylindersBasedCameraResectioning.Report.multiple_seeds_multiple_configuration'

if [ -z "$seed" ] && [ -z "$noise" ] && [ -z "$output" ]; then
    julia -e "${command_common}(;output='./tmp/reports/calibration/');"
else
    julia -e "${command_common}(;seed_index=$seed, noises=[$noise], output=$output);"
fi
