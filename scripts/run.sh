#!/bin/bash

seed=$1
noise=$2
save_in_folder=$3

command_common='using Pkg; Pkg.activate("./"); using CylindersBasedCameraResectioning; CylindersBasedCameraResectioning.Report.multiple_seeds_multiple_configuration'

if [ -z "$seed" ] && [ -z "$noise" ] && [ -z "$save_in_folder" ]; then
    julia -e "${command_common}();"
else
    julia -e "${command_common}(;seed_index=$seed, noises=[$noise], save_in_folder=$save_in_folder);"
fi
