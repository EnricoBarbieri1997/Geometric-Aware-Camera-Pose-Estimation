#!/bin/bash

seed=$1
noise=$2

command_common='using Pkg; Pkg.activate("./"); using CylindersBasedCameraResectioning; CylindersBasedCameraResectioning.Report.multiple_seeds_multiple_configuration'

if [ -z "$seed" ] && [ -z "$noise" ]; then
    julia -e "${command_common}();"
else
    julia -e "${command_common}(;seed_index=$seed, noises=[$noise]);"
fi
