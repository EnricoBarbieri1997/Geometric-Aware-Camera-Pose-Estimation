#!/bin/bash

julia -e 'using Pkg; Pkg.activate("./"); using CylindersBasedCameraResectioning; CylindersBasedCameraResectioning.Report.multiple_seeds_multiple_configuration()'