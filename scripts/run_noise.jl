using Pkg
Pkg.activate(".")

using CylindersBasedCameraResectioning

noise = parse(Float64, ARGS[1])
@info "Running with noise=$noise"

CylindersBasedCameraResectioning.Report.multiple_seeds_multiple_configuration(;
    noises=[noise],
    output=true,
)
