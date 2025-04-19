module Plotting
    export initfigure, add_2d_axis!, add_slider!, clean_plots!, get_or_add_2d_axis!, get_or_add_camera_rotation_axis!, plot_2dpoints, plot_line_2d, Plot3dCameraInput, plot_3dcamera, plot_3dcamera_rotation, Plot3dCylindersInput, plot_cylinders_contours, plot_3dcylinders, plot_2dcylinders

    colors = [:red, :green, :blue, :yellow, :purple, :orange, :pink, :brown]
    linestyles = Dict(
        "dash" => :dash,
        "solid" => :solid,
    )

    struct Plot3dCameraInput
        cameraRotation::Vector{<:Number}
        cameraTranslation::Vector{<:Number}
    end

    struct Plot3dCylindersInput
        transforms::Vector{Matrix{Float64}}
        radiuses::Vector{Vector{Float64}}
        numberOfCylinders::Int
        cameraProjectionMatrix::Union{Matrix{<:Number}, UndefInitializer}

        function Plot3dCylindersInput(transforms::Vector{Matrix{Float64}}, radiuses::Vector{Vector{Float64}}, numberOfCylinders::Int, cameraProjectionMatrix::Union{Matrix{<:Number}, UndefInitializer} = undef)
            new(transforms, radiuses, numberOfCylinders, cameraProjectionMatrix)
        end
    end

    """ Dummy implementations (no-op or error) for environments without plotting support """
    function initfigure()
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function add_2d_axis!()
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function add_slider!(; start=0.0, stop=1.0, step=0.01)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function get_or_add_2d_axis!(index)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function get_or_add_camera_rotation_axis!(index)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function clean_plots!()
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_2dpoints(points; axindex=1)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_line_2d(line; color=:black, linestyle=:solid, axindex=1)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_3dcamera(info::Any, color=:black)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_3dcamera_rotation(info::Any; color=:black, axindex=nothing)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_3dcylinders(cylindersInfo::Any; axindex=1)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_cylinders_contours(contours::Any; linestyle=:solid)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_2dcylinders(conic_contours; linestyle=:solid, alpha=1, axindex=1)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end

    const ENABLE_PLOTTING = get(ENV, "ENABLE_PLOTTING", "false") == "true"

    if ENABLE_PLOTTING && Base.find_package("GLMakie") !== nothing
        include("plotting-real.jl")
    end
end