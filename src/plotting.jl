module Plotting
    export initfigure, add_2d_axis!, add_slider!, clean_plots!, get_or_add_2d_axis!, get_or_add_camera_rotation_axis!, plot_2dpoints, plot_line_2d, Plot3dCameraInput, plot_3dcamera, plot_3dcamera_rotation, plot_cylinders_contours, plot_3dcylinders, plot_2dcylinders, plot_image_background, save_2d_figures, Figure, Axis, lift, lines!, image!, on, events, Point2f, Observable, ispressed, Mouse

    using ..CylindersBasedCameraResectioning: GUI_ENABLED

    colors = [:red, :green, :blue, :yellow, :purple, :orange, :pink, :brown]
    linestyles = Dict(
        "dash" => :dash,
        "solid" => :solid,
    )

    struct Plot3dCameraInput
        cameraRotation::Vector{<:Number}
        cameraTranslation::Vector{<:Number}
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
    function plot_3dcamera(camera::Any, color=:black)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end
    function plot_3dcamera_rotation(camera::Any; color=:black, axindex=nothing)
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
    function plot_image_background(image; axindex=1)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end

    function save_2d_figures(path, scene, problems)
        # Base.showerror(stdout, error("Plotting not available in this environment."))
    end

    if GUI_ENABLED && Base.find_package("GLMakie") !== nothing
        include("plotting-real.jl")
    end

    const Figure = if GUI_ENABLED
        Makie.Figure
    else
        function (; kwargs...)
            # error("Plotting not available in this environment.")
        end
    end

    const Axis = if GUI_ENABLED
        Makie.Axis
    else
        function (; kwargs...)
            # error("Plotting not available in this environment.")
        end
    end

    const lift = if GUI_ENABLED
        Makie.lift
    else
        function (v::Any)
            # error("Plotting not available in this environment.")
        end
    end

    const lines! = if GUI_ENABLED
        Makie.lines!
    else
        function (ax::Any, x::Any, y::Any; kwargs...)
            # error("Plotting not available in this environment.")
        end
    end

    const image! = if GUI_ENABLED
        Makie.image!
    else
        function (ax::Any, img::Any; kwargs...)
            # error("Plotting not available in this environment.")
        end
    end

    const on = if GUI_ENABLED
        Makie.on
    else
        function (ax::Any, event::Any; kwargs...)
            # error("Plotting not available in this environment.")
        end
    end

    const events = if GUI_ENABLED
        Makie.events
    else
        function (ax::Any, event::Any; kwargs...)
            # error("Plotting not available in this environment.")
        end
    end

    const Point2f = if GUI_ENABLED
        Makie.Point2f
    else
        function (x::Float64, y::Float64)
            # error("Plotting not available in this environment.")
        end
    end

    const Observable = if GUI_ENABLED
        Makie.Observable
    else
        function (x::Any)
            # error("Plotting not available in this environment.")
        end
    end

    const ispressed = if GUI_ENABLED
        Makie.ispressed
    else
        function (x::Any)
            # error("Plotting not available in this environment.")
        end
    end
    
    const Mouse = if GUI_ENABLED
        Makie.Mouse
    else
        function (x::Any)
            # error("Plotting not available in this environment.")
        end
    end
end