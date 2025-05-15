module Printing
    export print_camera_differences, print_error_analysis, print_relative_motion_errors, create_single_noise_result, save_results_to_json

    using ..Utils: normalized_diff, vector_difference, matrix_difference, rotations_difference, translations_difference, intrinsic_difference
    using PrettyTables, JSON
    using Rotations: params as rotations_params
    using LinearAlgebra: norm

    transparent_first_col = Highlighter(
            (data, i, j) -> (j == 1),
            Crayon(
                foreground = :dark_gray,
            )
        );
    low_value_good = Highlighter(
            (data, i, j) -> data[i, j] < 0.01,
            Crayon(
                foreground = :green,
                bold = :true
            )
        );
    high_value_bad = Highlighter(
            (data, i, j) -> data[i, j] > 0.1,
            Crayon(
                foreground = :red,
                bold = :true
            )
        );
    
    function comparisondata(data1, data2)
        return vcat(vec(data1)', vec(data2)')
    end
    function originalcalculateddata(original, calculated)
        return hcat([
            "original",
            "calculated",
        ], comparisondata(original, calculated))
    end
    function pretty_table_withdefaults(data; header, highlighters = (transparent_first_col))
        return pretty_table(
            data;
            formatters    = ft_printf("%5.2f", 2:4),
            header        = header,
            header_crayon = crayon"yellow bold",
            text_crayon   = crayon"white",
            highlighters  = highlighters,
            tf            = tf_unicode_rounded
        )
    end

    function print_camera_differences(original_camera, calculated_camera; verbose = false)
        header = (
            ["Intrinsic", "f_x", "f_y", "skew", "c_x", "c_y"],
        )
        intrinsic_indexes = [1, 5, 4, 7, 8]
        original_intrinsic = original_camera.intrinsic ./ original_camera.intrinsic[3, 3]
        calculated_intrinsic = calculated_camera.intrinsic ./ calculated_camera.intrinsic[3, 3]
        pretty_table_withdefaults(
            originalcalculateddata(
                original_intrinsic[intrinsic_indexes],
                calculated_intrinsic[intrinsic_indexes]
            );
            header = header,
        )

        header = (
            ["Rotation (Quat)", "w", "x", "y", "z"],
        )
        pretty_table_withdefaults(
            originalcalculateddata(
                rotations_params(original_camera.quaternion_rotation),
                rotations_params(calculated_camera.quaternion_rotation),
            );
            header = header,
        )

        header = (
            ["Rotation (grad)", "x", "y", "z"],
        )
        pretty_table_withdefaults(
            originalcalculateddata(
                original_camera.euler_rotation,
                calculated_camera.euler_rotation
            );
            header = header,
        )

        header = (
            ["Translation", "x", "y", "z"],
        )
        pretty_table_withdefaults(
            originalcalculateddata(
                original_camera.position,
                calculated_camera.position
            );
            header = header,
        )

        header = (
            ["Parameters", "Error"],
        )
        data = [
            "Projection matrix" matrix_difference(
                calculated_camera.matrix,
                original_camera.matrix
            );
            "Quaternion" vector_difference(
                convert(Vector{Float64}, rotations_params(calculated_camera.quaternion_rotation)),
                convert(Vector{Float64}, rotations_params(original_camera.quaternion_rotation)),
            );
            "Intrinsic" matrix_difference(
                calculated_camera.intrinsic,
                original_camera.intrinsic
            );
            "Rotation" rotations_difference(
                calculated_camera.quaternion_rotation,
                original_camera.quaternion_rotation
            );
            "Translation" translations_difference(
                calculated_camera.position,
                original_camera.position
            ) / norm(original_camera.position);
        ]
        pretty_table_withdefaults(data;
            header = header,
            highlighters = (transparent_first_col, low_value_good, high_value_bad)
        )

        header = (
            ["Metric", "Error"],
        )
        data = hcat([
            "Δf",
            "Δc",
            "Δskew",
        ], intrinsic_difference(
            calculated_camera.intrinsic,
            original_camera.intrinsic
        ))
        pretty_table_withdefaults(data;
            header = header,
            highlighters = (transparent_first_col, low_value_good, high_value_bad)
        )

        if verbose
            header = (
                ["Projection", "$j$i" for i in 1:4, j in 1:4],
            )
            pretty_table_withdefaults(
                originalcalculateddata(
                    original_camera.matrix,
                    calculated_camera.matrix
                );
                header = header,
            )
        end
    end

    function print_error_analysis(errors::Matrix{<:Number}; header=nothing, noise_steps=(0.0:1.0:10.0))
        header = vcat(["Metric/Noise"], something(header, noise_steps))
        data_rows = [
            "Δf",
            "Δc",
            "Δskew",
            "Camera matrix",
        ]
        if (size(errors)[1] > 4)
            data_rows = vcat(data_rows, [
                "ΔR",
                "ΔT",
            ])
        end
        data = hcat(data_rows, errors)
        pretty_table_withdefaults(data;
            header = header,
            highlighters = (transparent_first_col, low_value_good, high_value_bad)
        )
    end

    function print_relative_motion_errors(scene, problems)
        views_header = ["Views $(i) - $(i+1)" for i in 1:(size(scene.instances, 1)-1)]
        header = vcat(["Metric/View pair"], views_header)
        data_rows = [
            "ΔR_gt",
            "ΔR",
            "ΔR diff",
            "ΔT_gt",
            "ΔT",
            "ΔT diff",
        ]
        errors = Matrix{Float64}(undef, length(data_rows), length(views_header))
        for i in 1:(length(scene.instances)-1)
            view1 = scene.instances[i]
            view2 = scene.instances[i + 1]
            errors[1, i] = rotations_difference(
                view1.camera.quaternion_rotation,
                view2.camera.quaternion_rotation
            )
            errors[4, i] = norm(view2.camera.position - view1.camera.position)

            problem1 = problems[i]
            problem2 = problems[i + 1]
            errors[2, i] = rotations_difference(
                problem1.camera.quaternion_rotation,
                problem2.camera.quaternion_rotation
            )
            errors[5, i] = norm(problem2.camera.position - problem1.camera.position)

            errors[3, i] = normalized_diff(errors[1, i], errors[2, i])
            errors[6, i] = normalized_diff(errors[4, i], errors[5, i])
        end
        data = hcat(data_rows, errors)
        pretty_table_withdefaults(data;
            header = header,
            highlighters = (transparent_first_col, low_value_good, high_value_bad)
        )
    end

    """
    Ensure value is always a Vector.
    If already a vector, return it.
    If scalar, wrap in a 1-element vector.
    """
    function ensure_vector(x)
        return x isa AbstractVector ? x : [x]
    end

    """
    Creates a single noise-level calibration result in the format:
    {
        "noise": 0.1,
        "method": "ours",
        "delta_f": [1.2],
        "delta_uv": [0.2, 0.3],
        "delta_skew": []
    }
    Handles scalar or vector inputs and nulls.
    """
    function create_single_noise_result(method::String, noise::Float64,
        delta_f=nothing, delta_uv=nothing, delta_skew=nothing, success_rate=nothing, delta_r=nothing, delta_t=nothing)

        result = Dict(
            "noise" => noise,
            "method" => method,
            "delta_f" => delta_f === nothing ? [] : [x === missing || x === nothing ? nothing : x for x in ensure_vector(delta_f)],
            "delta_uv" => delta_uv === nothing ? [] : [x === missing || x === nothing ? nothing : x for x in ensure_vector(delta_uv)],
            "delta_skew" => delta_skew === nothing ? [] : [x === missing || x === nothing ? nothing : x for x in ensure_vector(delta_skew)],
            "delta_r" => delta_r === nothing ? [] : [x === missing || x === nothing ? nothing : x for x in ensure_vector(delta_r)],
            "delta_t" => delta_t === nothing ? [] : [x === missing || x === nothing ? nothing : x for x in ensure_vector(delta_t)],
            "success_rate" => success_rate === nothing ? [] : success_rate,
        )

        return result
    end

    """
    Saves a list of results (vector of Dicts) to a JSON file.
    """
    function save_results_to_json(filename::String, results::Vector{Dict})
        open(filename, "w") do io
            JSON.print(io, results, 2)
        end
    end

end