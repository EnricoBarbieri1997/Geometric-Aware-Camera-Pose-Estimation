module Printing
    export print_camera_differences, print_error_analysis

    using ..Utils: vector_difference, matrix_difference, rotations_difference, translations_difference, intrinsic_difference
    using PrettyTables
    using Rotations: params as rotations_params

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
        pretty_table_withdefaults(
            originalcalculateddata(
                original_camera.intrinsic[intrinsic_indexes],
                calculated_camera.intrinsic[intrinsic_indexes]
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
            );
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
        if (length(errors) > 4)
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
end