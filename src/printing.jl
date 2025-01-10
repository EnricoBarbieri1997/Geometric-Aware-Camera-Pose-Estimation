module Printing
    using using PrettyTables

    function print_camera_differences(original_camera, calculated_camera)
        header = (
            ["Time", "Acceleration", "Velocity", "Distance"],
            [   "s",     "[m / s²]",  "[m / s]",      "[m]"]
        )
        (["Time", "Acceleration", "Velocity", "Distance"], ["s", "[m / s²]", "[m / s]", "[m]"])

        julia> hl_p = Highlighter(
                (data, i, j) -> (j == 4) && (data[i, j] > 9),
                crayon"blue bold"
            );

        julia> hl_v = Highlighter(
                (data, i, j) -> (j == 3) && (data[i, j] > 9),
                crayon"red bold"
            );

        julia> hl_10 = Highlighter(
                (data, i, j) -> (i == 10),
                crayon"fg:white bold bg:dark_gray"
            );
        display("Original quaternion:$(round.(original_camera.quaternion_rotation, digits=2))")
        display("Calculated quaternion:$(round.(calculated_camera.quaternion_rotation, digits=2))")

        display("Actual rotation: $(original_camera.euler_rotation)")
        display("Best solution for rotation: $(calculated_camera.euler_rotation)")
        display("Difference between rotations: $(rotations_difference(calculated_camera.quaternion_rotation, original_camera.quaternion_rotation))")

        display("Actual intrinsic: $(original_camera.intrinsic)")
        display("Calculated intrinsic: $(calculated_camera.intrinsic)")

        display("Calculated translation: $(calculated_camera.position)")
        display("Actual translation: $(original_camera.position)")
        display("Difference between translations: $(translations_difference(calculated_camera.position, original_camera.position))")

        display("Camera projection matrix: $(original_camera.matrix ./ original_camera.matrix[3, 4])")
        display("Calculated projection camera matrix: $(calculated_camera.matrix ./ calculated_camera.matrix[3, 4])")
    end
end