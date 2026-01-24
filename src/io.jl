module IO
export read_axis_rig_lines, read_camera

using ..Camera: CameraProperties

using JSON, Rotations

"Reads a json and extract the specified object path."
function read_json_object(filepath::String; object_path::String = "")
    file = open(filepath, "r")
    json_data = JSON.parse(file)
    close(file)

    if object_path != ""
        for key in split(object_path, ".")
            json_data = json_data[key]
        end
    end

    return json_data
end

" Reads camera properties from a file. "
function read_camera(filepath::String; object_path::String = "")::CameraProperties
    json_data = read_json_object(filepath; object_path=object_path)

    camera = CameraProperties()
    camera.position = json_data["position"]
    camera.euler_rotation = json_data["euler_rotation"]
    camera.quaternion_rotation = QuatRotation(0.311379, 0.222926, -0.527197, 0.758558) * QuatRotation(json_data["quaternion_rotation"])
    camera.intrinsic = Float64.((hcat(json_data["intrinsics"]...)'))

    return camera
end

function read_axis_rig_lines(filepath::String; object_path::String = "")::Vector{Vector{Vector{Float64}}}
    json_data = read_json_object(filepath; object_path=object_path)

    lines::Vector{Vector{Vector{Float64}}} = []
    push!(lines, json_data["x"])
    push!(lines, json_data["y"])
    push!(lines, json_data["z"])

    return lines
end

end # module IO
