module IO
export read_axis_rig_lines, read_camera, read_point_cloud_yup_to_zup, read_cylinders_yup_to_zup

using ..Camera: CameraProperties

using JSON, Rotations

function _yup_to_zup(vec::Vector{<:Number})::Vector{Float64}
    x, y, z = vec
    return Float64[x, -z, y]
end

function _read_vec3(point)::Vector{Float64}
    if point isa AbstractVector
        return Float64.(point[1:3])
    end
    if point isa Dict
        if haskey(point, "x")
            return Float64[point["x"], point["y"], point["z"]]
        end
    end
    throw(ArgumentError("Point must be either [x, y, z] or {x, y, z}."))
end

function _extract_collection(json_data, keys::Vector{String})
    if json_data isa AbstractVector
        return json_data
    end
    if json_data isa Dict
        for key in keys
            if haskey(json_data, key)
                return json_data[key]
            end
        end
    end
    throw(ArgumentError("Unable to find a valid collection in json data."))
end

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

function read_point_cloud_yup_to_zup(filepath::String; object_path::String = "")::Vector{Vector{Float64}}
    json_data = read_json_object(filepath; object_path=object_path)
    points = _extract_collection(json_data, ["point_cloud", "points", "pointcloud"])
    return [_yup_to_zup(_read_vec3(point)) for point in points]
end

function read_cylinders_yup_to_zup(filepath::String; object_path::String = "")
    json_data = read_json_object(filepath; object_path=object_path)
    cylinders = _extract_collection(json_data, ["cylinders"])

    return [(
        center = _yup_to_zup(_read_vec3(get(cylinder, "center", get(cylinder, "position", nothing)))),
        axis = _yup_to_zup(_read_vec3(cylinder["axis"])),
        radius = Float64(cylinder["radius"]),
    ) for cylinder in cylinders]
end

end # module IO
