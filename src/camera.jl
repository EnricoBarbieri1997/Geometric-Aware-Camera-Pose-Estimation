module Camera
	export CameraProperties, build_camera_matrix, build_intrinsic_matrix, build_camera_matrix, lookat_rotation

	using Rotations
	using LinearAlgebra
	using ..Space: Vec3Tuple

	@kwdef mutable struct CameraProperties
		position::Vector{Number} = [0, 0, 0]
		euler_rotation::Vector{Number} = [0, 0, 0]
		quaternion_rotation::QuatRotation{Float64} = one(QuatRotation)
		focal_length::Number = 1
		matrix::Matrix{<:Number} = Matrix(I, 3, 4)
	end

	function build_camera_matrix(
		translation::Union{Array{<:Number}, Vector{<:Number}} = [0, 0, 0],
		rotation::Matrix{<:Number} = Matrix(I, 3, 3),
		focalLength::Number = 1,
		pixelSize::Number = 1
	)
		f = focalLength / pixelSize
		K = [f 0 0 0; 0 f 0 0; 0 0 1 0]

		M = zeros(4, 4)

		r₁ = rotation
	
		t = [translation[1], translation[2], translation[3]]
		t₁ = -r₁ * t

		M[1:3, 1:3] = r₁
		M[1:3, 4] = t₁
		M[4, 4] = 1

		return K * M
	end

	function build_camera_matrix(
		translation::Union{Array{<:Number}, Vector{<:Number}} = [0, 0, 0],
		rotation::QuatRotation{Float64} = QuatRotation(1, 0, 0, 0),
		focal_length::Number = 1,
		pixel_size::Number = 1
	)
		K = build_intrinsic_matrix(focal_length, pixel_size)

		M = zeros(4, 4)

		r₁ = rotation' # inv(r)
	
		t₁ = -r₁ * translation

		M[1:3, 1:3] = r₁
		M[1:3, 4] = t₁
		M[4, 4] = 1

		return K * M
	end

	function build_camera_matrix(
		translation::Union{Array{<:Number}, Vector{<:Number}} = [0, 0, 0],
		rotation::Union{Array{<:Number}, Vector{<:Number}} = [0, 0, 0],
		focal_length::Number = 1,
		pixel_size::Number = 1
	)
		K = build_intrinsic_matrix(focal_length, pixel_size)

		M = zeros(4, 4)

		r = RotXYZ(deg2rad.(rotation)...)
		r₁ = r' # inv(r)
	
		t₁ = -r₁ * translation

		M[1:3, 1:3] = r₁
		M[1:3, 4] = t₁
		M[4, 4] = 1

		return K * M
	end

	function build_intrinsic_matrix(focal_length::Number, pixel_size::Number = 1)
		f = focal_length / pixel_size
		return [
			f 0 0 0;
			0 f 0 0;
			0 0 1 0
		]
	end

	function build_camera_matrix(intrinsic, rotation, translation)
		r₁ = rotation' # inv(r)
		t₁ = -r₁ * translation
		return intrinsic * vcat(hcat(r₁, t₁), [0 0 0 1])
	end

	function lookat_axis(eye, at, up)
		zaxis = normalize(at - eye)
		xaxis = normalize(cross(up, zaxis))
		yaxis = cross(zaxis, xaxis)

		return xaxis, yaxis, zaxis
	end

	function lookat_rotation(eye, at, up)
		xaxis, yaxis, zaxis = lookat_axis(eye, at, up)

		return [xaxis[1] yaxis[1] zaxis[1];
			xaxis[2] yaxis[2] zaxis[2];
			xaxis[3] yaxis[3] zaxis[3]] * RotZ(π)
	end

	function lookat_matrix(eye, at, up)
		xaxis, yaxis, zaxis = lookat_axis(eye, at, up)

		return [xaxis[1] yaxis[1] zaxis[1] 0;
			xaxis[2] yaxis[2] zaxis[2] 0;
			xaxis[3] yaxis[3] zaxis[3] 0;
			-dot(xaxis, eye) -dot(yaxis, eye) -dot(zaxis, eye) 1]
	end
end