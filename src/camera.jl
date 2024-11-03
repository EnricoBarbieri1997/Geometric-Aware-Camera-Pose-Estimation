module Camera
	export CameraMatrix, build_intrinsic_matrix, build_camera_matrix

	using Rotations
	using LinearAlgebra
	using ..Space: Vec3Tuple

	function CameraMatrix(
		translation::Vec3Tuple = (0, 0, 0),
		rotation::QuatRotation{Float64} = QuatRotation(1, 0, 0, 0),
		focalLength::Number = 1,
		pixelSize::Number = 1
	)
		f = focalLength / pixelSize
		K = [f 0 0 0; 0 f 0 0; 0 0 1 0]

		P = zeros(4, 4)

		r₁ = rotation
	
		t = [translation[1], translation[2], translation[3]]
		t₁ = -r₁ * t

		P[1:3, 1:3] = r₁
		P[1:3, 4] = t₁
		P[4, 4] = 1

		return K * P
	end

	function CameraMatrix(
		translation::Vec3Tuple = (0, 0, 0),
		rotation::Vec3Tuple = (0, 0, 0),
		focalLength::Number = 1,
		pixelSize::Number = 1
	)
		f = focalLength / pixelSize
		K = [f 0 0 0; 0 f 0 0; 0 0 1 0]

		P = zeros(4, 4)

		r = RotXYZ(deg2rad.(rotation)...)
		r₁ = r' # inv(r)
	
		t = [translation[1], translation[2], translation[3]]
		t₁ = -r₁ * t

		P[1:3, 1:3] = r₁
		P[1:3, 4] = t₁
		P[4, 4] = 1

		return K * P
	end

	function build_intrinsic_matrix(focalLength::Number)
		return [focalLengthCalculated 0 0 0;
		0 focalLengthCalculated 0 0;
		0 0 focalLengthCalculated 0]
	end

	function build_camera_matrix(intrinsic, rotation, translation)
		return intrinsic * vcat(hcat(rotation, translation), [0 0 0 1])
	end
end