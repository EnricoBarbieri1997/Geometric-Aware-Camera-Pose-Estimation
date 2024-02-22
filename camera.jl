module Camera
	export CameraMatrix

	using Rotations
	using LinearAlgebra
	using ..Space: Vec3Tuple

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
end