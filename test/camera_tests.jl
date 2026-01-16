using CylindersBasedCameraResectioning.Camera: lookat_quaternion

using Rotations

@testset "lookat_quaternion" begin
    camera_pos = [
        24.906270720338103
        -9.562569965674594
        -4.350406017145329
    ]
    target_pos = [
        0.0
        0.0
        0.0
    ]
    quaternion = lookat_quaternion(camera_pos, target_pos; up=[0.0, 0.0, 1.0], use_model_front=true)
    q_params = Rotations.params(quaternion)
    target_q_params = [0.6279046546655354, -0.5338085283310561, -0.3668496443264594, 0.43151539739380895]

    @test isapprox(q_params, target_q_params)
    euler_angles = rad2deg.(eulerangles_from_rotationmatrix(Matrix{Float64}(quaternion)))
    target_euler_angles = [
        -80.73856258963846
        0.0
        68.99609355451608
    ]

    @test isapprox(euler_angles, target_euler_angles)
end