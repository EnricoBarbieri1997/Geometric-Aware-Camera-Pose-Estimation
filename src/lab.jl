module Lab
    using HomotopyContinuation
    using DelimitedFiles
    using StatsBase
    using LinearAlgebra

    function shuttercameras()
        link1 = "https://gist.githubusercontent.com/PBrdng/e17d0e3bc4d983734238b9cb8386d560/raw/07272b125a6ad03c791fdf99e741318f1d85149b/3Dpoints"
        link2 = "https://gist.githubusercontent.com/PBrdng/e17d0e3bc4d983734238b9cb8386d560/raw/07272b125a6ad03c791fdf99e741318f1d85149b/2Dpoints"
        points3D = readdlm(download(link1)) |> transpose
        points2D = readdlm(download(link2)) |> transpose
        τ = 0.01
        m = 5
        @var x[1:3, 1:m], u[1:2, 1:m]
        @var R[1:3, 1:3], t[1:3], v[1:3], n[1:m]

        cams = [[R t-(n[i]*τ).*v] for i in 1:m] # m = 5 cameras

        y = [cams[i] * [x[:,i]; 1] for i in 1:m] # m = 5 images

        g = [cross(y[i],  [u[:,i]; 1]) for i in 1:m]

        k = R * R' - diagm(ones(3))
        Rconstraints = [k[i,j] for i in 1:3, j in 1:3 if i<=j]
        vconstraints = transpose(v) * v - 1

        @var l1[1:6], l2       # Lagrange multipliers 

        L = transpose(g) * g - transpose(l1) * Rconstraints - l2 * vconstraints
        Lag = differentiate(L, [vec(R); t; v; l1; l2]);

        F = System(Lag, variables = [vec(R); t; v; l1; l2], parameters = [vec(x); vec([u; n'])]);
        display(F)

        p0 = randn(ComplexF64, 30) 
        S0 = solve(F, target_parameters = p0)
        start = solutions(S0);

        N = size(points2D, 2)
        s = StatsBase.sample(collect(1:N), m; replace=false)

        X = points2D[1:3, s] 
        Y = points2D[:, s]

        p1 = [vec(X); vec(Y)]
        S1 = solve(F, start, start_parameters = p0, target_parameters = p1)
        G = System(vcat(g...), variables = [vec(R); t; v], parameters = [vec(x); vec([u; n'])]);
        function find_min(sols, p)
            sols = [r[1:15] for r in sols]
            a = map(sols) do r 
                norm(G(r, p))
            end
            i = findmin(a)
            sols[i[2]]
        end
        recovery = find_min(real_solutions(S1), p1)
        R1, t1, v1 = reshape(recovery[1:9], 3, 3), recovery[10:12], recovery[13:15]
    end

    function dinosaur()
        include(download("https://gist.githubusercontent.com/PBrdng/46436855f3755c5a959a7c5d6ba7e32b/raw/5e6be8f4c9673f0dd26010b5e0cefc8e953fc1c0/cameras.jl"))
        include(download("https://gist.githubusercontent.com/PBrdng/46436855f3755c5a959a7c5d6ba7e32b/raw/5e6be8f4c9673f0dd26010b5e0cefc8e953fc1c0/pictures.jl"))
        n₁, n₂ = 1, 2
        camera_numbers = [n₁, n₂]

        @var x[1:3]
        @var t[1:2]
        @var p[1:2,1:2]

        y = [cams[i][1:2,:]*[x; 1] for i in camera_numbers]
        z = [cams[i][3,:] ⋅ [x; 1] for i in camera_numbers] .* 648

        g = [t[i] .* y[i] - p[:,i] for i in 1:2]
        G = sum(gᵢ ⋅ gᵢ for gᵢ in g)
        # G is the function that we want to optimize

        @var λ[1:2] # the Lagrance multipliers
        L = G - sum(λ[i] * (t[i] * z[i] - 1) for i in 1:2)
        ∇L = differentiate(L, [x;t;λ])
        F = System(∇L; variables = [x;t;λ], parameters = vec(p))
        display(F)

        p₀ = randn(ComplexF64, 4)
        start = solutions(solve(F; target_parameters = p₀))

        #first, we need to preprocess the photo data
        photos = [ps[i, [2 * n₁ - 1, 2 * n₁, 2 * n₂ - 1, 2 * n₂]] for i in 1:4983]

        # the data from the dataset is incomplete.
        # the cameras did not take pictures of all world points.
        # if a world point was not captured, the entry in the data set is -1.
        filter!(pᵢ -> all(pᵢ .> 0), photos)

        # we divide the photo coordinates by 648
        # to work with coordinates between 0 and 1
        # (as explained above)
        map!(pᵢ -> pᵢ./648, photos, photos)

        function reconstruction_from_critical_points(R, pᵢ, G)
            r = real_solutions(R)
            N = map(r) do rᵢ
                G([x;t] => rᵢ[1:5], vec(p) => pᵢ)
            end
            i = findmin(N)

            return r[i[2]][1:3]
        end

        reconstructed_points = solve(
            F,
            start;
            start_parameters =  p₀,
            target_parameters = photos,
            transform_result = (R,pᵢ) -> reconstruction_from_critical_points(R,pᵢ,G)
        )
    end

    function ellipses_meet()
        @polyvar Q₁[1:2, 1:2] Q₂[1:2, 1:2] p₁[1:2] p₂[1:2]
        @polyvar x[1:2] r
        z₁ = x - p₁
        z₂ = x - p₂
        # initialize the equations for E₁ and E₂
        f₁ = (Q₁ * z₁) ⋅ (Q₁ * z₁) - r^2
        f₂ = (Q₂ * z₂) ⋅ (Q₂ * z₂) - r^2
        # initialize the equation for E₁ and E₂ being tangent
        @polyvar λ
        g = (Q₁' * Q₁) * z₁ - λ .* (Q₂' * Q₂) * z₂
        # gather everything in one system
        params = [vec(Q₁); vec(Q₂); p₁; p₂]
        F = System([f₁; f₂; g]; variables=[x; r; λ], parameters=params)
        display(F)

        q = [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, -1, 0]
        p = [vec([1 1; 1 0]); vec([0 2; 2 1]); [3, 0]; [1, 2]]
        solve(F, [[0, 0, 1, -1]], start_parameters=q, target_parameters=p)
    end

    function tritangents()
        @var h[1:3] # variables for the plane
        @var x[1:3] y[1:3] z[1:3] #variables for the contact points

        #the quadric
        Q = x[3] - x[1] * x[2]
        #the cubic C with coefficients c
        C, c = dense_poly(x, 3, coeff_name = :c)

        #generate the system P for the contact point x
        P_x = [
        h ⋅ x - 1;
        Q;
        C;
        det([h differentiate(Q, x) differentiate(C, x)])
        ]

        #generate a copy of P for the other contact points y,z
        P_y = [p([h; x; c] => [h; y; c]) for p in P_x]
        P_z = [p([h; x; c] => [h; z; c]) for p in P_x]

        #define F
        F = System([P_x; P_y; P_z]; variables = [h;x;y;z], parameters = c)
        display(F)

        c₁ = randn(ComplexF64, 20)
        #solve the system for c₁
        S = solve(F; target_parameters = c₁)

        sols = solutions(S)

        #define the coefficients for C
        c₀ = coeffs_as_dense_poly(x[1]^3+x[2]^3+x[3]^3-1, x, 3)
        #track the solutions from c₁ to c₀
        R = solve(F, sols, start_parameters = c₁, target_parameters = c₀)
        display(R)
    end

    function circles_tangent_to_conics()
        @var a[1:2] r #variables for the circle center and radius
        @var x y #variables of the circle
        @var B[1:3,1:3] #coefficients of the conics
        @var v[1:2, 1:3] #variables of the 3 points at which the circle is tangent

        circle = ([x; y] - a) ⋅ ([x; y] - a) - r
        conic  = [x; y; 1] ⋅ (B * [x; y; 1]);
        tangential_condition = det([differentiate(circle, [x, y]) differentiate(conic, [x, y])])

        conditions = [circle; conic; tangential_condition]

        #define coefficients of the three conics
        C1 = randn(3,3)
        C2 = randn(3,3)
        C3 = randn(3,3)

        #Plug in the variables of the 3 points
        #and coefficients of the 3 conics
        F = System([
            f([x; y; a; r; vec(B)] => [v[:,i]; a; r; vec(C)])
            for f in conditions
            for (i,C) in enumerate([C1, C2, C3])
            ])
        display(F)
        sol = solve(F)
        display(sol)
    end
end