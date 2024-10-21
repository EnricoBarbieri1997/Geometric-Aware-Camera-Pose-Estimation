module Debug
    markers = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        # [-1, 0, 0],
        # [0, -1, 0],
        # [0, 0, -1]
    ]

    markers = [[(marker * 2)..., 1] + [0, 0, 10, 0] for marker in markers]

    function plot3DMarkers()
        for (i, marker) in enumerate(markers)
            scatter!(ax3, (marker[1], marker[2], marker[3]), color=colors[i])
        end
    end

    function plot2DMarkers()
#     for (i, marker) in enumerate(markers)
#         projectedMarker = cameraProjectionMatrix * marker
#         projectedMarker = projectedMarker ./ projectedMarker[3]
#         scatter!(ax2,(projectedMarker[1], -projectedMarker[2]), color=colors[i])
#     end
    end
end