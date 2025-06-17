import scipy.io
import json
import numpy as np

def mat_to_json(mat_path, json_path_view):
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, bytes):
            return obj.decode('utf-8')
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        else:
            return obj

    # Load the .mat file
    mat_data = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)

    # Remove MATLAB metadata
    mat_data = {k: v for k, v in mat_data.items() if not k.startswith("__")}
    
    cylinders = np.transpose(mat_data["CCs"], (2, 0, 1))
    singular_points = np.transpose(mat_data["VPs"], (1, 0))
    lines = np.transpose(mat_data["ll"], (1, 0,))
    # Pair lines 2 by 2
    paired_lines = [lines[i:i+2, :] for i in range(0, lines.shape[0], 2)]

    lines_object = {}
    for i, line_pair in enumerate(paired_lines):
        lines_object[f"cylinder_{i}"] = line_pair.tolist()

    view_file = [
      {
        "lines": lines_object
      }
    ]

    scene_file = {
        "configuration": 27,
        "cylinders": [
            {
              "dual_matrix": cylinders[i].tolist(),
              "singular_point": singular_points[i].tolist(),
            } for i in range(cylinders.shape[0])
        ],
        "intrinsics": mat_data["K"].tolist(),
        "cameras": [
            {
                "width": 1920,
                "height": 1080,
                "matrix": mat_data["Pgt"].tolist(),
                "image": "assets/test_scenes/roller_coaster/roller_coaster.jpg"
            }
        ]
    }

    # Convert to JSON-serializable format
    json_data_view = convert(view_file)

    # Write to JSON file
    with open(json_path_view, 'w') as f:
        json.dump(json_data_view, f, indent=2)

    json_data_scene = convert(scene_file)
    json_path_scene = json_path_view.replace("views", "scene")
    # Write to JSON file
    with open(json_path_scene, 'w') as f:
        json.dump(json_data_scene, f, indent=2)

# Example usage:
mat_to_json("/Users/enricobarbieri/Documents/GitHub/minimal_solvers_conics/images/roller_coaster_data.mat", "../test_scenes/roller_coaster/views.json")
