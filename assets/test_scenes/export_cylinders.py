import bpy
import json
import os
from mathutils import Vector


def cylinder_to_dict(obj):
    center = obj.matrix_world @ Vector((0.0, 0.0, 0.0))
    axis = obj.matrix_world.to_3x3() @ Vector((0.0, 0.0, 1.0))
    if axis.length != 0.0:
        axis.normalize()
    radius = (obj.dimensions.x + obj.dimensions.y) * 0.25
    return {
        "center": [center.x, center.y, center.z],
        "axis": [axis.x, axis.y, axis.z],
        "radius": radius,
    }


def main():
    cylinders = [
        obj for obj in bpy.context.scene.objects
        if obj.type == "MESH" and obj.name.startswith("Cylinder")
    ]

    if bpy.data.filepath:
        output_dir = os.path.dirname(bpy.data.filepath)
    else:
        output_dir = bpy.path.abspath("//")

    output_path = os.path.join(output_dir, "cylinders.json")
    payload = {"cylinders": [cylinder_to_dict(obj) for obj in cylinders]}

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=4)

    print(f"Exported {len(cylinders)} cylinders to {output_path}")


if __name__ == "__main__":
    main()
