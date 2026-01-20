import bpy
import math
import numpy as np
import cv2 as cv
import json

camera_utils = bpy.data.texts["camera_utils"].as_module()

cameras = [
    bpy.data.objects["Camera"],
    bpy.data.objects["Camera2"]
]

X = bpy.data.objects["x"]
Y = bpy.data.objects["y"]
Z = bpy.data.objects["z"]

axis = [X, Y, Z]
scene = bpy.context.scene

scene_data = {}

for (k, camera) in enumerate(cameras):
    cylinder_contours = {}
    for (i, ax) in enumerate(axis):
        filename = bpy.path.abspath(f"//renders/{k}/{i}.png")
        src = cv.imread(filename, cv.IMREAD_GRAYSCALE)
        _, src = cv.threshold(
            src,
            123,
            255,              # value for "white"
            cv.THRESH_BINARY
        )
    #    erode_elem = cv.getStructuringElement(cv.MORPH_RECT, (2, 2))
    #    src = cv.erode(src, erode_elem)
    #    src = cv.Canny(src, 50, 200, None, 3)
        cdst = cv.cvtColor(src, cv.COLOR_GRAY2BGR)

        threshold = 50
        increase = 50
        hlines = None
        while True:
            new_threshold = threshold + increase
            current_hlines = cv.HoughLines(src, 1, np.pi / 720, new_threshold)
            if current_hlines is None or len(current_hlines) < 2:
                increase = increase * 2 / 3
                if increase < 1:
                    break
                else:
                    increase = int(increase)
                    continue
            hlines = current_hlines
            threshold = new_threshold

        print(threshold)

        lines = []
        prev_rho = 0
        prev_theta = 0

        for hline in hlines:
            rho = hline[0][0]
            theta = hline[0][1]
            print(abs(prev_rho - rho))
            print(abs(prev_theta - theta))
            if abs(prev_rho - rho) < 10 and abs(prev_theta - theta) < 1e-6:
                continue
            prev_rho = rho
            prev_theta = theta
            pdirx = math.cos(theta)
            pdiry = math.sin(theta)
            dirx = pdiry
            diry = -pdirx
            x0 = rho * pdirx
            y0 = rho * pdiry
            x1 = x0 + dirx
            y1 = y0 + diry
            line = np.array([-diry, dirx, x0*y1-x1*y0])
            
            cv.line(cdst, (int(x0+1000*dirx), int(y0+1000*diry)), (int(x0-1000*dirx), int(y0-1000*diry)), (0,0,255), 3, cv.LINE_AA)

    #        print(line.dot(np.array([x0,y0,1])))
    #        print(line.dot(np.array([x1+dirx, y1+diry, 1])))
    #        print(line.dot(np.array([x0-dirx, y0-diry, 1])))

            lines.append(line.tolist())
            
            if len(lines) > 1:
                print("end")
                break

        # cv.imshow("found", cdst)
        # cv.waitKey()
        cylinder_contours[ax.name] = lines

    scene_data[f"{k}"] = {
        "contours": cylinder_contours,
        "camera": {
            "position": list(camera.matrix_world.translation),
            "euler_rotation": [math.degrees(rad) for rad in camera.matrix_world.to_euler('XYZ')],
            "quaternion_rotation": list(camera_utils.get_R_and_T_matrix_from_blender(camera)[0].to_quaternion()),
            "intrinsics": camera_utils.scene_camera_intrinsics(camera).tolist()
        }
    }

with open(bpy.path.abspath("//scene.json"), "w") as f:
        json.dump(scene_data, f, indent=4)