import bpy
import numpy as np

def scene_camera_intrinsics():

    scene = bpy.context.scene
    render = scene.render
    cam = bpy.data.objects["Camera"]
    camd = cam.data  # bpy.types.Camera

    # Resolution in pixels (with percentage scale)
    scale = render.resolution_percentage / 100.0
    res_x = render.resolution_x * scale
    res_y = render.resolution_y * scale

    # Pixel aspect ratio
    par_x = render.pixel_aspect_x
    par_y = render.pixel_aspect_y

    # Sensor fit determines which dimension drives the FOV
    # Blender uses: HORIZONTAL / VERTICAL / AUTO
    sensor_fit = camd.sensor_fit
    sensor_w = camd.sensor_width
    sensor_h = camd.sensor_height

    if sensor_fit == 'AUTO':
        sensor_fit = 'HORIZONTAL' if res_x >= res_y else 'VERTICAL'

    # Focal length in mm
    f_mm = camd.lens

    # Convert focal length to pixels (accounting for pixel aspect)
    if sensor_fit == 'HORIZONTAL':
        # fx from horizontal sensor size
        fx = (res_x * f_mm) / sensor_w
        fy = fx * (par_x / par_y)
    else:
        # fy from vertical sensor size
        fy = (res_y * f_mm) / sensor_h
        fx = fy * (par_y / par_x)

    # Principal point: center + shift (shift is in sensor units, normalized)
    # In Blender, shift_x, shift_y are relative to the sensor, not pixels.
    cx = res_x * (0.5 - camd.shift_x)
    cy = res_y * (0.5 + camd.shift_y)

    K = np.matrix((
        (fx, 0.0, cx),
        (0.0, fy,  cy),
        (0.0, 0.0, 1.0),
    ))
    
    return K