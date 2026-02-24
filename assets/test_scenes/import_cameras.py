import bpy
from mathutils import Matrix, Vector

# --- Settings ---
RES_X = 1936
RES_Y = 1296
FOCAL_DISTANCE_LENS = 29.0 # mm

# OpenCV -> Blender camera coords conversion (x right, y down, z forward)
CV_TO_BLENDER = Matrix(((1, 0, 0),
                        (0, -1, 0),
                        (0, 0, -1)))


def apply_camera(name, K, R, t):
    cam = bpy.data.objects.get(name)
    if cam is None or cam.type != 'CAMERA':
        raise ValueError(f"Camera '{name}' not found or not a camera")

    fx = K[0][0]
    fy = K[1][1]
    cx = K[0][2]
    cy = K[1][2]
    
    SENSOR_WIDTH = FOCAL_DISTANCE_LENS * RES_X / fx

    # Render resolution and pixel aspect to match fx/fy
    scene = bpy.context.scene
    scene.render.resolution_x = RES_X
    scene.render.resolution_y = RES_Y
    scene.render.pixel_aspect_x = 1.0
    scene.render.pixel_aspect_y = fx / fy

    # Intrinsics -> Blender camera settings
    cam.data.sensor_fit = 'HORIZONTAL'
    cam.data.sensor_width = SENSOR_WIDTH
    cam.data.lens = fx * SENSOR_WIDTH / RES_X
    cam.data.shift_x = -(cx - RES_X / 2) / RES_X
    cam.data.shift_y = -(RES_Y / 2 - cy) / RES_Y

    # Extrinsics: OpenCV world->camera (R,t) -> Blender camera->world
    R_cv = Matrix(R)
    t_cv = Vector(t)
    R_bl = CV_TO_BLENDER @ R_cv
    t_bl = CV_TO_BLENDER @ t_cv

    world_to_blender_cam = R_bl.to_4x4()
    world_to_blender_cam.translation = t_bl
    cam.matrix_world = world_to_blender_cam.inverted()


# --- Camera1 ---
K1 = [
    [2393.9521661194594, 0, 932.3821770809016],
    [0.0, 2398.118540286655, 628.2649953288085],
    [0.0, 0.0, 1.0],
]
R1 = [
    [0.991948325623958, -0.12281116751614761, 0.03091822157154779],
    [0.11931476770293795, 0.8244240836579798, -0.553253031164488],
    [0.04245592420170071, 0.5524874183345906, 0.8324392753121147],
]
t1 = [0.0, 0.0, 0.0]

# --- Camera2 ---
K2 = [
    [2393.952166119456, 0, 932.3821770809021],
    [0.0, 2398.118540286656, 628.2649953288068],
    [0.0, 0.0, 1.0],
]
R2 = [
    [0.99962457261941, -0.02691925995181757, 0.0051056105519794205],
    [0.026964770035306725, 0.9995952831979801, -0.00906482130508122],
    [-0.004859525944476208, 0.009199089737387098, 0.9999458794132804],
]
t2 = [-2.1085287241533286, -13.872979092132999, 8.090162817243641]

apply_camera("Camera1", K1, R1, t1)
apply_camera("Camera2", K2, R2, t2)

print("Cameras updated.")
