import numpy as np
from scipy.spatial.transform import Rotation as R
import os

# === CONFIGURATION ===
TRANSFORM_PATH = "./assets/test_scenes/markers/colmap/transform.txt"
INPUT_TXT_DIR = "./assets/test_scenes/markers/colmap/"
OUTPUT_TXT_DIR = "./assets/test_scenes/markers/colmap/recentered/"

# Create output directory
os.makedirs(OUTPUT_TXT_DIR, exist_ok=True)

# === LOAD TRANSFORMATION MATRIX ===
T = np.loadtxt(TRANSFORM_PATH)  # 4x4 transformation matrix

# === FUNCTION: TRANSFORM A POINT ===
def transform_point(p):
    p_h = np.append(p, 1)  # make homogeneous
    return (T @ p_h)[:3]

# === PROCESS POINTS3D ===
with open(os.path.join(INPUT_TXT_DIR, "points3D.txt"), "r") as fin, \
     open(os.path.join(OUTPUT_TXT_DIR, "points3D.txt"), "w") as fout:
    for line in fin:
        if line.startswith("#") or line.strip() == "":
            fout.write(line)
            continue
        elems = line.strip().split()
        xyz = np.array(list(map(float, elems[1:4])))
        new_xyz = transform_point(xyz)
        fout.write(f"{elems[0]} {new_xyz[0]} {new_xyz[1]} {new_xyz[2]} " + " ".join(elems[4:]) + "\n")

# === FUNCTION: Convert qvec and tvec to 4x4 pose matrix ===
def qvec_tvec_to_pose(qvec, tvec):
    R_mat = R.from_quat([qvec[1], qvec[2], qvec[3], qvec[0]]).as_matrix()
    pose = np.eye(4)
    pose[:3, :3] = R_mat
    pose[:3, 3] = tvec
    return pose

# === FUNCTION: Convert pose matrix back to qvec and tvec ===
def pose_to_qvec_tvec(pose):
    R_mat = pose[:3, :3]
    tvec = pose[:3, 3]
    q = R.from_matrix(R_mat).as_quat()
    # COLMAP uses (qw, qx, qy, qz)
    qvec = [q[3], q[0], q[1], q[2]]
    return qvec, tvec

R_fix = np.array([
    [1, 0,  0],
    [0, 0, -1],
    [0, 1,  0]
])

def fix_camera_pose_after_transform(pose):
    print("fixing")
    fixed = np.eye(4)
    fixed[:3, :3] = R_fix @ pose[:3, :3]
    fixed[:3, 3] = R_fix @ pose[:3, 3]
    return fixed

def quat_to_matrix(qvec):
    return R.from_quat([qvec[1], qvec[2], qvec[3], qvec[0]]).as_matrix()

def matrix_to_quat(R_mat):
    q = R.from_matrix(R_mat).as_quat()
    return [q[3], q[0], q[1], q[2]]  # COLMAP: qw, qx, qy, qz

def transform_camera_pose(qvec, tvec):
    R_wc = quat_to_matrix(qvec)
    C = -R_wc.T @ tvec  # camera center in world coordinates

    # Transform camera center with Blender transform
    C_new = transform_point(C)

    # Rotate camera using T and axis fix
    R_cw_new = T[:3, :3] @ R_wc.T
    R_wc_new = R_cw_new.T
    tvec_new = -R_wc_new @ C_new
    qvec_new = matrix_to_quat(R_wc_new)
    return qvec_new, tvec_new

# === PROCESS IMAGES.TXT ===
with open(os.path.join(INPUT_TXT_DIR, "images.txt"), "r") as fin, \
     open(os.path.join(OUTPUT_TXT_DIR, "images.txt"), "w") as fout:
    for line in fin:
        if line.startswith("#") or line.strip() == "":
            fout.write(line)
            continue
        elems = line.strip().split()
        image_id = elems[0]
        qvec = list(map(float, elems[1:5]))
        tvec = np.array(list(map(float, elems[5:8])))
        qvec_new, tvec_new = transform_camera_pose(qvec, tvec)
        cam_id = elems[8]
        name = elems[9]
        new_line = f"{image_id} {' '.join(map(str, qvec_new))} {' '.join(map(str, tvec_new))} {cam_id} {name}\n"
        fout.write(new_line)

        # Second line (2D points) is copied as-is
        fout.write(fin.readline())
