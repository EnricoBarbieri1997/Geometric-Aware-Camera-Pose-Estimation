import numpy as np
import os
from scipy.spatial.transform import Rotation as R

# --- USER SETTINGS ---
# Set your input and output folders here
input_sparse_dir = "../assets/test_scenes/markers/colmap/"    # Folder containing cameras.txt, images.txt, points3D.txt
output_sparse_dir = "../assets/test_scenes/markers/colmap/recentered"  # New output folder

# Set your desired new origin (x, y, z)
new_origin = np.array([-1.63415, 1.73131, 2.84306])  # <--- change this to your object's coordinates

# Optional: define axis alignment using known points
# Set to None if you only want translation
x1 = new_origin  # point along desired X axis
x2 = np.array([-1.37493, 1.31547, 3.65643])  # second point along desired X axis
z1 = new_origin  # point along desired Z axis
z2 = np.array([-1.82585, 0.908253, 2.56153])  # second point along desired Z axis

# --- SCRIPT START ---

def parse_images_txt(filepath):
    images = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
        idx = 0
        while idx < len(lines):
            line = lines[idx].strip()
            if line.startswith('#') or line == '':
                idx += 1
                continue
            elems = line.split()
            image_id = int(elems[0])
            qw, qx, qy, qz = map(float, elems[1:5])
            tx, ty, tz = map(float, elems[5:8])
            camera_id = int(elems[8])
            image_name = elems[9]
            images[image_id] = {
                'qvec': np.array([qw, qx, qy, qz]),
                'tvec': np.array([tx, ty, tz]),
                'camera_id': camera_id,
                'image_name': image_name,
                'lines': [lines[idx], lines[idx+1]]  # save original lines for later
            }
            idx += 2
    return images

def parse_points3D_txt(filepath):
    points = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            elems = line.strip().split()
            point3D_id = int(elems[0])
            xyz = np.array(list(map(float, elems[1:4])))
            rest = elems[4:]  # color, error, track info
            points[point3D_id] = {'xyz': xyz, 'rest': rest}
    return points

def save_images_txt(filepath, images):
    with open(filepath, 'w') as f:
        f.write("# Image list with two lines of data per image:\n")
        f.write("#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, IMAGE_NAME\n")
        f.write("#   POINTS2D[] as (X, Y, POINT3D_ID)\n")
        for image_id in sorted(images.keys()):
            lines = images[image_id]['lines']
            elems = lines[0].strip().split()
            qw, qx, qy, qz = map(float, elems[1:5])
            tx, ty, tz = images[image_id]['tvec']  # updated
            new_line = f"{image_id} {qw} {qx} {qy} {qz} {tx} {ty} {tz} {images[image_id]['camera_id']} {images[image_id]['image_name']}\n"
            f.write(new_line)
            f.write(lines[1])  # unchanged 2D points

def save_points3D_txt(filepath, points):
    with open(filepath, 'w') as f:
        f.write("# 3D point list with one line of data per point:\n")
        f.write("#   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[] as (IMAGE_ID, POINT2D_IDX)\n")
        for point3D_id in sorted(points.keys()):
            xyz = points[point3D_id]['xyz']
            rest = points[point3D_id]['rest']
            line = f"{point3D_id} {xyz[0]} {xyz[1]} {xyz[2]} {' '.join(rest)}\n"
            f.write(line)

def build_alignment_matrix(x1, x2, z1, z2):
    x_axis = x2 - x1
    x_axis /= np.linalg.norm(x_axis)

    z_axis = z2 - z1
    z_axis /= np.linalg.norm(z_axis)

    # Orthogonalize x_axis to z_axis
    x_axis = x_axis - np.dot(x_axis, z_axis) * z_axis
    x_axis /= np.linalg.norm(x_axis)

    y_axis = np.cross(z_axis, x_axis)
    y_axis /= np.linalg.norm(y_axis)

    R_align = np.vstack([x_axis, y_axis, z_axis]).T
    return R_align

# --- MAIN PROCESS ---

os.makedirs(output_sparse_dir, exist_ok=True)

images = parse_images_txt(os.path.join(input_sparse_dir, "images.txt"))
points = parse_points3D_txt(os.path.join(input_sparse_dir, "points3D.txt"))

# Compute rotation alignment if axes are provided
if x1 is not None and x2 is not None and z1 is not None and z2 is not None:
    R_align = build_alignment_matrix(x1 - new_origin, x2 - new_origin, z1 - new_origin, z2 - new_origin)
else:
    R_align = np.eye(3)  # No rotation

# Update cameras
for image in images.values():
    qvec = image['qvec']
    tvec = image['tvec']
    C_old = tvec

    C_centered = C_old - new_origin
    C_new = R_align @ C_centered

    rotation = R.from_quat([qvec[1], qvec[2], qvec[3], qvec[0]])  # (qx, qy, qz, qw)
    R_wc = rotation.as_matrix()

    R_new = R_align @ R_wc
    t_new = -R_new @ C_new

    new_rotation = R.from_matrix(R_new)
    new_qvec = new_rotation.as_quat()
    new_qvec = [new_qvec[3], new_qvec[0], new_qvec[1], new_qvec[2]]  # (qw, qx, qy, qz)

    image['qvec'] = np.array(new_qvec)
    image['tvec'] = t_new

# Update points3D
for point in points.values():
    X_centered = point['xyz'] - new_origin
    X_new = R_align @ X_centered
    point['xyz'] = X_new

# Save outputs
save_images_txt(os.path.join(output_sparse_dir, "images.txt"), images)
os.system(f"cp {os.path.join(input_sparse_dir, 'cameras.txt')} {output_sparse_dir}/cameras.txt")
save_points3D_txt(os.path.join(output_sparse_dir, "points3D.txt"), points)

print("✅ Done! Model recentered and realigned. New sparse model saved to:", output_sparse_dir)
