import numpy as np
import cv2
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from calib_results import create_single_noise_result, save_results_to_json, intrinsic_difference, iterations_results_to_metrics, generate_K
import argparse
import json
from scipy.spatial.transform import Rotation

np.random.seed(42)

def generate_camera(K):
    # Camera position (random, but looking at world origin)
    camera_pos = np.random.uniform(-500, 500, 3)
    look_at = np.array([0, 0, 0])
    up_vector = np.array([0, 1, 0])

    def look_at_matrix(pos, target, up):
        z = (pos - target)
        z /= np.linalg.norm(z)
        x = np.cross(up, z)
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        R = np.vstack([x, y, z])
        return R

    R = look_at_matrix(camera_pos, look_at, up_vector)
    t = -R @ camera_pos.reshape(3, 1)
    Rt = np.hstack([R, t])
    P = K @ Rt
    return P, R, t

def generate_cylinder(radius=100, length=300, center=np.zeros(3)):
    n = 100
    theta = np.linspace(0, 2*np.pi, n)
    circle = np.array([radius * np.cos(theta), radius * np.sin(theta), np.zeros(n)])

    # Tangent reference points at π/4 and 3π/4 and opposite
    q1_idx = int(n * 1/8)
    q2_idx = int(n * 3/8)
    q3_idx = int(n * 5/8)
    q4_idx = int(n * 7/8)

    base1 = circle + center.reshape(3,1)
    base2 = circle + center.reshape(3,1) + np.array([[0], [0], [length]])
    
    # Tangent 3D points Q1 to Q4 on base edges
    Q1 = base1[:, q1_idx]
    Q2 = base1[:, q2_idx]
    Q3 = base2[:, q3_idx]
    Q4 = base2[:, q4_idx]

    return base1, base2, Q1, Q2, Q3, Q4

def project_points(P, points3D):
    if points3D.shape[0] == 3:
        points3D = np.vstack([points3D, np.ones((1, points3D.shape[1]))])
    points2D = P @ points3D
    return points2D[:2] / points2D[2]

def plot_projection(img_points, base1, base2):
    plt.figure(figsize=(8, 8))
    plt.plot(img_points[0], img_points[1], 'b.')
    n = base1.shape[1]
    for i in range(0, n, n // 20):  # Draw fewer lines for clarity
      pt1 = project_points(P, base1[:, i:i+1])
      pt2 = project_points(P, base2[:, i:i+1])
      plt.plot([pt1[0, 0], pt2[0, 0]], [pt1[1, 0], pt2[1, 0]], 'k-', linewidth=0.5)
    # plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title("Projected Cylinder")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.grid(True)
    plt.show()

def plot_3d_scene(base1, base2, cameras):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Cylinder
    ax.plot(base1[0], base1[1], base1[2], 'b')
    ax.plot(base2[0], base2[1], base2[2], 'b')
    for i in range(base1.shape[1]):
        ax.plot([base1[0, i], base2[0, i]],
                [base1[1, i], base2[1, i]],
                [base1[2, i], base2[2, i]], 'c', linewidth=0.5)

    for camera in cameras:
      _, R, t = camera
      camera_center = -R.T @ t.flatten()
      # Camera center
      ax.scatter(*camera_center, color='red', s=50)
      
      # Camera orientation axes
      scale = 100
      for i, color in zip(range(3), ['r', 'g', 'b']):
          ax.quiver(*camera_center,
                    *R[i]*scale, color=color, linewidth=2)

    ax.set_title("3D Scene with Cylinder and Camera")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.grid(True)
    ax.set_box_aspect([1,1,1])
    ax.set_aspect('equal')
    plt.show()

def fit_line(p1, p2):
    """Fit a homogeneous line through two 2D points."""
    p1_h = np.array([p1[0], p1[1], 1.0])
    p2_h = np.array([p2[0], p2[1], 1.0])
    return np.cross(p1_h, p2_h)

def line_intersection(l1, l2):
    """Intersect two lines in homogeneous coordinates."""
    pt_h = np.cross(l1, l2)
    return pt_h / pt_h[2]

def fit_ellipse(points2D):
    """Fit a conic to 2D points using OpenCV."""
    pts = points2D.T.astype(np.float32)
    ellipse = cv2.fitEllipse(pts)
    return ellipse

def print_intrinsics_comparison(K_true, A_est, debug = 1):
    K_true = K_true / K_true[2, 2]
    A_est = A_est / A_est[2, 2]

    if debug >= 1:
      print("\nGround Truth Intrinsic Matrix K:")
      print(np.round(K_true, 2))

      print("\nEstimated Intrinsic Matrix A:")
      print(np.round(A_est, 2))

      print("\nIntrinsic Difference:")
    delta = intrinsic_difference(A_est, K_true)
    print(f"Δf: {delta[0]:.4f}, Δuv: {delta[1]:.4f}, Δskew: {delta[2]:.4f}")

def print_intrinsics_comparison(K_true, A_est, debug = 1):
    if debug >= 1:
      print("\nGround Truth Intrinsic Matrix K:")
      print(np.round(K_true, 2))

      print("\nEstimated Intrinsic Matrix A:")
      print(np.round(A_est, 2))

      print("\nIntrinsic Difference:")
    delta = intrinsic_difference(A_est, K_true)
    print(f"Δf: {delta[0]:.4f}, Δuv: {delta[1]:.4f}, Δskew: {delta[2]:.4f}")

def estimate_extrinsics(A, v1_h, v2_h, p1_2D, p2_2D):
    A_inv = np.linalg.inv(A)

    # Step 1–2: Estimate r3 and r2
    r3 = A_inv @ v1_h
    r2 = A_inv @ v2_h

    r3 = r3 / np.linalg.norm(r3)
    r2 = r2 / np.linalg.norm(r2)

    # Step 3: r1 = r2 × r3
    r1 = np.cross(r2, r3)
    r1 = r1 / np.linalg.norm(r1)

    # Re-orthogonalize r2 using cross product to ensure right-handed frame
    r2 = np.cross(r3, r1)

    # Form full rotation matrix
    R_est = np.vstack([r1, r2, r3]).T  # Columns = r1 r2 r3

    # Step 4: Estimate translation vector from projected ellipse centers
    # p1_2D and p2_2D are 2D projections of P1 and P2
    p1_h = np.append(p1_2D, 1)
    p2_h = np.append(p2_2D, 1)

    P1_est = A_inv @ p1_h
    P2_est = A_inv @ p2_h

    # Now compute translation (assuming world origin at base P1)
    # The paper uses direction of cylinder axis from vanishing point
    # We compute t as position of P1 in camera coordinates
    t_est = P1_est / np.linalg.norm(P1_est)  # scale is up to choice

    return R_est, t_est

def normalized_diff(calculated, truth):
    if calculated == 0 and truth == 0:
        return 0.0
    denominator = truth if truth != 0 else calculated
    return abs(calculated - truth) / denominator

def intrinsic_difference(calculated, truth):
    calculated_est = calculated / calculated[2, 2]
    truth_est = truth / truth[2, 2]
    fx, fy = calculated_est[0, 0], calculated_est[1, 1]
    cx, cy = calculated_est[0, 2], calculated_est[1, 2]
    skew = calculated_est[0, 1]

    fx_t, fy_t = truth_est[0, 0], truth_est[1, 1]
    cx_t, cy_t = truth_est[0, 2], truth_est[1, 2]
    skew_t = truth_est[0, 1]

    deltaF = (normalized_diff(fx, fx_t) + normalized_diff(fy, fy_t)) / 2
    deltaUV = (normalized_diff(cx, cx_t) + normalized_diff(cy, cy_t)) / 2
    deltaSkew = abs(skew - skew_t)

    return [deltaF, deltaUV, deltaSkew]

def homogeneous_line_intercept(x, line):
    # Assuming line is [a, b, c] representing ax + by + c = 0
    a, b, c = line
    if b == 0:
        return 0  # Avoid division by zero; you can handle this differently
    return -(a * x + c) / b

def normalize(v):
    return v / np.linalg.norm(v)

def homogeneous_line_from_points(p1, p2):
    return np.cross(p1, p2)

def add_noise_to_line(line, noise):
    x1 = 0
    p1 = np.array([x1, homogeneous_line_intercept(x1, line), 1.0]) + normalize(np.append(np.random.rand(2) - 0.5, 0)) * noise

    x2 = 1000
    p2 = np.array([x2, homogeneous_line_intercept(x2, line), 1.0]) + normalize(np.append(np.random.rand(2) - 0.5, 0)) * noise

    return homogeneous_line_from_points(p1, p2)

def add_noise_to_lines(line1, line2, noise):
    noisy_line1 = add_noise_to_line(line1, noise)
    noisy_line2 = add_noise_to_line(line2, noise)
    return noisy_line1, noisy_line2

def load_cylinder_data(path):
    with open(path, "r") as f:
        data = json.load(f)["cylinders"][0]

        center = np.array(data["position"])
        radius = data["radius"]
        length = data["length"]

        return generate_cylinder(radius=radius, length=length, center=center)

def load_camera_data(path):
    with open(path, "r") as f:
        data = json.load(f)
        K = data["intrinsics"]
        cameras = []
        for cam in data.get("cameras", []):
            R = Rotation.from_quat(np.array(cam["R"])).as_matrix()
            t = np.array(cam["t"]).reshape(3, 1)
            P = np.multiply(np.hstack((K, np.array([0, 0, 1]).reshape(3, 1))), np.hstack((R, t)))
            cameras.append((P, R, t))

        return K, cameras

def load_camera_ellipse_projection(path):
    views = []
    with open(path, "r") as f:
        data = json.load(f)
        for view in data:
          centers = []
          Qs = []
          for ellipse in view["ellipses"]:
              (cx, cy), (w, h), angle = cv2.fitEllipse(np.array(ellipse, dtype=np.int32))

              # Convert angle and direction into a vector
              radians = np.deg2rad(angle)
              dx = np.sin(radians)
              dy = -np.cos(radians)

              # Vertical direction vector (along height)
              vertical = np.array([dx, dy])

              # Half the vertical height
              offset = (h / 2.0) * vertical

              center = np.array([cx, cy])
              top = center - offset
              bottom = center + offset

              centers.append(center)
              Qs.append(top)
              Qs.append(bottom)

          views.append((centers, Qs))

    return views

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cylinder Calibration Experiment")
    parser.add_argument("--views", type=str, help="Path to saved JSON data instead of generating")
    parser.add_argument("--scene", type=str, help="Path to saved JSON data instead of generating")
    args = parser.parse_args()

    use_saved = args.views != None and args.scene != None

    np.random.seed(42)
    debug = 0
    iterations_count = 50 if not use_saved else 1
    results = []
    noises = np.arange(0.0, 0.04, 0.0005) if not use_saved else [0.0]
    print(use_saved)
    for noise in noises:
        print(f"Noise: {noise:.2f}")
        iterations_results = []

        for i in range(iterations_count):
            if use_saved:
                base1, base2, Q1, Q2, Q3, Q4 = load_cylinder_data(args.scene)
                K, cameras = load_camera_data(args.scene)
                ellipses_per_view = load_camera_ellipse_projection(args.views)
            else:
                base1, base2, Q1, Q2, Q3, Q4 = generate_cylinder(center=np.random.uniform(-100, 100, 3))
                K = generate_K()
                cameras = [generate_camera(K) for _ in range(4)]
            
            P1 = base1[:, 0]
            P2 = base2[:, 0]
            Qs = [Q1, Q2, Q3, Q4]
            views = []

            for camera_index, camera in enumerate(cameras):
                P, R, t = camera

                if use_saved:
                    qs_2D = ellipses_per_view[camera_index][1]
                else:
                    qs_2D = [project_points(P, q.reshape(3, 1))[:, 0] for q in Qs]

                l1 = fit_line(qs_2D[0], qs_2D[2])
                l2 = fit_line(qs_2D[1], qs_2D[3])

                if noise > 0:
                    l1, l2 = add_noise_to_lines(l1, l2, noise)

                v1 = line_intersection(l1, l2)

                q1q2 = fit_line(qs_2D[0], qs_2D[1])
                q3q4 = fit_line(qs_2D[2], qs_2D[3])
                v2 = line_intersection(q1q2, q3q4)

                if use_saved:
                    p1_2D = ellipses_per_view[camera_index][0][0]
                    p2_2D = ellipses_per_view[camera_index][0][1]
                else:
                    p1_2D = project_points(P, P1.reshape(3,1))[:, 0]
                    p2_2D = project_points(P, P2.reshape(3,1))[:, 0]

                views.append((v1, v2, p1_2D, p2_2D))

            rows = []
            for (v1, v2, _, _) in views:
                x1, y1, w1 = v1
                x2, y2, w2 = v2
                row = [
                    x1 * x2,
                    0,
                    x1 * w2 + w1 * x2,
                    y1 * y2,
                    y1 * w2 + w1 * y2,
                    w1 * w2
                ]
                rows.append(row)

            A = np.array(rows)
            _, _, Vt = np.linalg.svd(A)
            omega_vec = Vt[-1]

            omega = np.array([
                [omega_vec[0], 0, omega_vec[2]/2],
                [0, omega_vec[3], omega_vec[4]/2],
                [omega_vec[2]/2, omega_vec[4]/2, omega_vec[5]]
            ])

            try:
                A_inv = np.linalg.cholesky(omega).T
                A_est = np.linalg.inv(A_inv)

                iterations_results.append(intrinsic_difference(A_est, K) + [1/iterations_count])
                print_intrinsics_comparison(K, A_est, debug)

            except np.linalg.LinAlgError as e:
                print("⚠️ Cholesky decomposition failed:", e)

        results.append(create_single_noise_result("right_cylinder", noise, *iterations_results_to_metrics(iterations_results)))
        print("\n" + "="*50 + "\n")

    save_results_to_json("right_cylinder_results.json", results)