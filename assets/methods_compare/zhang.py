import numpy as np
import cv2
import matplotlib.pyplot as plt
from calib_results import create_single_noise_result, save_results_to_json

def plot_img_points(image_points, image_size=(640, 480)):
    img = np.zeros((image_size[1], image_size[0], 3), dtype=np.uint8)
    for p in image_points:
        cv2.circle(img, tuple(p.ravel().astype(int)), 5, (0, 255, 0), -1)
    plt.imshow(img)
    plt.title('Image Points')
    plt.show()

# 1. Generate synthetic 3D object points (checkerboard corners)
def generate_object_points(board_size=(9, 6), square_size=1.0):
    objp = np.zeros((np.prod(board_size), 3), np.float32)
    objp[:, :2] = np.indices(board_size).T.reshape(-1, 2)
    objp *= square_size
    return objp

# 2. Simulate a known camera
def simulate_camera(fx=800, fy=800, cx=320, cy=240, skew=0, noise_std=0.0, num_views=4, debug=1):
    # Intrinsics
    K = np.array([[fx, skew, cx],
                  [0, fy, cy],
                  [0,  0,  1]])
    dist_coeffs = np.zeros(5)  # No distortion
    object_points = []
    image_points = []

    objp = generate_object_points()

    for i in range(num_views):
        # Random rotation and translation
        rvec = np.random.uniform(-0.2, 0.2, 3)
        tvec = np.random.uniform(-10, 10, 3)

        # Project 3D points to image plane
        imgp, _ = cv2.projectPoints(objp, rvec, tvec, K, dist_coeffs)

        # Add noise to image points
        if noise_std > 0:
            imgp += np.random.normal(0, noise_std, imgp.shape)

        object_points.append(objp)
        image_points.append(imgp)
        if debug >= 2:
            plot_img_points(imgp)

    return object_points, image_points, K

# 3. Run Zhang calibration
def run_calibration(object_points, image_points, image_size=(640, 480)):
    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(
        object_points, image_points, image_size, None, None
    )
    return mtx, dist, ret

def normalized_diff(calculated, truth):
    if calculated == 0 and truth == 0:
        return 0.0
    denominator = truth if truth != 0 else calculated
    return abs(calculated - truth) / denominator

def intrinsic_difference(calculated, truth):
    truth = truth / truth[2, 2]
    calculated = calculated / calculated[2, 2]
    fx, fy = calculated[0, 0], calculated[1, 1]
    cx, cy = calculated[0, 2], calculated[1, 2]
    skew = calculated[0, 1]

    fx_t, fy_t = truth[0, 0], truth[1, 1]
    cx_t, cy_t = truth[0, 2], truth[1, 2]
    skew_t = truth[0, 1]

    deltaF = (normalized_diff(fx, fx_t) + normalized_diff(fy, fy_t)) / 2
    deltaUV = (normalized_diff(cx, cx_t) + normalized_diff(cy, cy_t)) / 2
    deltaSkew = 2 * abs(skew - skew_t)

    return [deltaF, deltaUV, deltaSkew]

def print_intrinsics_comparison(K_true, A_est, debug = 1):
    if debug >= 1:
      print("\nGround Truth Intrinsic Matrix K:")
      print(np.round(K_true, 2))

      print("\nEstimated Intrinsic Matrix A:")
      print(np.round(A_est, 2))

      print("\nIntrinsic Difference:")
    delta = intrinsic_difference(A_est, K_true)
    print(f"Δf: {delta[0]:.4f}, Δuv: {delta[1]:.4f}, Δskew: {delta[2]:.4f}")

debug = 1

fx = 3500 + np.random.normal(0, 50)
fy = 3500 + np.random.normal(0, 50)
cx = 0 + np.random.normal(0, 10)
cy = 0 + np.random.normal(0, 10)
skew = 0 + np.random.normal(0, 0.1)

results = []

for noise in np.arange(0, 5.5, 0.5):
    obj_pts, img_pts, gt_K = simulate_camera(
        fx=fx,
        fy=fy,
        cx=cx,
        cy=cy,
        skew=skew,
        noise_std=noise,
        num_views=4,
        debug=debug
    )
    est_K, _, error = run_calibration(obj_pts, img_pts)
    delta_f, delta_uv, delta_skew = intrinsic_difference(est_K, gt_K)
    results.append(
        create_single_noise_result(
            "zhang",
            noise / 10.0,
            delta_f,
            delta_uv,
            delta_skew,
        )
    )
    print(f"Noise STD: {noise:.1f} - Reprojection Error: {error:.4f}")
    print_intrinsics_comparison(gt_K, est_K)
    print("\n========================================\n")

save_results_to_json("zhang_results.json", results)