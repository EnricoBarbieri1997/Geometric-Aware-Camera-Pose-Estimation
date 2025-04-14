import numpy as np
import cv2
import matplotlib.pyplot as plt
from calib_results import create_single_noise_result, save_results_to_json, intrinsic_difference, iterations_results_to_metrics, generate_K

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
def simulate_camera(noise_std=0.0, num_views=4, debug=1):
    # Intrinsics
    K = generate_K()
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
cx = 1750 + np.random.normal(0, 10)
cy = 1750 + np.random.normal(0, 10)
skew = 0 + np.random.normal(0, 0.1)

results = []

zhang_variants = [
    "zhang_4",
    "zhang_30",
]

views_for_variant = {
    "zhang_4": 4,
    "zhang_30": 30,
}

for noise in np.arange(0, 5.5, 0.5):
    for variant in zhang_variants:
        iteration_results = []
        for i in range(50):
            obj_pts, img_pts, gt_K = simulate_camera(
                noise_std=noise / 10.0,
                num_views=views_for_variant[variant],
                debug=debug
            )
            est_K, _, error = run_calibration(obj_pts, img_pts)
            iteration_results.append(intrinsic_difference(est_K, gt_K))
            print(f"Noise STD: {noise:.1f} - Reprojection Error: {error:.4f}")
            print_intrinsics_comparison(gt_K, est_K)
            print("\n========================================\n")
        

        results.append(
            create_single_noise_result(
                variant,
                noise / 10.0,
                *iterations_results_to_metrics(iteration_results),
            )
        )

save_results_to_json("zhang_results.json", results)