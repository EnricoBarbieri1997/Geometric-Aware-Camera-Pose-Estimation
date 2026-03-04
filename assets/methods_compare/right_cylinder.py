import numpy as np
import cv2
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from calib_results import create_single_noise_result, save_results_to_json, intrinsic_difference, iterations_results_to_metrics, generate_K
import argparse
import json
from scipy.spatial.transform import Rotation

np.random.seed(42)

if __name__ == "__main__":
    np.random.seed(42)
    debug = 0
    iterations_count = 50
    results = []
    noises = np.arange(0.0, 0.04, 0.0005)
    for noise in noises:
        print(f"Noise: {noise:.2f}")
        iterations_results = []

        for i in range(iterations_count):
            cylinder = generate_cylinder(center=np.random.uniform(-100, 100, 3))            
            K = generate_K()
            cameras = [generate_camera(K) for _ in range(4)]

            projected_items_views = []
            for camera in cameras:
                projected_items = project_cylinder(cylinder)
                projected_items_views.append(projected_items)

            rows = []
            for projected_items in projected_items_views:
                q1 = projected_items.q1
                q2 = projected_items.q2
                q3 = projected_items.q3
                q4 = projected_items.q4
                c1 = projected_items.c1
                c2 = projected_items.c2
                v1 = projected_items.v1
                v2 = projected_items.v2

                # Eq. (23): v1^T * omega * v2 = 0
                # omega is parameterized with zero-skew as:
                # [[w11,    0, w13/2],
                #  [  0,  w22, w23/2],
                #  [w13/2,w23/2,w33 ]]
                v1 = np.asarray(v1).reshape(-1)
                v2 = np.asarray(v2).reshape(-1)
                if v1.size == 2:
                    v1 = np.array([v1[0], v1[1], 1.0])
                if v2.size == 2:
                    v2 = np.array([v2[0], v2[1], 1.0])

                rows.append([
                    v1[0] * v2[0],
                    0.0,
                    0.5 * (v1[0] * v2[2] + v1[2] * v2[0]),
                    v1[1] * v2[1],
                    0.5 * (v1[1] * v2[2] + v1[2] * v2[1]),
                    v1[2] * v2[2],
                ])

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

    save_results_to_json("./synthetic/right_cylinder_results.json", results)
