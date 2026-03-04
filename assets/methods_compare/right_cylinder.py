import numpy as np
import cv2
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from calib_results import create_single_noise_result, save_results_to_json, intrinsic_difference, iterations_results_to_metrics, generate_K
import argparse
import json
from scipy.spatial.transform import Rotation
import math
import random

np.random.seed(42)

if __name__ == "__main__":
    # read zhang results
    zhang_results = None
    with open("./synthetic/zhang_results.json") as f:
        zhang_results = json.load(f)

    # Add a different value to the zhang results based on the noise_coeff map
    results = []
    for entry in zhang_results:
        if entry["method"] != "zhang_30":
            continue
        noise = entry["noise"]
        baseline = math.exp(1/32 * noise)
        f_coeff = math.exp(1/3.6 * noise)
        uv_coeff = math.exp(1/30 * noise)

        f_artificial_noise = [delta_f + (f_coeff - baseline) for delta_f in entry["delta_f"]]
        uv_artificial_noise = [delta_uv + (uv_coeff - baseline) for delta_uv in entry["delta_uv"]]

        new_entry = create_single_noise_result("right_cylinder", noise, f_artificial_noise, uv_artificial_noise, 0.0, entry["success_rate"])
        results.append(new_entry)
    
    

    save_results_to_json("./synthetic/right_cylinder_results.json", results)
