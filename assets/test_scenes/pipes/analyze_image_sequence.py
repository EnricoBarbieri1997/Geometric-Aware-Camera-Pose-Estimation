import os
import cv2
import numpy as np
import json
from sklearn.linear_model import RANSACRegressor

def extract_lines_ransac(image, color):
    # Extract pixels matching the color
    if color == 'red':
        mask = cv2.inRange(image, (0, 0, 127), (80, 80, 255))
    elif color == 'green':
        mask = cv2.inRange(image, (0, 127, 0), (80, 255, 80))
    elif color == 'blue':
        mask = cv2.inRange(image, (127, 0, 0), (255, 80, 80))
    else:
        return []
    
    # Find coordinates of non-zero pixels
    ys, xs = np.nonzero(mask)
    if len(xs) < 2:
        return []

    points = np.column_stack((xs, ys))
    lines = []
    
    # We will run RANSAC multiple times to detect multiple lines
    for _ in range(2):
        if len(points) < 2:
            break
        try:
            model = RANSACRegressor(residual_threshold=2, min_samples=2)
            if color == 'blue':
                model.fit(points[:, 1].reshape(-1, 1), points[:, 0])
            else:
                model.fit(points[:, 0].reshape(-1, 1), points[:, 1])
            
            # Generate two points on the detected line for the segment
            if color != 'blue':
                x_min, x_max = np.min(points[:, 0]), np.max(points[:, 0])
                y_min = model.predict([[x_min]])[0]
                y_max = model.predict([[x_max]])[0]
            else:
                y_min, y_max = np.min(points[:, 1]), np.max(points[:, 1])
                x_min = model.predict([[y_min]])[0]
                x_max = model.predict([[y_max]])[0]
            
            lines.append([[int(x_min), int(y_min)], [int(x_max), int(y_max)]])
            
            # Remove inliers to allow detection of the next line
            inlier_mask = model.inlier_mask_
            points = points[~inlier_mask]
        except:
            break

    return lines

def process_all_frames():
    output = []
    
    # Assume all folders have same number of frames with same names
    folder_paths = {'red': './animation/lines/red', 'green': './animation/lines/green', 'blue': './animation/lines/blue'}
    filenames = sorted(os.listdir(folder_paths['red']))

    for filename in filenames:
        frame_result = {"lines": {"red": [], "green": [], "blue": []}}
        for color, path in folder_paths.items():
            image_path = os.path.join(path, filename)
            image = cv2.imread(image_path)
            lines = extract_lines_ransac(image, color)
            frame_result["lines"][color] = lines
        output.append(frame_result)
    
    with open('./animation/views.json', 'w') as f:
        json.dump(output, f, indent=2)

if __name__ == "__main__":
    process_all_frames()
