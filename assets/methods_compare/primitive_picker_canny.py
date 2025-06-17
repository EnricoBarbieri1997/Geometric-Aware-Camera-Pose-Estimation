# filepath: /Users/enricobarbieri/Documents/GitHub/Geometric-Aware-Camera-Pose-Estimation/assets/methods_compare/primitive_picker_canny.py
import cv2
import numpy as np
import json
import os

# Globals
points_by_line = []
current_line = []
colors = ["red", "green", "blue", "yellow", "magenta", "cyan"]
color_bgr = [(0,0,255), (0,255,0), (255,0,0), (0,255,255), (255,0,255), (255,255,0)]
drawing = False
current_color = 0
img_original = None
img_display = None
canny_edges = None

OUTPUT_FILE = "../test_scenes/pipes/views_2.json"


def on_mouse(event, x, y, flags, param):
    global drawing, current_color, current_line, img_display, canny_edges
    shift_pressed = (flags & cv2.EVENT_FLAG_SHIFTKEY) != 0
    ctrl_pressed = (flags & cv2.EVENT_FLAG_CTRLKEY) != 0

    temp_display = img_display.copy()
    cv2.circle(temp_display, (x, y), 5, (200, 200, 200), 1)
    cv2.imshow("Image", temp_display)

    if ctrl_pressed:
        if current_line:
            dists = [np.hypot(px - x, py - y) for px, py in current_line]
            idx = int(np.argmin(dists))
            removed = current_line.pop(idx)
            print(f"CTRL: Removed point {removed} from current line.")
            img_display[:] = img_original.copy()
            color_idx = current_color
            for line in points_by_line:
                for px, py in line:
                    cv2.circle(img_display, (px, py), 3, color_bgr[color_idx % len(color_bgr)], -1)
                color_idx += 1
            for px, py in current_line:
                cv2.circle(img_display, (px, py), 3, color_bgr[current_color % len(color_bgr)], -1)
            cv2.imshow("Image", img_display)
        return

    if event == cv2.EVENT_LBUTTONDOWN or (event == cv2.EVENT_MOUSEMOVE and (drawing or shift_pressed)):
        if event == cv2.EVENT_LBUTTONDOWN:
            drawing = True
        if canny_edges is not None and 0 <= y < canny_edges.shape[0] and 0 <= x < canny_edges.shape[1]:
            for dy in range(-2, 3):
                for dx in range(-2, 3):
                    nx, ny = x + dx, y + dy
                    if 0 <= ny < canny_edges.shape[0] and 0 <= nx < canny_edges.shape[1]:
                        if canny_edges[ny, nx] != 0:
                            if (nx, ny) not in current_line:
                                current_line.append((nx, ny))
                                cv2.circle(img_display, (nx, ny), 3, color_bgr[current_color % len(color_bgr)], -1)
    elif event == cv2.EVENT_LBUTTONUP:
        drawing = False
    elif event == cv2.EVENT_RBUTTONDOWN:
        if current_line:
            points_by_line.append(current_line.copy())
            current_line.clear()
            print(f"Line {len(points_by_line)} finalized.")
            current_color += 1
        img_display[:] = img_original.copy()
        color_idx = 0
        for line in points_by_line:
            for px, py in line:
                cv2.circle(img_display, (px, py), 3, color_bgr[color_idx % len(color_bgr)], -1)
            color_idx += 1
    cv2.imshow("Image", img_display)

def on_trackbar(val):
    global img_display, canny_edges
    low = cv2.getTrackbarPos("Canny Low", "Canny")
    high = cv2.getTrackbarPos("Canny High", "Canny")
    canny_edges = cv2.Canny(img_original, low, high)
    cv2.imshow("Canny", canny_edges)

def fit_line_through_points(points):
    points = np.array(points)
    x = points[:, 0]
    y = points[:, 1]
    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    # Find endpoints for the fitted line segment (project original points onto the line)
    x_min, x_max = x.min(), x.max()
    pt1 = (float(x_min), float(m * x_min + c))
    pt2 = (float(x_max), float(m * x_max + c))
    return pt1, pt2

def save_lines_primitive_picker(points_by_line, filename):
    # Group lines: 2 red, 2 blue, 2 green
    color_order = ["red", "red", "green", "green", "blue", "blue"]
    data = {"red": [], "blue": [], "green": []}
    for i, line_points in enumerate(points_by_line):
        if i >= 6:
            break  # Only save first 6 lines
        color = color_order[i]
        pt1, pt2 = fit_line_through_points(line_points)
        data[color].append([list(map(int, pt1)), list(map(int, pt2))])
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved fitted lines to {filename}")

def main():
    global img_original, img_display, current_color
    path = input("Enter path to image: ").strip()
    img_original = cv2.imread(path)
    if img_original is None:
        print("Could not load image.")
        return
    img_display = img_original.copy()
    current_color = 0
    cv2.namedWindow("Image")
    cv2.setMouseCallback("Image", on_mouse)
    cv2.namedWindow("Canny")
    cv2.createTrackbar("Canny Low", "Canny", 50, 255, on_trackbar)
    cv2.createTrackbar("Canny High", "Canny", 150, 255, on_trackbar)
    on_trackbar(0)
    print("Left click to select points along a line.")
    print("Right click to finish the current line.")
    print("Press ESC when done.")
    while True:
        cv2.imshow("Image", img_display)
        key = cv2.waitKey(1) & 0xFF
        if key == 27:  # ESC
            if current_line:
                points_by_line.append(current_line.copy())
            break
    if len(points_by_line) < 1:
        print("No lines to save.")
        return
    save_lines_primitive_picker(points_by_line, OUTPUT_FILE)
    cv2.destroyAllWindows()

if __name__ == "__main__":
    main()
