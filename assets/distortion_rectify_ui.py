import cv2
import numpy as np
from scipy.optimize import minimize
import os

# Globals
points_by_line = []
current_line = []
colors = [(255, 0, 0), (0, 255, 0), (0, 0, 255),
          (255, 255, 0), (255, 0, 255), (0, 255, 255)]
drawing = False
current_color = 0
img_original = None
img_display = None
canny_edges = None  # Store the current Canny edge map
LINES_FILE = "lines.npy"

def on_mouse(event, x, y, flags, param):
    global drawing, current_color, current_line, img_display, canny_edges
    shift_pressed = (flags & cv2.EVENT_FLAG_SHIFTKEY) != 0
    ctrl_pressed = (flags & cv2.EVENT_FLAG_CTRLKEY) != 0

    # Always show a small circle at the cursor position
    temp_display = img_display.copy()
    cv2.circle(temp_display, (x, y), 5, (200, 200, 200), 1)
    cv2.imshow("Image", temp_display)

    if ctrl_pressed:
        # Remove nearest point from current_line if any
        if current_line:
            dists = [np.hypot(px - x, py - y) for px, py in current_line]
            idx = int(np.argmin(dists))
            removed = current_line.pop(idx)
            print(f"CTRL: Removed point {removed} from current line.")
            # Redraw all points
            img_display[:] = img_original.copy()
            color_idx = current_color
            for line in points_by_line:
                for px, py in line:
                    cv2.circle(img_display, (px, py), 3, colors[color_idx % len(colors)], -1)
                color_idx += 1
            for px, py in current_line:
                cv2.circle(img_display, (px, py), 3, colors[current_color % len(colors)], -1)
            cv2.imshow("Image", img_display)
        return

    # Start drawing on left button down or shift+move
    if event == cv2.EVENT_LBUTTONDOWN or (event == cv2.EVENT_MOUSEMOVE and (drawing or shift_pressed)):
        if event == cv2.EVENT_LBUTTONDOWN:
            drawing = True
        # Only allow selection if (x, y) is on a Canny edge
        if canny_edges is not None and 0 <= y < canny_edges.shape[0] and 0 <= x < canny_edges.shape[1]:
            # Check for nonzero Canny edge pixels within a 2 px radius
            for dy in range(-2, 3):
              for dx in range(-2, 3):
                nx, ny = x + dx, y + dy
                if 0 <= ny < canny_edges.shape[0] and 0 <= nx < canny_edges.shape[1]:
                  if canny_edges[ny, nx] != 0:
                    if (nx, ny) not in current_line:
                      current_line.append((nx, ny))
                      cv2.circle(img_display, (nx, ny), 3, colors[current_color % len(colors)], -1)
    elif event == cv2.EVENT_LBUTTONUP:
        drawing = False
    elif event == cv2.EVENT_RBUTTONDOWN:
        # Finalize current line and start a new one
        if current_line:
            points_by_line.append(current_line.copy())
            current_line.clear()
            print(f"Line {len(points_by_line)} finalized.")
            current_color += 1
            save_lines()  # Save after each line
        img_display[:] = img_original.copy()
        color_idx = 0
        for line in points_by_line:
            for px, py in line:
                cv2.circle(img_display, (px, py), 3, colors[color_idx % len(colors)], -1)
            color_idx += 1
    
    cv2.imshow("Image", img_display)

def apply_radial_distortion(points, k1, center):
    cx, cy = center
    undistorted = []
    for x, y in points:
        dx, dy = x - cx, y - cy
        r2 = dx*dx + dy*dy
        factor = 1 + k1 * r2
        undistorted.append((cx + dx * factor, cy + dy * factor))
    return np.array(undistorted)

def line_fit_error(k1):
    k1 = float(k1)  # Unpack from array to scalar for optimizer compatibility
    total_error = 0
    h, w = img_original.shape[:2]
    center = (w / 2, h / 2)
    for line in points_by_line:
        undistorted = apply_radial_distortion(line, k1, center)
        x = undistorted[:, 0]
        y = undistorted[:, 1]
        A = np.vstack([x, np.ones_like(x)]).T
        m, c = np.linalg.lstsq(A, y, rcond=None)[0]
        y_fit = m * x + c
        total_error += np.sum((y - y_fit) ** 2)
    return total_error

def undistort_image(img, k1):
    h, w = img.shape[:2]
    cx, cy = w / 2, h / 2
    map_x = np.zeros((h, w), dtype=np.float32)
    map_y = np.zeros((h, w), dtype=np.float32)
    for y in range(h):
        for x in range(w):
            dx, dy = x - cx, y - cy
            r2 = dx * dx + dy * dy
            factor = 1 + k1 * r2
            u = cx + dx / factor
            v = cy + dy / factor
            map_x[y, x] = u
            map_y[y, x] = v
    return cv2.remap(img, map_x, map_y, cv2.INTER_LINEAR)

def on_trackbar(val):
    global img_display, canny_edges
    low = cv2.getTrackbarPos("Canny Low", "Canny")
    high = cv2.getTrackbarPos("Canny High", "Canny")
    canny_edges = cv2.Canny(img_original, low, high)
    cv2.imshow("Canny", canny_edges)

def save_lines():
    np.save(LINES_FILE, np.array(points_by_line, dtype=object))

def load_lines():
    global points_by_line
    if os.path.exists(LINES_FILE):
        points_by_line = np.load(LINES_FILE, allow_pickle=True).tolist()
        print(f"Loaded {len(points_by_line)} lines from {LINES_FILE}.")
    else:
        points_by_line.clear()

def main():
    global img_original, img_display, current_color

    path = input("Enter path to image: ").strip()
    img_original = cv2.imread(path)
    if img_original is None:
        print("Could not load image.")
        return

    img_display = img_original.copy()

    load_lines()  # Load lines at start
    current_color = len(points_by_line)
    # Redraw loaded lines
    color_idx = 0
    for line in points_by_line:
        for px, py in line:
            cv2.circle(img_display, (px, py), 3, colors[color_idx % len(colors)], -1)
        color_idx += 1

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
                save_lines()  # Save on exit if line in progress
            break

    if len(points_by_line) < 1:
        print("Not enough data to optimize.")
        return

    print("Optimizing distortion parameter...")
    res = minimize(line_fit_error, x0=0.0, bounds=[(-1e-9, 1e-9)])
    print(res)
    k1_opt = res.x[0]
    print(f"Optimal k1: {k1_opt:.2e}")

    corrected = undistort_image(img_original, k1_opt)
    cv2.imshow("Corrected", corrected)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

if __name__ == "__main__":
    main()
