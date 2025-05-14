import cv2
import json
import os
import numpy as np

WINDOW_NAME = "Primitive Selector"

STATE_FILE = "./real/markers.json"  # Change this to select different state files

IMAGE_PATH = "../test_scenes/markers/"
IMAGES = sorted([
    os.path.join(IMAGE_PATH, f)
    for f in os.listdir(IMAGE_PATH)
    if f.endswith(".png")
])
IMAGE_INDEX = 0

# Initial default positions (can be replaced by loading from file)
default_lines = {
    "red":    [[(100, 100), (200, 100)], [(100, 150), (200, 150)]],
    "green":  [[(300, 100), (400, 100)], [(300, 150), (400, 150)]],
    "blue":   [[(500, 100), (600, 100)], [(500, 150), (600, 150)]],
}
default_ellipses = [[(150, 300), (200, 280), (250, 300), (200, 320), (150, 310)],  # Ellipse 1
                    [(400, 300), (450, 280), (500, 300), (450, 320), (400, 280)]]  # Ellipse 2

dragging_point = None
selected_shape = None

def draw_shapes(img, lines, ellipses):
    for color, line_pairs in lines.items():
        for pt1, pt2 in line_pairs:
            cv2.line(img, pt1, pt2, color_to_bgr(color), 2)
            cv2.circle(img, pt1, 5, color_to_bgr(color), -1)
            cv2.circle(img, pt2, 5, color_to_bgr(color), -1)
    for ellipse in ellipses:
        for pt in ellipse:
            cv2.circle(img, pt, 5, (0, 0, 255), -1)
        if len(ellipse) >= 4:
            cv2.ellipse(img, cv2.fitEllipse(np.array(ellipse, dtype=np.int32)), (0, 0, 255), 1)

def color_to_bgr(color):
    return {"red": (0, 0, 255), "green": (0, 255, 0), "blue": (255, 0, 0)}[color]

def save_shapes(lines, ellipses):
    current_data = [
        {}, {}, {}, {}
    ]
    try:
        with open(STATE_FILE, "r") as f:
            current_data = json.load(f)
    except Exception as e:
        print(e)
    current_data[IMAGE_INDEX] = {
        "lines": {color: [[list(pt1), list(pt2)] for pt1, pt2 in pairs] for color, pairs in lines.items()},
        "ellipses": [[list(pt) for pt in e] for e in ellipses]
    }
    with open(STATE_FILE, "w") as f:
        json.dump(current_data, f, indent=2)
    print("Saved to", STATE_FILE)

def load_shapes():
    if not os.path.exists(STATE_FILE):
        return default_lines.copy(), default_ellipses.copy()
    with open(STATE_FILE, "r") as f:
        try:
            data = json.load(f)[IMAGE_INDEX]
            lines = {color: [(tuple(pt1), tuple(pt2)) for pt1, pt2 in pairs]
                    for color, pairs in data["lines"].items()}
            ellipses = [[tuple(pt) for pt in e] for e in data["ellipses"]]
            return lines, ellipses
        except:
            print("Error reading state file. Using default shapes.")
            return default_lines.copy(), default_ellipses.copy()

def mouse_event(event, x, y, flags, param):
    global dragging_point, selected_shape
    if event == cv2.EVENT_LBUTTONDOWN:
        for color, pairs in param["lines"].items():
            for i, (pt1, pt2) in enumerate(pairs):
                if np.linalg.norm(np.array(pt1) - [x, y]) < 10:
                    dragging_point = (color, i, 0)
                    return
                if np.linalg.norm(np.array(pt2) - [x, y]) < 10:
                    dragging_point = (color, i, 1)
                    return
        for i, ellipse in enumerate(param["ellipses"]):
            for j, pt in enumerate(ellipse):
                if np.linalg.norm(np.array(pt) - [x, y]) < 10:
                    dragging_point = ("ellipse", i, j)
                    return
    elif event == cv2.EVENT_MOUSEMOVE and dragging_point:
        shape, i, j = dragging_point
        if shape == "ellipse":
            param["ellipses"][i][j] = (x, y)
        else:
            pt1, pt2 = param["lines"][shape][i]
            if j == 0:
                param["lines"][shape][i] = ((x, y), pt2)
            else:
                param["lines"][shape][i] = (pt1, (x, y))
    elif event == cv2.EVENT_LBUTTONUP:
        dragging_point = None

def main(image_path):
    lines, ellipses = load_shapes()
    image = cv2.imread(image_path)
    if image is None:
        print("Could not load image:", image_path)
        return
    cv2.namedWindow(WINDOW_NAME)
    cv2.setMouseCallback(WINDOW_NAME, mouse_event, {"lines": lines, "ellipses": ellipses})

    while True:
        display = image.copy()
        draw_shapes(display, lines, ellipses)
        cv2.imshow(WINDOW_NAME, display)
        key = cv2.waitKey(20)
        if key == 27:  # ESC
            break
        elif key == ord("s"):
            save_shapes(lines, ellipses)
        elif key == ord("r"):
            lines, ellipses = load_shapes()

    cv2.destroyAllWindows()

if __name__ == "__main__":
    for index, image_path in enumerate(IMAGES):
        IMAGE_INDEX = index
        print(f"Processing: {image_path} with STATE_FILE: {STATE_FILE}")
        main(image_path)
