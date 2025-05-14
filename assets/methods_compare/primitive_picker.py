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

ZOOM_STEP = 0.1
MIN_ZOOM = 0.2
MAX_ZOOM = 5.0
zoom = 1.0
offset = np.array([0, 0], dtype=np.float32)
pan_start = None

def to_image_coords(x, y):
    return int((x - offset[0]) / zoom), int((y - offset[1]) / zoom)

def to_display_coords(x, y):
    return int(x * zoom + offset[0]), int(y * zoom + offset[1])

def draw_shapes_zoomed(img, lines, ellipses):
    for color, line_pairs in lines.items():
        lines = []
        points = []
        for pt1, pt2 in line_pairs:
            pt1_disp = to_display_coords(*pt1)
            pt2_disp = to_display_coords(*pt2)
            lines.append(np.cross(pt1 + (1.0, ), pt2 + (1.0, )))
            points.append(pt1)
            cv2.line(img, pt1_disp, pt2_disp, color_to_bgr(color), 2)
            cv2.circle(img, pt1_disp, 5, color_to_bgr(color), -1)
            cv2.circle(img, pt2_disp, 5, color_to_bgr(color), -1)
        vanishing_point = np.cross(lines[0], lines[1])
        vanishing_point /= vanishing_point[2]
        vanishing_point = (vanishing_point[0], vanishing_point[1])
        cv2.circle(img, to_display_coords(*vanishing_point), 5, color_to_bgr(color), -1)
        for pt in points:
            cv2.line(img, to_display_coords(*pt), to_display_coords(*vanishing_point), color_to_bgr(color), 2, lineType=4)
    # for ellipse in ellipses:
    #     pts_disp = np.array([to_display_coords(*pt) for pt in ellipse])
    #     for pt in pts_disp:
    #         cv2.circle(img, pt, 5, (0, 0, 255), -1)
    #     if len(pts_disp) >= 4:
    #         cv2.ellipse(img, cv2.fitEllipse(pts_disp.astype(np.int32)), (0, 0, 255), 1)

def get_viewport(image_shape, offset, zoom, window_size):
    # Convert screen offset to image space
    x0 = int(max(0, (0 - offset[0]) / zoom))
    y0 = int(max(0, (0 - offset[1]) / zoom))

    # Compute bottom-right corner in image space
    x1 = int(min(image_shape[1], (window_size[0] - offset[0]) / zoom))
    y1 = int(min(image_shape[0], (window_size[1] - offset[1]) / zoom))

    return x0, y0, x1, y1

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
    global dragging_point, selected_shape, zoom, offset, pan_start

    if flags & cv2.EVENT_FLAG_CTRLKEY and event == cv2.EVENT_MOUSEWHEEL:
        old_zoom = zoom
        zoom = np.clip(zoom + (ZOOM_STEP if flags > 0 else -ZOOM_STEP), MIN_ZOOM, MAX_ZOOM)
        # Adjust offset so zoom focuses on cursor
        cursor_img_x, cursor_img_y = to_image_coords(x, y)
        offset = np.array([x - cursor_img_x * zoom, y - cursor_img_y * zoom], dtype=np.float32)

    elif event == cv2.EVENT_RBUTTONDOWN:
        pan_start = np.array([x, y], dtype=np.float32)

    elif event == cv2.EVENT_MOUSEMOVE and pan_start is not None:
        delta = np.array([x, y], dtype=np.float32) - pan_start
        offset += delta
        pan_start = np.array([x, y], dtype=np.float32)

    elif event == cv2.EVENT_RBUTTONUP:
        pan_start = None

    elif event == cv2.EVENT_LBUTTONDOWN:
        ix, iy = to_image_coords(x, y)
        for color, pairs in param["lines"].items():
            for i, (pt1, pt2) in enumerate(pairs):
                if np.linalg.norm(np.array(pt1) - [ix, iy]) < 10:
                    dragging_point = (color, i, 0)
                    return
                if np.linalg.norm(np.array(pt2) - [ix, iy]) < 10:
                    dragging_point = (color, i, 1)
                    return
        for i, ellipse in enumerate(param["ellipses"]):
            for j, pt in enumerate(ellipse):
                if np.linalg.norm(np.array(pt) - [ix, iy]) < 10:
                    dragging_point = ("ellipse", i, j)
                    return

    elif event == cv2.EVENT_MOUSEMOVE and dragging_point:
        ix, iy = to_image_coords(x, y)
        shape, i, j = dragging_point
        if shape == "ellipse":
            param["ellipses"][i][j] = (ix, iy)
        else:
            pt1, pt2 = param["lines"][shape][i]
            if j == 0:
                param["lines"][shape][i] = ((ix, iy), pt2)
            else:
                param["lines"][shape][i] = (pt1, (ix, iy))

    elif event == cv2.EVENT_LBUTTONUP:
        dragging_point = None

def main(image_path):
    global zoom, offset
    lines, ellipses = load_shapes()
    image = cv2.imread(image_path)
    if image is None:
        print("Could not load image:", image_path)
        return
    zoom = 1.0
    offset = np.array([0, 0], dtype=np.float32)

    cv2.namedWindow(WINDOW_NAME)
    cv2.setMouseCallback(WINDOW_NAME, mouse_event, {"lines": lines, "ellipses": ellipses})

    while True:
        window_size = [image.shape[1], image.shape[0]]
        x0, y0, x1, y1 = get_viewport(image.shape, offset, zoom, window_size)
        viewport = image[y0:y1, x0:x1]

        # Resize to window size for display
        canvas = cv2.resize(viewport, window_size, interpolation=cv2.INTER_LINEAR)
        draw_shapes_zoomed(canvas, lines, ellipses)
        cv2.imshow(WINDOW_NAME, canvas)
        key = cv2.waitKey(20)
        if key == 27:  # ESC
            break
        elif key == ord("s"):
            save_shapes(lines, ellipses)
        elif key == ord("r"):
            lines, ellipses = load_shapes()
        elif key == ord('+') or key == ord('='):
            # Zoom in
            zoom = min(zoom + ZOOM_STEP, MAX_ZOOM)
        elif key == ord('-') or key == ord('_'):
            # Zoom out
            zoom = max(zoom - ZOOM_STEP, MIN_ZOOM)

    cv2.destroyAllWindows()

if __name__ == "__main__":
    for index, image_path in enumerate(IMAGES):
        IMAGE_INDEX = index
        print(f"Processing: {image_path} with STATE_FILE: {STATE_FILE}")
        main(image_path)
