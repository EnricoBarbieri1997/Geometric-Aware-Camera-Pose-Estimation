import sys
import struct
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QWidget, QSlider, QPushButton, QLabel, QVBoxLayout, QHBoxLayout, QFileDialog
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class ColmapTransformer(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('COLMAP Points Recenter')

        self.points = None
        self.colors = None
        self.original_points = None

        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111, projection='3d')
        self.canvas = FigureCanvas(self.figure)

        self.load_points()
        self.setup_ui()
        self.plot_points()

    def load_points(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Open COLMAP points3D.bin', '', 'COLMAP files (*.bin *.txt)')
        if path.endswith('.bin'):
            self.points, self.colors = self.read_points3d_bin(path)
        elif path.endswith('.txt'):
            self.points, self.colors = self.read_points3d_txt(path)
        else:
            raise ValueError("Unsupported file format.")

        self.original_points = np.copy(self.points)

    def read_points3d_bin(self, path):
        points = []
        colors = []
        with open(path, 'rb') as f:
            num_points = struct.unpack('<Q', f.read(8))[0]
            for _ in range(num_points):
                point3D_id = struct.unpack('<Q', f.read(8))[0]
                xyz = struct.unpack('<3d', f.read(24))
                rgb = struct.unpack('<3B', f.read(3))
                f.read(1)  # padding
                error = struct.unpack('<d', f.read(8))
                track_length = struct.unpack('<Q', f.read(8))[0]
                f.read(8 * track_length * 2)  # image_id + point2D_idx
                points.append(xyz)
                colors.append(rgb)
        return np.array(points), np.array(colors)

    def read_points3d_txt(self, path):
        points = []
        colors = []
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                elems = line.split()
                xyz = list(map(float, elems[1:4]))
                rgb = list(map(int, elems[4:7]))
                points.append(xyz)
                colors.append(rgb)
        return np.array(points), np.array(colors)

    def setup_ui(self):
        layout = QVBoxLayout()

        layout.addWidget(self.canvas)

        self.sliders = {}
        for axis in ['rx', 'ry', 'rz', 'tx', 'ty', 'tz']:
            hlayout = QHBoxLayout()
            label = QLabel(axis.upper())
            slider = QSlider(Qt.Horizontal)
            slider.setMinimum(-180 if axis.startswith('r') else -10000)
            slider.setMaximum(180 if axis.startswith('r') else 10000)
            slider.setValue(0)
            slider.valueChanged.connect(self.update_transformation)
            hlayout.addWidget(label)
            hlayout.addWidget(slider)
            layout.addLayout(hlayout)
            self.sliders[axis] = slider

        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_points)
        layout.addWidget(self.save_button)

        self.setLayout(layout)

    def rotation_matrix(self, rx, ry, rz):
        rx, ry, rz = np.radians([rx, ry, rz])
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(rx), -np.sin(rx)],
            [0, np.sin(rx), np.cos(rx)]
        ])
        Ry = np.array([
            [np.cos(ry), 0, np.sin(ry)],
            [0, 1, 0],
            [-np.sin(ry), 0, np.cos(ry)]
        ])
        Rz = np.array([
            [np.cos(rz), -np.sin(rz), 0],
            [np.sin(rz), np.cos(rz), 0],
            [0, 0, 1]
        ])
        return Rz @ Ry @ Rx

    def update_transformation(self):
        rx = self.sliders['rx'].value()
        ry = self.sliders['ry'].value()
        rz = self.sliders['rz'].value()
        tx = self.sliders['tx'].value()
        ty = self.sliders['ty'].value()
        tz = self.sliders['tz'].value()

        R = self.rotation_matrix(rx, ry, rz)
        t = np.array([tx/100.0, ty/100.0, tz/100.0])

        self.points = (R @ self.original_points.T).T + t
        self.plot_points()

    def set_axes_equal(self):
        x_limits = self.ax.get_xlim3d()
        y_limits = self.ax.get_ylim3d()
        z_limits = self.ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        plot_radius = 0.5 * max([x_range, y_range, z_range])

        self.ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        self.ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        self.ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    def plot_points(self):
        self.ax.clear()
        normalized_colors = self.colors / 255.0
        self.ax.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2], c=normalized_colors, s=1)
        self.ax.scatter(0, 0, 0, c='red', s=50)
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_title('Points (Z Up)')
        self.ax.view_init(elev=30, azim=45)
        self.set_axes_equal()
        self.canvas.draw()

    def save_points(self):
        path, _ = QFileDialog.getSaveFileName(self, 'Save Transformed Points', '', 'TXT files (*.txt)')
        if path:
            with open(path, 'w') as f:
                for p, c in zip(self.points, self.colors):
                    f.write(f"{p[0]} {p[1]} {p[2]} {int(c[0])} {int(c[1])} {int(c[2])}\n")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    transformer = ColmapTransformer()
    transformer.show()
    sys.exit(app.exec_())
