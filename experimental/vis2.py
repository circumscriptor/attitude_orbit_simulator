import taichi as ti
import pyvista as pv
from pyvista import examples
import pandas as pd
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse
import sys
import math

# Initialize Taichi
try:
    ti.init(arch=ti.vulkan)
except:
    ti.init(arch=ti.cpu)

# --- CONSTANTS ---
EARTH_RADIUS_METERS = 6371000.0
SCALE = 1.0 / EARTH_RADIUS_METERS
AXIS_LEN = 0.3  # Length of the colored axes lines
SAT_SIZE = 0.08  # Physical size of the satellite box (relative to Earth=1.0)
EARTH_OMEGA = 7.2921159e-5


@ti.data_oriented
class AOCSVisualizer:
    def __init__(self, df, initial_angle_deg=0.0, rotation_speed=1.0):
        self.n_frames = len(df)
        self.initial_angle = math.radians(initial_angle_deg)
        self.rotation_speed = rotation_speed

        print(f"Initializing Visualizer: {self.n_frames} frames")

        # --- 1. DATA PROCESSING ---
        self.timestamps = df["time"].to_numpy().astype(np.float32)

        # Position (Scaled)
        pos_np = df[["r_x", "r_y", "r_z"]].to_numpy().astype(np.float32) * SCALE

        # Orientation
        quats = df[["q_x", "q_y", "q_z", "q_w"]].to_numpy()
        self.rotations = R.from_quat(quats)

        # Pre-compute Body Axes (Scaled for Line Drawing)
        # We also need normalized vectors for generating the box geometry
        self.axis_x = self.rotations.apply([1, 0, 0]).astype(np.float32)
        self.axis_y = self.rotations.apply([0, 1, 0]).astype(np.float32)
        self.axis_z = self.rotations.apply([0, 0, 1]).astype(np.float32)

        # Store scaled versions for the "Gizmo" (colored lines)
        self.axis_x_scaled = self.axis_x * AXIS_LEN
        self.axis_y_scaled = self.axis_y * AXIS_LEN
        self.axis_z_scaled = self.axis_z * AXIS_LEN

        # Nadir Vector
        norms = np.linalg.norm(pos_np, axis=1, keepdims=True)
        self.nadir = (-pos_np / np.where(norms == 0, 1, norms)).astype(
            np.float32
        ) * AXIS_LEN

        # --- 2. TAICHI FIELDS (Trajectory) ---
        self.pos_field = ti.Vector.field(3, dtype=ti.f32, shape=self.n_frames)
        self.pos_field.from_numpy(pos_np)

        if self.n_frames > 1:
            self.orbit_indices = ti.field(dtype=ti.i32, shape=(self.n_frames - 1) * 2)
            indices = np.repeat(np.arange(self.n_frames - 1), 2)
            indices[1::2] += 1
            self.orbit_indices.from_numpy(indices.astype(np.int32))

        # --- 3. SATELLITE GEOMETRY (Box & Axes) ---

        # A. Axes Lines (8 verts)
        self.vec_verts = ti.Vector.field(3, dtype=ti.f32, shape=8)
        self.vec_colors = ti.Vector.field(3, dtype=ti.f32, shape=8)
        self.vec_indices = ti.field(dtype=ti.i32, shape=8)

        # Colors: X=Red, Y=Green, Z=Cyan, Nadir=Magenta
        axis_colors = np.array(
            [
                [1, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [0, 1, 0],
                [0, 1, 1],
                [0, 1, 1],
                [1, 0, 1],
                [1, 0, 1],
            ],
            dtype=np.float32,
        )
        self.vec_colors.from_numpy(axis_colors)
        self.vec_indices.from_numpy(np.arange(8, dtype=np.int32))

        self.box_size = 0.05

        # 1. Static Local Vertices (Computed ONCE)
        self.box_verts_local = ti.Vector.field(3, dtype=ti.f32, shape=8)
        # 2. Dynamic World Vertices (Updated every frame on GPU)
        self.box_verts_world = ti.Vector.field(3, dtype=ti.f32, shape=8)

        self.box_indices = ti.field(dtype=ti.i32, shape=36)

        # Initialize Local Box (Centered at 0,0,0)
        hs = self.box_size * 0.5
        local_corners = np.array(
            [
                [-hs, -hs, -hs],
                [hs, -hs, -hs],
                [hs, hs, -hs],
                [-hs, hs, -hs],
                [-hs, -hs, hs],
                [hs, -hs, hs],
                [hs, hs, hs],
                [-hs, hs, hs],
            ],
            dtype=np.float32,
        )

        # Indices (Same as before)
        indices_np = np.array(
            [
                0,
                2,
                1,
                0,
                3,
                2,
                4,
                5,
                6,
                4,
                6,
                7,
                0,
                1,
                5,
                0,
                5,
                4,
                1,
                2,
                6,
                1,
                6,
                5,
                2,
                3,
                7,
                2,
                7,
                6,
                3,
                0,
                4,
                3,
                4,
                7,
            ],
            dtype=np.int32,
        )

        self.box_verts_local.from_numpy(local_corners)
        self.box_indices.from_numpy(indices_np)

        # --- 4. EARTH MESH ---
        self._init_earth_mesh()

    def _init_earth_mesh(self):
        print("Generating Earth mesh...")
        earth = examples.planets.load_earth(radius=1.0)
        earth = earth.triangulate()
        verts = np.asarray(earth.points, dtype=np.float32)
        faces = earth.faces.reshape(-1, 4)[:, 1:].astype(np.int32)

        tex_obj = examples.load_globe_texture()
        tex_data = np.asarray(tex_obj.to_array())
        if tex_data.ndim == 3 and tex_data.shape[2] == 4:
            tex_data = tex_data[:, :, :3]
        tex_data = tex_data.astype(np.float32) / 255.0
        tex_data = np.flipud(tex_data)
        tex_data = np.transpose(tex_data, (1, 0, 2))

        # UV generation
        x, y, z = verts[:, 0], verts[:, 1], verts[:, 2]
        u = 0.5 + np.arctan2(y, x) / (2 * np.pi)
        v = 0.5 + np.arcsin(np.clip(z, -1.0, 1.0)) / np.pi
        uvs = np.stack([u, v], axis=-1).astype(np.float32)

        self.n_earth_verts = len(verts)
        self.earth_verts_ref = ti.Vector.field(
            3, dtype=ti.f32, shape=self.n_earth_verts
        )
        self.earth_verts_render = ti.Vector.field(
            3, dtype=ti.f32, shape=self.n_earth_verts
        )
        self.earth_indices = ti.field(dtype=ti.i32, shape=faces.size)
        self.earth_uvs = ti.Vector.field(2, dtype=ti.f32, shape=len(uvs))
        self.earth_colors = ti.Vector.field(3, dtype=ti.f32, shape=len(verts))
        self.earth_tex = ti.Vector.field(
            3, dtype=ti.f32, shape=(tex_data.shape[0], tex_data.shape[1])
        )

        self.earth_verts_ref.from_numpy(verts)
        self.earth_verts_render.from_numpy(verts)
        self.earth_indices.from_numpy(faces.flatten())
        self.earth_uvs.from_numpy(uvs)
        self.earth_tex.from_numpy(tex_data)
        self.bake_texture()

    @ti.func
    def quat_rotate(self, q: ti.types.vector(4, float), v: ti.types.vector(3, float)):
        # Standard Quaternion Rotation Formula
        q_vec = ti.Vector([q.x, q.y, q.z])
        uv = q_vec.cross(v)
        uuv = q_vec.cross(uv)
        return v + 2.0 * (q.w * uv + uuv)

    @ti.kernel
    def update_box_kernel(
        self, pos: ti.types.vector(3, float), q: ti.types.vector(4, float)
    ):
        # Parallelize over the 8 vertices
        for i in self.box_verts_local:
            local_v = self.box_verts_local[i]
            rotated_v = self.quat_rotate(q, local_v)
            self.box_verts_world[i] = rotated_v + pos

    @ti.kernel
    def bake_texture(self):
        for i in self.earth_verts_ref:
            uv = self.earth_uvs[i]
            w, h = float(self.earth_tex.shape[0]), float(self.earth_tex.shape[1])
            tx, ty = int(uv.x * (w - 1)), int(uv.y * (h - 1))
            self.earth_colors[i] = self.earth_tex[tx, ty]

    @ti.kernel
    def rotate_earth_kernel(self, angle: float):
        c, s = ti.cos(angle), ti.sin(angle)
        for i in self.earth_verts_render:
            p = self.earth_verts_ref[i]
            self.earth_verts_render[i] = ti.Vector(
                [p.x * c - p.y * s, p.x * s + p.y * c, p.z]
            )

    def update(self, idx):
        # 1. Get raw data
        pos_vec = self.pos_field[idx]

        # Get quaternion from your DataFrame (assumes you stored it or access it here)
        # Note: You might need to create a Field for quats if you haven't,
        # or just pass it as a value if it's single frame.
        # Assuming you have the dataframe 'df' available or stored quats in a field:
        # For efficiency, best to store quats in a Taichi field like pos_field.
        # But here is the Python-to-Kernel pass for simplicity:

        q_np = self.rotations[idx].as_quat()  # SciPy returns [x,y,z,w]
        q_vec = ti.Vector([q_np[0], q_np[1], q_np[2], q_np[3]])

        # 2. Run GPU Kernel
        self.update_box_kernel(pos_vec, q_vec)

        # 3. Update Earth Rotation
        if self.rotation_speed != 0.0:
            t = self.timestamps[idx]
            angle = self.initial_angle + (EARTH_OMEGA * t * self.rotation_speed)
            self.rotate_earth_kernel(angle)

    def run(self):
        window = ti.ui.Window("AOCS Visualizer", (1280, 720), vsync=True)
        canvas = window.get_canvas()
        scene = window.get_scene()
        camera = ti.ui.Camera()

        cam_dist, cam_yaw, cam_pitch = 3.5, 0.0, 1.0
        last_mouse = window.get_cursor_pos()
        current_frame, is_playing, speed = 0.0, True, 1.0

        while window.running:
            # Camera Input
            curr_mouse = window.get_cursor_pos()
            if window.is_pressed(ti.ui.LMB):
                dx, dy = curr_mouse[0] - last_mouse[0], curr_mouse[1] - last_mouse[1]
                cam_yaw -= dx * 5.0
                cam_pitch = max(0.01, min(cam_pitch + dy * 5.0, 3.13))
            if window.is_pressed(ti.ui.RMB):
                cam_dist = max(
                    1.1, min(cam_dist + (curr_mouse[1] - last_mouse[1]) * 5.0, 50.0)
                )
            last_mouse = curr_mouse

            camera.position(
                cam_dist * math.sin(cam_pitch) * math.cos(cam_yaw),
                cam_dist * math.sin(cam_pitch) * math.sin(cam_yaw),
                cam_dist * math.cos(cam_pitch),
            )
            camera.lookat(0, 0, 0)
            camera.up(0, 0, 1)
            scene.set_camera(camera)

            # Update Logic
            if is_playing:
                current_frame += speed
                if current_frame >= self.n_frames:
                    current_frame = 0

            idx = int(max(0, min(current_frame, self.n_frames - 1)))
            self.update(idx)

            # Rendering
            scene.ambient_light((0.3, 0.3, 0.3))
            scene.point_light(pos=(10, 10, 10), color=(1, 1, 1))

            # 1. Earth
            scene.mesh(
                self.earth_verts_render,
                indices=self.earth_indices,
                per_vertex_color=self.earth_colors,
            )

            # 2. Orbit Trajectory
            if self.n_frames > 1:
                scene.lines(
                    self.pos_field,
                    indices=self.orbit_indices,
                    color=(0.7, 0.7, 0.7),
                    width=1.0,
                )

            # 3. Satellite Box (Gold Color)
            scene.mesh(
                self.box_verts_world, indices=self.box_indices, color=(1.0, 0.8, 0.2)
            )

            # 4. Satellite Axes (Gizmo)
            scene.lines(
                self.vec_verts,
                indices=self.vec_indices,
                per_vertex_color=self.vec_colors,
                width=3.0,
            )

            # GUI
            window.GUI.begin("Simulation", 0.05, 0.05, 0.25, 0.25)
            if window.GUI.button("Pause/Play"):
                is_playing = not is_playing
            speed = window.GUI.slider_float("Speed", speed, 0.1, 5.0)
            prev = current_frame
            current_frame = window.GUI.slider_float(
                "Time", float(current_frame), 0.0, float(self.n_frames - 1)
            )
            if abs(current_frame - prev) > speed + 0.1:
                is_playing = False
            window.GUI.text(f"Time: {self.timestamps[idx]:.1f} s")
            window.GUI.end()

            canvas.scene(scene)
            window.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", nargs="?", default="output.csv")
    parser.add_argument("--t-start", type=float)
    parser.add_argument("--t-end", type=float)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--initial-angle", type=float, default=0.0)
    parser.add_argument("--rotation-speed", type=float, default=1.0)
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.filename)
        df.columns = df.columns.str.strip()
        if args.t_start:
            df = df[df["time"] >= args.t_start]
        if args.t_end:
            df = df[df["time"] <= args.t_end]
        if len(df) == 0:
            sys.exit(1)
        if len(df) > 20000 and args.stride == 1:
            args.stride = len(df) // 10000
        df = df.iloc[:: args.stride].reset_index(drop=True)

        vis = AOCSVisualizer(
            df, initial_angle_deg=args.initial_angle, rotation_speed=args.rotation_speed
        )
        vis.run()
    except Exception as e:
        print(f"Error: {e}")
