import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.transform import Rotation as R
import argparse
import sys

# Constants
EARTH_RADIUS = 6371000.0  # meters


def load_data(filename, t_start=None, t_end=None, stride=1):
    try:
        df = pd.read_csv(filename)
        df.columns = df.columns.str.strip()
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        sys.exit(1)

    # Filter Time
    if t_start is not None:
        df = df[df["time"] >= t_start]
    if t_end is not None:
        df = df[df["time"] <= t_end]

    if df.empty:
        print("No data in selected timeframe.")
        sys.exit(1)

    # Downsample for smooth animation
    if len(df) > 1000 and stride == 1:
        stride = len(df) // 500
        print(f"Auto-striding data by factor of {stride} for performance.")

    return df.iloc[::stride].reset_index(drop=True)


def generate_earth_sphere():
    phi = np.linspace(0, 2 * np.pi, 50)
    theta = np.linspace(0, np.pi, 50)
    phi, theta = np.meshgrid(phi, theta)

    x = EARTH_RADIUS * np.sin(theta) * np.cos(phi)
    y = EARTH_RADIUS * np.sin(theta) * np.sin(phi)
    z = EARTH_RADIUS * np.cos(theta)

    return go.Surface(
        x=x, y=y, z=z, colorscale="Earth", showscale=False, opacity=0.8, name="Earth"
    )


def get_rotation_matrices(df):
    quats = df[["q_x", "q_y", "q_z", "q_w"]].to_numpy()
    rotations = R.from_quat(quats)
    return rotations.inv().as_matrix()


def get_nadir_vector(r_eci):
    """
    Returns the normalized vector pointing from the satellite to Earth's center.
    Nadir = -Position / |Position|
    """
    rn = np.linalg.norm(r_eci)
    if rn < 1e-9:
        return np.array([0, 0, -1])  # Fallback
    return -r_eci / rn


def get_cube_mesh(rot_matrix, scale=0.3):
    """
    Returns x, y, z vertices for a rotated cube.
    Scale is relative to the axis length (1.0).
    """
    # Define standard cube vertices (size 2*scale)
    local_pts = (
        np.array(
            [
                [-1, -1, -1],
                [1, -1, -1],
                [1, 1, -1],
                [-1, 1, -1],
                [-1, -1, 1],
                [1, -1, 1],
                [1, 1, 1],
                [-1, 1, 1],
            ]
        )
        * scale
    )

    # Rotate vertices
    rotated_pts = (rot_matrix @ local_pts.T).T

    return rotated_pts[:, 0], rotated_pts[:, 1], rotated_pts[:, 2]


def visualize(filename, t_start, t_end):
    df = load_data(filename, t_start, t_end)

    times = df["time"].values
    pos = df[["r_x", "r_y", "r_z"]].to_numpy()
    rot_mats = get_rotation_matrices(df)

    # --- SETUP FIGURES ---
    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "scene"}, {"type": "scene"}]],
        subplot_titles=("Orbital Trajectory (ECI)", "Attitude (Body Frame vs Nadir)"),
        horizontal_spacing=0.05,
    )

    # -------------------------
    # LEFT PANEL: ORBIT
    # -------------------------
    fig.add_trace(generate_earth_sphere(), row=1, col=1)

    fig.add_trace(
        go.Scatter3d(
            x=pos[:, 0],
            y=pos[:, 1],
            z=pos[:, 2],
            mode="lines",
            line=dict(color="gray", width=2),
            opacity=0.5,
            name="Orbit Path",
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter3d(
            x=[pos[0, 0]],
            y=[pos[0, 1]],
            z=[pos[0, 2]],
            mode="markers",
            marker=dict(size=5, color="cyan"),
            name="Satellite",
        ),
        row=1,
        col=1,
    )

    # -------------------------
    # RIGHT PANEL: ATTITUDE
    # -------------------------
    axis_length = 1.0
    rm0 = rot_mats[0]
    nadir_0 = get_nadir_vector(pos[0]) * axis_length

    # 1. CubeSat Body (Mesh)
    cx, cy, cz = get_cube_mesh(rm0)

    # Define triangles for a cube (0-7 vertices)
    i_idx = [0, 0, 4, 4, 0, 0, 1, 1, 0, 0, 3, 3]
    j_idx = [1, 2, 5, 6, 1, 5, 2, 6, 3, 7, 2, 6]
    k_idx = [2, 3, 6, 7, 5, 4, 6, 5, 7, 4, 6, 2]

    fig.add_trace(
        go.Mesh3d(
            x=cx,
            y=cy,
            z=cz,
            i=i_idx,
            j=j_idx,
            k=k_idx,
            color="gray",
            opacity=0.3,
            name="Chassis",
            flatshading=True,
        ),
        row=1,
        col=2,
    )

    # 2. Axes Helper
    def make_axis(vec, color, name):
        return go.Scatter3d(
            x=[0, vec[0]],
            y=[0, vec[1]],
            z=[0, vec[2]],
            mode="lines+markers",
            line=dict(color=color, width=5),
            marker=dict(size=2),
            name=name,
        )

    # 3. Body Axes
    fig.add_trace(make_axis(rm0 @ [axis_length, 0, 0], "red", "Body X"), row=1, col=2)
    fig.add_trace(make_axis(rm0 @ [0, axis_length, 0], "green", "Body Y"), row=1, col=2)
    fig.add_trace(
        make_axis(rm0 @ [0, 0, axis_length], "cyan", "Body Z (Magnet)"), row=1, col=2
    )

    # 4. Nadir Vector (Replaces B-Field)
    fig.add_trace(make_axis(nadir_0, "magenta", "Nadir (Down)"), row=1, col=2)

    # -------------------------
    # ANIMATION FRAMES
    # -------------------------
    frames = []
    print("Generating animation frames...")

    anim_stride = 1

    for i in range(0, len(df), anim_stride):
        p = pos[i]
        rm = rot_mats[i]
        nadir_vec = get_nadir_vector(p) * axis_length

        bx = rm @ [axis_length, 0, 0]
        by = rm @ [0, axis_length, 0]
        bz = rm @ [0, 0, axis_length]

        # Calculate new Cube vertices
        ncx, ncy, ncz = get_cube_mesh(rm)

        frames.append(
            go.Frame(
                data=[
                    # -- LEFT PANEL --
                    # Trace 2: Sat Position
                    go.Scatter3d(x=[p[0]], y=[p[1]], z=[p[2]]),
                    # -- RIGHT PANEL --
                    # Trace 3: Cube Mesh
                    go.Mesh3d(x=ncx, y=ncy, z=ncz),
                    # Trace 4,5,6: Body Axes
                    go.Scatter3d(x=[0, bx[0]], y=[0, bx[1]], z=[0, bx[2]]),
                    go.Scatter3d(x=[0, by[0]], y=[0, by[1]], z=[0, by[2]]),
                    go.Scatter3d(x=[0, bz[0]], y=[0, bz[1]], z=[0, bz[2]]),
                    # Trace 7: Nadir
                    go.Scatter3d(
                        x=[0, nadir_vec[0]], y=[0, nadir_vec[1]], z=[0, nadir_vec[2]]
                    ),
                ],
                name=str(times[i]),
                traces=[2, 3, 4, 5, 6, 7],
            )
        )

    fig.frames = frames

    # -------------------------
    # LAYOUT
    # -------------------------
    fig.update_layout(
        title=f"AOCS Simulation Visualization ({t_start}s - {t_end}s)",
        template="plotly_dark",
        scene=dict(aspectmode="data"),
        scene2=dict(
            aspectmode="cube",
            xaxis=dict(range=[-1, 1], visible=False),
            yaxis=dict(range=[-1, 1], visible=False),
            zaxis=dict(range=[-1, 1], visible=False),
            annotations=[
                dict(
                    showarrow=False,
                    text="Attitude (Magenta = Nadir)",
                    x=0,
                    y=0,
                    z=1.3,
                    opacity=0.7,
                )
            ],
        ),
        updatemenus=[
            dict(
                type="buttons",
                showactive=False,
                y=0,
                x=0.5,
                xanchor="left",
                yanchor="bottom",
                pad=dict(t=45, r=10),
                buttons=[
                    dict(
                        label="Play",
                        method="animate",
                        args=[
                            None,
                            dict(
                                frame=dict(duration=30, redraw=True), fromcurrent=True
                            ),
                        ],
                    ),
                    dict(
                        label="Pause",
                        method="animate",
                        args=[
                            [None],
                            dict(
                                frame=dict(duration=0, redraw=False),
                                mode="immediate",
                                transition=dict(duration=0),
                            ),
                        ],
                    ),
                ],
            )
        ],
        sliders=[
            dict(
                steps=[
                    dict(
                        method="animate",
                        args=[
                            [str(t)],
                            dict(mode="immediate", frame=dict(duration=0, redraw=True)),
                        ],
                        label=f"{t:.0f}s",
                    )
                    for t in times[::anim_stride]
                ],
                currentvalue=dict(prefix="Time: "),
                pad=dict(t=0),
            )
        ],
    )

    print("Opening browser...")
    fig.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="High-Fidelity 3D AOCS Visualization")
    parser.add_argument("filename", nargs="?", default="output.csv", help="CSV file")
    parser.add_argument("--t-start", type=float, help="Start time [s]")
    parser.add_argument("--t-end", type=float, help="End time [s]")

    args = parser.parse_args()
    visualize(args.filename, args.t_start, args.t_end)
