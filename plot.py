import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse


def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def plot_data(filename):
    try:
        # Load data
        df = pd.read_csv(filename)
        # Strip whitespace from headers just in case
        df.columns = df.columns.str.strip()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return

    print(f"Loaded {len(df)} rows from {filename}")

    # ---------------------------------------------------------
    # 1. TIME SERIES PLOTS (2D)
    # ---------------------------------------------------------
    fig, axes = plt.subplots(5, 1, figsize=(12, 18), sharex=True)
    t = df["time"]

    # --- Plot 1: Position ---
    ax = axes[0]
    if "r_x" in df.columns:
        ax.plot(t, df["r_x"], label="x", linewidth=1)
        ax.plot(t, df["r_y"], label="y", linewidth=1)
        ax.plot(t, df["r_z"], label="z", linewidth=1)
    if "r" in df.columns:
        ax.plot(
            t,
            df["r"],
            label="|r|",
            color="black",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )
    ax.set_ylabel("Position [m]")
    ax.set_title("Orbital Position (ECI)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize="small")

    # --- Plot 2: Velocity ---
    ax = axes[1]
    if "v_x" in df.columns:
        ax.plot(t, df["v_x"], label="vx", linewidth=1)
        ax.plot(t, df["v_y"], label="vy", linewidth=1)
        ax.plot(t, df["v_z"], label="vz", linewidth=1)
    if "v" in df.columns:
        ax.plot(
            t,
            df["v"],
            label="|v|",
            color="black",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )
    ax.set_ylabel("Velocity [m/s]")
    ax.set_title("Orbital Velocity")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize="small")

    # --- Plot 3: Attitude (Quaternion) ---
    ax = axes[2]
    # Check for quaternion columns
    quat_cols = ["q_w", "q_x", "q_y", "q_z"]
    if all(col in df.columns for col in quat_cols):
        ax.plot(t, df["q_w"], label="w", linewidth=1.5)
        ax.plot(t, df["q_x"], label="x", linewidth=1)
        ax.plot(t, df["q_y"], label="y", linewidth=1)
        ax.plot(t, df["q_z"], label="z", linewidth=1)
    ax.set_ylabel("Quaternion")
    ax.set_title("Attitude (Body -> ECI)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize="small")

    # --- Plot 4: Angular Velocity ---
    ax = axes[3]
    if "w_x" in df.columns:
        ax.plot(t, df["w_x"], label="wx", linewidth=1)
        ax.plot(t, df["w_y"], label="wy", linewidth=1)
        ax.plot(t, df["w_z"], label="wz", linewidth=1)
    if "w" in df.columns:
        ax.plot(
            t,
            df["w"],
            label="|w|",
            color="black",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )
    ax.set_ylabel("Ang. Vel [rad/s]")
    ax.set_title("Body Angular Velocity")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize="small")

    # --- Plot 5: Magnetization (Hysteresis) ---
    ax = axes[4]
    # Find all columns starting with M_
    mag_cols = [col for col in df.columns if col.startswith("M_")]
    if mag_cols:
        for col in mag_cols:
            ax.plot(t, df[col], label=col, linewidth=1)
        ax.set_ylabel("Magnetization [A/m]")
        ax.set_title("Rod Magnetization (M_irr)")
        ax.legend(loc="upper right", fontsize="small")
    else:
        ax.text(0.5, 0.5, "No Rod Magnetization Data", ha="center", va="center")

    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Time [s]")

    plt.tight_layout()

    # ---------------------------------------------------------
    # 2. 3D ORBIT VISUALIZATION
    # ---------------------------------------------------------
    if all(col in df.columns for col in ["r_x", "r_y", "r_z"]):
        fig3d = plt.figure(figsize=(10, 8))
        ax3d = fig3d.add_subplot(111, projection="3d")

        # Draw Earth Sphere (approximate representation)
        R_earth = 6378137.0  # meters
        u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
        x_e = R_earth * np.cos(u) * np.sin(v)
        y_e = R_earth * np.sin(u) * np.sin(v)
        z_e = R_earth * np.cos(v)

        # Plot Earth wireframe
        ax3d.plot_wireframe(x_e, y_e, z_e, color="b", alpha=0.1)

        # Plot Orbit
        ax3d.plot(
            df["r_x"], df["r_y"], df["r_z"], color="r", label="Trajectory", linewidth=2
        )

        # Plot Start and End points
        ax3d.scatter(
            df["r_x"].iloc[0],
            df["r_y"].iloc[0],
            df["r_z"].iloc[0],
            c="g",
            marker="o",
            s=50,
            label="Start",
        )
        ax3d.scatter(
            df["r_x"].iloc[-1],
            df["r_y"].iloc[-1],
            df["r_z"].iloc[-1],
            c="black",
            marker="x",
            s=50,
            label="End",
        )

        ax3d.set_xlabel("X [m]")
        ax3d.set_ylabel("Y [m]")
        ax3d.set_zlabel("Z [m]")
        ax3d.set_title("3D Orbit Trajectory (ECI Frame)")
        ax3d.legend()

        # Force aspect ratio to be equal so Earth looks like a sphere
        set_axes_equal(ax3d)

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot AOS Simulation Results")
    parser.add_argument(
        "filename", nargs="?", default="output.csv", help="CSV file to plot"
    )
    args = parser.parse_args()

    plot_data(args.filename)
