import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys


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

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def plot_data(filename, plot_types, display_mode, t_start=None, t_end=None):
    try:
        df = pd.read_csv(filename)
        df.columns = df.columns.str.strip()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    # ---------------------------------------------------------
    # 0. DATA FILTERING (TIMEFRAME)
    # ---------------------------------------------------------
    original_count = len(df)

    if t_start is not None:
        df = df[df["time"] >= t_start]

    if t_end is not None:
        df = df[df["time"] <= t_end]

    if df.empty:
        print(
            f"Error: No data found in the specified time range ({t_start} to {t_end})."
        )
        sys.exit(1)

    print(f"Loaded {original_count} rows from {filename}")
    print(
        f"Plotting {len(df)} rows ({df['time'].iloc[0]:.2f}s to {df['time'].iloc[-1]:.2f}s)"
    )

    t = df["time"]

    # Filter Plot Types
    show_all = "all" in plot_types
    do_pos = show_all or "pos" in plot_types
    do_vel = show_all or "vel" in plot_types
    do_att = show_all or "att" in plot_types
    do_omega = show_all or "omega" in plot_types
    do_mag = show_all or "mag" in plot_types
    do_3d = show_all or "3d" in plot_types

    # Filter Display Mode (Components vs Magnitude)
    show_comp = display_mode in ["comp", "both"]
    show_mag = display_mode in ["mag", "both"]

    # Calculate number of 2D subplots needed
    active_2d_plots = [do_pos, do_vel, do_att, do_omega, do_mag]
    n_rows = sum(active_2d_plots)

    # ---------------------------------------------------------
    # 1. TIME SERIES PLOTS (2D)
    # ---------------------------------------------------------
    if n_rows > 0:
        # Create figure with dynamic height
        fig, axes = plt.subplots(n_rows, 1, figsize=(12, 3.5 * n_rows), sharex=True)

        # Ensure axes is iterable even if there is only 1 subplot
        if n_rows == 1:
            axes = [axes]

        current_ax_idx = 0

        # --- Plot: Position ---
        if do_pos:
            ax = axes[current_ax_idx]
            if show_comp and "r_x" in df.columns:
                ax.plot(t, df["r_x"], label="x", linewidth=1)
                ax.plot(t, df["r_y"], label="y", linewidth=1)
                ax.plot(t, df["r_z"], label="z", linewidth=1)
            if show_mag and "r" in df.columns:
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
            current_ax_idx += 1

        # --- Plot: Velocity ---
        if do_vel:
            ax = axes[current_ax_idx]
            if show_comp and "v_x" in df.columns:
                ax.plot(t, df["v_x"], label="vx", linewidth=1)
                ax.plot(t, df["v_y"], label="vy", linewidth=1)
                ax.plot(t, df["v_z"], label="vz", linewidth=1)
            if show_mag and "v" in df.columns:
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
            current_ax_idx += 1

        # --- Plot: Attitude ---
        if do_att:
            ax = axes[current_ax_idx]
            quat_cols = ["q_w", "q_x", "q_y", "q_z"]
            if show_comp and all(col in df.columns for col in quat_cols):
                ax.plot(t, df["q_w"], label="w", linewidth=1.5)
                ax.plot(t, df["q_x"], label="x", linewidth=1)
                ax.plot(t, df["q_y"], label="y", linewidth=1)
                ax.plot(t, df["q_z"], label="z", linewidth=1)

            if show_mag:
                q_mag = np.sqrt(
                    df["q_w"] ** 2 + df["q_x"] ** 2 + df["q_y"] ** 2 + df["q_z"] ** 2
                )
                ax.plot(
                    t,
                    q_mag,
                    label="|q|",
                    color="black",
                    linestyle="--",
                    linewidth=1,
                    alpha=0.7,
                )

            ax.set_ylabel("Quaternion")
            ax.set_title("Attitude (Body -> ECI)")
            ax.grid(True, alpha=0.3)
            ax.legend(loc="upper right", fontsize="small")
            current_ax_idx += 1

        # --- Plot: Angular Velocity ---
        if do_omega:
            ax = axes[current_ax_idx]
            if show_comp and "w_x" in df.columns:
                ax.plot(t, df["w_x"], label="wx", linewidth=1)
                ax.plot(t, df["w_y"], label="wy", linewidth=1)
                ax.plot(t, df["w_z"], label="wz", linewidth=1)
            if show_mag and "w" in df.columns:
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
            current_ax_idx += 1

        # --- Plot: Magnetization ---
        if do_mag:
            ax = axes[current_ax_idx]
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
            current_ax_idx += 1

        axes[-1].set_xlabel("Time [s]")
        plt.tight_layout()

    # ---------------------------------------------------------
    # 2. 3D ORBIT VISUALIZATION
    # ---------------------------------------------------------
    if do_3d:
        if all(col in df.columns for col in ["r_x", "r_y", "r_z"]):
            fig3d = plt.figure(figsize=(10, 8))
            ax3d = fig3d.add_subplot(111, projection="3d")

            # Draw Earth Sphere
            R_earth = 6378137.0
            u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
            x_e = R_earth * np.cos(u) * np.sin(v)
            y_e = R_earth * np.sin(u) * np.sin(v)
            z_e = R_earth * np.cos(v)

            ax3d.plot_wireframe(x_e, y_e, z_e, color="b", alpha=0.1)
            ax3d.plot(
                df["r_x"],
                df["r_y"],
                df["r_z"],
                color="r",
                label="Trajectory",
                linewidth=2,
            )

            # Markers relative to the filtered dataset
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
            ax3d.set_title(
                f"3D Orbit Trajectory ({df['time'].iloc[0]:.0f}s - {df['time'].iloc[-1]:.0f}s)"
            )
            ax3d.legend()
            set_axes_equal(ax3d)
        else:
            print("Warning: Cannot plot 3D Orbit. Missing r_x, r_y, or r_z columns.")

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot AOS Simulation Results")

    parser.add_argument(
        "filename",
        nargs="?",
        default="output.csv",
        help="CSV file to plot (default: output.csv)",
    )

    parser.add_argument(
        "-t",
        "--types",
        nargs="+",
        default=["all"],
        choices=["pos", "vel", "att", "omega", "mag", "3d", "all"],
        help="Select which graphs to display (default: all)",
    )

    parser.add_argument(
        "-d",
        "--display",
        default="both",
        choices=["comp", "mag", "both"],
        help="Display components (x,y,z), magnitude (|v|), or both (default: both)",
    )

    parser.add_argument(  #
        "-s", "--t-start", type=float, help="Start time for plotting [s]"
    )

    parser.add_argument(  #
        "-e", "--t-end", type=float, help="End time for plotting [s]"
    )

    args = parser.parse_args()
    plot_data(args.filename, args.types, args.display, args.t_start, args.t_end)
