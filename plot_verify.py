from pathlib import Path
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import platform

# Fix for "QSocketNotifier" and decoration warnings: Use the TkAgg backend
matplotlib.use('TkAgg')

def maximize_window():
    """Utility to maximize the matplotlib window based on the OS."""
    manager = plt.get_current_fig_manager()
    system = platform.system()
    try:
        if system == "Linux":
            manager.window.attributes('-zoomed', True)
        elif system == "Windows":
            manager.window.state('zoomed')
        elif system == "Darwin":  # macOS
            manager.frame.Maximize(True)
    except Exception:
        # Fallback if the specific backend manager doesn't support zooming
        pass

def plot_torques(t, df, save=False, prefix="output", dpi=150):
    """Figure 1: Comparison of all disturbance and control torques."""
    groups = [
        (['t_mag_x', 't_mag_y', 't_mag_z'], "Permanent Magnet Torque [Nm]"),
        (['t_rods_x', 't_rods_y', 't_rods_z'], "Hysteresis Rod Torque [Nm]"),
        (['t_grav_x', 't_grav_y', 't_grav_z'], "Gravity Gradient Torque [Nm]"),
        (['t_gyro_x', 't_gyro_y', 't_gyro_z'], "Gyroscopic Torque [Nm]"),
        (['t_face_x', 't_face_y', 't_face_z'], "Surface Effects Torque (Drag/SRP) [Nm]")
    ]

    fig, axs = plt.subplots(len(groups), 1, figsize=(14, 10), sharex=True)
    plt.suptitle("Spacecraft Torque Analysis (Body Frame)", fontsize=14, fontweight='bold')

    for ax, (cols, title) in zip(axs, groups):
        # Filter only columns that exist in the dataframe
        existing_cols = [c for c in cols if c in df.columns]
        if existing_cols:
            ax.plot(t, df[existing_cols], lw=0.8)

        ax.set_title(title, fontsize=10, pad=5)
        ax.grid(True, alpha=0.3)
        ax.legend(['X', 'Y', 'Z'], loc='upper right', fontsize='small', ncol=3)

    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save:
        filename = f"{prefix}_torques.png"
        plt.savefig(filename, dpi=dpi)
        plt.close()
        print(f"Saved: {filename} at {dpi} DPI")
    else:
        maximize_window()

def plot_attitude_state(t, df, save=False, prefix="output", dpi=150):
    """Figure 2: Satellite rotational state and rod magnetizations."""
    fig, axs = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("Attitude Dynamics and Internal State", fontsize=14, fontweight='bold')

    # 1. Angular Velocity
    axs[0].plot(t, df[['w_x', 'w_y', 'w_z']], lw=1)
    axs[0].set_title("Angular Velocity [rad/s]")
    axs[0].legend(['ωx', 'ωy', 'ωz'], loc='upper right', ncol=3)

    # 2. Quaternions
    axs[1].plot(t, df[['q_w', 'q_x', 'q_y', 'q_z']], lw=1)
    axs[1].set_title("Attitude Quaternions")
    axs[1].legend(['qw', 'qx', 'qy', 'qz'], loc='upper right', ncol=4)

    # 3. Rod Magnetization
    rod_cols = [c for c in df.columns if c.startswith('M_')]
    if rod_cols:
        axs[2].plot(t, df[rod_cols], lw=1)
        axs[2].set_title("Rod Magnetizations [Am^2]")
        # Limit columns in legend to keep it readable
        axs[2].legend(rod_cols, loc='upper right', ncol=min(len(rod_cols), 4))

    for ax in axs:
        ax.grid(True, alpha=0.3)

    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save:
        filename = f"{prefix}_attitude_state.png"
        plt.savefig(filename, dpi=dpi)
        plt.close()
        print(f"Saved: {filename} at {dpi} DPI")
    else:
        maximize_window()

def plot_environment(t, df, save=False, prefix="output", dpi=150):
    """Figure 3: External environment conditions."""
    fig, axs = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("Environmental Verification", fontsize=14, fontweight='bold')

    # 1. Magnetic Field Magnitude (Vectorized norm)
    mag_norm = np.linalg.norm(df[['mag_x', 'mag_y', 'mag_z']], axis=1)
    axs[0].plot(t, mag_norm, color='blue', label='|B| Total')
    axs[0].set_title("Magnetic Field Magnitude [Tesla]")
    axs[0].legend(loc='upper right')

    # 2. Shadow and Solar Pressure
    axs[1].fill_between(t, 0, df['shadow'], alpha=0.2, color='orange', label='In Sun (1=Yes)')
    axs[1].set_ylabel("Shadow Factor", color='orange')

    ax_srp = axs[1].twinx()
    ax_srp.plot(t, df['solar_p'], color='red', lw=1, label='Solar Pressure')
    ax_srp.set_ylabel("SRP [Pa]", color='red')
    axs[1].set_title("Solar Context")

    # Merge legends for twin axes
    h1, l1 = axs[1].get_legend_handles_labels()
    h2, l2 = ax_srp.get_legend_handles_labels()
    axs[1].legend(h1 + h2, l1 + l2, loc='upper right')

    # 3. Atmospheric Density (Log scale)
    axs[2].semilogy(t, df['rho'], color='green', label='Density (ρ)')
    axs[2].set_title("Atmospheric Density [kg/m^3]")
    axs[2].legend(loc='upper right')

    for ax in axs: ax.grid(True, alpha=0.3)
    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save:
        filename = f"{prefix}_environment.png"
        plt.savefig(filename, dpi=dpi)
        plt.close()
        print(f"Saved: {filename} at {dpi} DPI")
    else:
        maximize_window()

def plot_3d_orbit(df, save=False, prefix="output", dpi=150):
    """Figure 4: 3D Orbital Path and Nadir Visualization."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    plt.suptitle("3D Orbital Shape and Nadir Vectors", fontsize=14, fontweight='bold')

    rx, ry, rz = df['r_x'].values, df['r_y'].values, df['r_z'].values

    # 1. Plot the orbital trajectory
    ax.plot(rx, ry, rz, color='blue', lw=1.5, label='Orbital Path')

    # 2. Plot Earth wireframe
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    earth_r = 6371000
    x = earth_r * np.cos(u) * np.sin(v)
    y = earth_r * np.sin(u) * np.sin(v)
    z = earth_r * np.cos(v)
    ax.plot_wireframe(x, y, z, color="green", alpha=0.1, label='Earth')

    # 3. Draw Nadir vectors (Vectorized sampling)
    interval = max(1, len(df) // 20)
    # Using indexing to avoid the Python for-loop
    ax.quiver(rx[::interval], ry[::interval], rz[::interval],
              -rx[::interval], -ry[::interval], -rz[::interval],
              color='red', length=0.8, arrow_length_ratio=0.1,
              alpha=0.6, label='Nadir')

    ax.set_xlabel('ECI X [m]')
    ax.set_ylabel('ECI Y [m]')
    ax.set_zlabel('ECI Z [m]')

    # 4. Aspect Ratio and Legend
    ax.set_aspect('equal')
    ax.legend(loc='upper right')
    plt.tight_layout()

    if save:
        filename = f"{prefix}_orbit_3d.png"
        plt.savefig(filename, dpi=dpi)
        plt.close()
        print(f"Saved: {filename} at {dpi} DPI")
    else:
        maximize_window()

def plot_nadir_pointing_error(t, df, save=False, prefix="output", dpi=150):
    """Figure 5: Nadir Pointing Error using vectorized math."""
    plt.figure(figsize=(14, 6))

    # Vectorized normalization of Nadir ECI (-R)
    r_vec = df[['r_x', 'r_y', 'r_z']].values
    nadir_eci = -r_vec / np.linalg.norm(r_vec, axis=1)[:, np.newaxis]

    # Vectorized transformation of Body Z-axis [0,0,1] to ECI
    qw, qx, qy, qz = df['q_w'], df['q_x'], df['q_y'], df['q_z']

    z_body_eci = np.stack([
        2 * (qx*qz + qw*qy),
        2 * (qy*qz - qw*qx),
        qw**2 - qx**2 - qy**2 + qz**2
    ], axis=1)

    # Vectorized dot product and angle calculation
    cos_theta = np.einsum('ij,ij->i', z_body_eci, nadir_eci)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    errors_deg = np.degrees(np.arccos(cos_theta))

    plt.plot(t, errors_deg, color='purple', lw=1, label='Nadir Error')
    plt.axhline(y=np.mean(errors_deg), color='r', linestyle='--', label=f'Mean: {np.mean(errors_deg):.2f}°')

    plt.title("Nadir Pointing Accuracy")
    plt.ylabel("Error [deg]")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.grid(True, alpha=0.3)

    if save:
        filename = f"{prefix}_nadir_error.png"
        plt.savefig(filename, dpi=dpi)
        plt.close()
        print(f"Saved: {filename} at {dpi} DPI")

def plot_dynamics_context(t, df, save=False, prefix="output", dpi=150):
    """Figure 6: Relative Velocity and Gravity Acceleration."""
    fig, axs = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    plt.suptitle("Orbital Dynamics Context", fontsize=14, fontweight='bold')

    # 1. Relative Velocity (v_rel)
    v_rel_norm = np.linalg.norm(df[['v_rel_x', 'v_rel_y', 'v_rel_z']], axis=1)
    axs[0].plot(t, df[['v_rel_x', 'v_rel_y', 'v_rel_z']], lw=1, alpha=0.7)
    axs[0].plot(t, v_rel_norm, 'k--', lw=1.2, label='Total Magnitude')
    axs[0].set_title("Velocity Relative to Atmosphere [m/s]")
    axs[0].legend(['vx_rel', 'vy_rel', 'vz_rel', 'Norm'], loc='upper right', ncol=4)

    # 2. Gravity Acceleration (grav)
    axs[1].plot(t, df[['grav_x', 'grav_y', 'grav_z']], lw=1)
    axs[1].set_title("Gravity Acceleration Vector [m/s²]")
    axs[1].set_ylabel("Accel [m/s²]")

    for ax in axs: ax.grid(True, alpha=0.3)
    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save:
        filename = f"{prefix}_dynamics.png"
        plt.savefig(filename, dpi=dpi)
        plt.close()
        print(f"Saved: {filename} at {dpi} DPI")
    else:
        maximize_window()

def main():
    parser = argparse.ArgumentParser(description="Satellite Telemetry Plotter")
    parser.add_argument("csv_file", nargs="?", default="output.csv", help="Path to CSV")
    parser.add_argument("--save", action="store_true", help="Export to PNG")
    parser.add_argument("--dpi", type=int, default=150, help="Resolution for export (default: 150)")
    args = parser.parse_args()

    # Get the filename without extension
    base_name = Path(args.csv_file).stem

    try:
        df = pd.read_csv(args.csv_file)
        df.columns = df.columns.str.strip()
        t = df['time']

        # Pass the base_name and dpi to each function
        plot_torques(t, df, save=args.save, prefix=base_name, dpi=args.dpi)
        plot_attitude_state(t, df, save=args.save, prefix=base_name, dpi=args.dpi)
        plot_environment(t, df, save=args.save, prefix=base_name, dpi=args.dpi)
        plot_3d_orbit(df, save=args.save, prefix=base_name, dpi=args.dpi)
        plot_nadir_pointing_error(t, df, save=args.save, prefix=base_name, dpi=args.dpi)
        plot_dynamics_context(t, df, save=args.save, prefix=base_name, dpi=args.dpi)

        if not args.save:
            plt.show()
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
