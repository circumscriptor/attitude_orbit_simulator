import os
import sys
import toml
import copy
import argparse
import subprocess
from pathlib import Path
import platform

# Force matplotlib to not use any Xwindows/Qt backend (Headless mode for saving)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Standardized colors for ECI/Body components
C_X, C_Y, C_Z, C_NORM = '#d62728', '#2ca02c', '#1f77b4', '#000000'

# --- SIMULATION MANAGEMENT ---

def prepare_and_run_simulation(args, t_end):
    """Prepares the TOML configuration and runs the verification simulator."""
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    base_name = Path(args.input_toml).stem
    csv_filename = f"{base_name}_vs.csv"
    toml_filename = f"{base_name}_vs.toml"

    csv_path = os.path.join(args.output_dir, csv_filename)
    run_toml_path = os.path.join(args.output_dir, toml_filename)

    # 1. Check if output already exists
    if os.path.exists(csv_path) and not args.force:
        print(f"Skipping simulation: Output CSV already exists at '{csv_path}'.")
        print("Use --force to override and re-run.")
        return csv_path

    # 2. Load the base TOML
    try:
        with open(args.input_toml, "r") as f:
            base_config = toml.load(f)
    except Exception as e:
        print(f"Error reading base TOML file '{args.input_toml}': {e}")
        sys.exit(1)

    # 3. Create a modified config for this run
    config = copy.deepcopy(base_config)
    config["t_end"] = t_end
    config["checkpoint_interval"] = args.checkpoint

    # Save the modified TOML to the output directory
    with open(run_toml_path, "w") as f:
        toml.dump(config, f)

    # 4. Run the simulation
    cmd = ["./build/pmaos_vs", "-o", csv_path, run_toml_path]
    print(f"Running verification simulation: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"[SUCCESS] {run_toml_path} -> {csv_path}")
            return csv_path
        else:
            print(f"[FAILED] Error running simulator:\n{result.stderr}")
            sys.exit(1)
    except FileNotFoundError:
        print("[ERROR] './build/pmaos_vs' command not found. Ensure it is built.")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Simulation failed: {str(e)}")
        sys.exit(1)


# --- VISUALIZATION UTILITIES ---

def quat_to_euler(qw, qx, qy, qz):
    """Convert quaternions to Euler angles (Yaw-Pitch-Roll / Z-Y-X sequence)."""
    sinr_cosp = 2.0 * (qw * qx + qy * qz)
    cosr_cosp = 1.0 - 2.0 * (qx**2 + qy**2)
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    sinp = 2.0 * (qw * qy - qz * qx)
    sinp = np.clip(sinp, -1.0, 1.0)
    pitch = np.arcsin(sinp)

    siny_cosp = 2.0 * (qw * qz + qx * qy)
    cosy_cosp = 1.0 - 2.0 * (qy**2 + qz**2)
    yaw = np.arctan2(siny_cosp, cosy_cosp)

    return np.degrees(roll), np.degrees(pitch), np.degrees(yaw)

# --- PLOTTING FUNCTIONS ---

def plot_attitude_kinematics(t, df, prefix="output", dpi=150):
    fig, axs = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("Attitude Kinematics & Rotation Rates", fontsize=14, fontweight='bold')

    roll, pitch, yaw = quat_to_euler(df['q_w'], df['q_x'], df['q_y'], df['q_z'])
    axs[0].plot(t, roll, color=C_X, lw=1, label='Roll (X)')
    axs[0].plot(t, pitch, color=C_Y, lw=1, label='Pitch (Y)')
    axs[0].plot(t, yaw, color=C_Z, lw=1, label='Yaw (Z)')
    axs[0].set_title("Euler Angles [deg]")
    axs[0].set_ylabel("Angle [deg]")
    axs[0].legend(loc='upper right', ncol=3)

    axs[1].plot(t, df['q_w'], color='black', lw=1, label='qw')
    axs[1].plot(t, df['q_x'], color=C_X, lw=1, label='qx')
    axs[1].plot(t, df['q_y'], color=C_Y, lw=1, label='qy')
    axs[1].plot(t, df['q_z'], color=C_Z, lw=1, label='qz')
    axs[1].set_title("Attitude Quaternions (ECI to Body)")
    axs[1].set_ylabel("Value")
    axs[1].legend(loc='upper right', ncol=4)

    w_norm = np.linalg.norm(df[['w_x', 'w_y', 'w_z']], axis=1)
    axs[2].plot(t, df['w_x'], color=C_X, lw=1, label='ωx')
    axs[2].plot(t, df['w_y'], color=C_Y, lw=1, label='ωy')
    axs[2].plot(t, df['w_z'], color=C_Z, lw=1, label='ωz')
    axs[2].plot(t, w_norm, color=C_NORM, lw=1.2, linestyle='--', label='|ω| Total')
    axs[2].set_title("Body Angular Velocity [rad/s]")
    axs[2].set_ylabel("Rate[rad/s]")
    axs[2].legend(loc='upper right', ncol=4)

    for ax in axs: ax.grid(True, alpha=0.3)
    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(f"{prefix}_01_attitude.png", dpi=dpi)
    plt.close()

def plot_adcs_pointing(t, df, prefix="output", dpi=150):
    fig, axs = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("ADCS Pointing & Alignment Errors", fontsize=14, fontweight='bold')

    r_vec = df[['r_x', 'r_y', 'r_z']].values
    nadir_eci = -r_vec / np.linalg.norm(r_vec, axis=1)[:, np.newaxis]
    qw, qx, qy, qz = df['q_w'], df['q_x'], df['q_y'], df['q_z']
    z_body_eci = np.stack([
        2*(qx*qz + qw*qy), 2*(qy*qz - qw*qx), qw**2 - qx**2 - qy**2 + qz**2
    ], axis=1)
    cos_theta_nadir = np.einsum('ij,ij->i', z_body_eci, nadir_eci)
    nadir_err = np.degrees(np.arccos(np.clip(cos_theta_nadir, -1.0, 1.0)))

    axs[0].plot(t, nadir_err, color='purple', lw=1, label='Angle to Nadir')
    axs[0].axhline(y=np.mean(nadir_err), color='r', linestyle='--', label=f'Mean: {np.mean(nadir_err):.1f}°')
    axs[0].set_title("Nadir Alignment (Body Z to Earth Center)")
    axs[0].set_ylabel("Error [deg]")
    axs[0].legend(loc='upper right')

    mag_vec = df[['mag_x', 'mag_y', 'mag_z']].values
    mag_norm = np.linalg.norm(mag_vec, axis=1)
    mag_norm[mag_norm == 0] = 1
    mag_dir_body = mag_vec / mag_norm[:, np.newaxis]
    mag_err = np.degrees(np.arccos(np.clip(mag_dir_body[:, 2], -1.0, 1.0)))

    axs[1].plot(t, mag_err, color='teal', lw=1, label='Angle to Local B-Field')
    axs[1].axhline(y=np.mean(mag_err), color='r', linestyle='--', label=f'Mean: {np.mean(mag_err):.1f}°')
    axs[1].set_title("Magnetic Alignment (Body Z to Local B-Field)")
    axs[1].set_ylabel("Error [deg]")
    axs[1].legend(loc='upper right')

    if all(col in df.columns for col in['sun_x', 'sun_y', 'sun_z']):
        sun_vec = df[['sun_x', 'sun_y', 'sun_z']].values
        sun_norm = np.linalg.norm(sun_vec, axis=1)
        sun_norm[sun_norm == 0] = 1
        sun_dir = sun_vec / sun_norm[:, np.newaxis]
        sun_err = np.degrees(np.arccos(np.clip(sun_dir[:, 2], -1.0, 1.0)))
        axs[2].plot(t, sun_err, color='orange', lw=1, label='Angle to Sun')
        axs[2].set_title("Sun Alignment (Body Z to Sun Vector)")
    else:
        axs[2].text(0.5, 0.5, "Sun Vector Data Unavailable", ha='center', va='center')

    axs[2].set_ylabel("Angle [deg]")
    axs[2].legend(loc='upper right')

    for ax in axs: ax.grid(True, alpha=0.3)
    plt.xlabel("Time[s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(f"{prefix}_02_pointing.png", dpi=dpi)
    plt.close()

def plot_torque_envelopes(t, df, prefix="output", dpi=150):
    fig, axs = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("Disturbance & Control Torque Envelopes", fontsize=14, fontweight='bold')

    torques = {
        'Permanent Mag':['t_mag_x', 't_mag_y', 't_mag_z'],
        'Hysteresis Rods':['t_rods_x', 't_rods_y', 't_rods_z'],
        'Gravity Grad':['t_grav_x', 't_grav_y', 't_grav_z'],
        'Gyroscopic':['t_gyro_x', 't_gyro_y', 't_gyro_z'],
        'Aero/SRP':['t_face_x', 't_face_y', 't_face_z']
    }

    for name, cols in torques.items():
        if all(c in df.columns for c in cols):
            mag = np.linalg.norm(df[cols], axis=1)
            axs[0].semilogy(t, np.clip(mag, 1e-13, None), lw=1, label=name)

    axs[0].set_title("Torque Magnitude Comparison (Log Scale)")
    axs[0].set_ylabel("Torque [Nm]")
    axs[0].legend(loc='upper right', ncol=3)
    axs[0].grid(True, which='both', alpha=0.3)

    total_t = np.zeros((len(df), 3))
    for cols in torques.values():
        if all(c in df.columns for c in cols):
            total_t += df[cols].values

    axs[1].plot(t, total_t[:, 0], color=C_X, lw=1, label='Total Tx')
    axs[1].plot(t, total_t[:, 1], color=C_Y, lw=1, label='Total Ty')
    axs[1].plot(t, total_t[:, 2], color=C_Z, lw=1, label='Total Tz')
    axs[1].set_title("Total Superimposed Torque (Linear Scale)")
    axs[1].set_ylabel("Torque[Nm]")
    axs[1].legend(loc='upper right', ncol=3)
    axs[1].grid(True, alpha=0.3)

    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(f"{prefix}_03_torques.png", dpi=dpi)
    plt.close()

def plot_passive_magnetic_system(t, df, prefix="output", dpi=150):
    fig, axs = plt.subplots(4, 1, figsize=(14, 12), sharex=True)
    plt.suptitle("Passive Magnetic System Analysis", fontsize=14, fontweight='bold')

    mag_norm = np.linalg.norm(df[['mag_x', 'mag_y', 'mag_z']], axis=1)
    axs[0].plot(t, df['mag_x'], color=C_X, lw=1, label='Bx')
    axs[0].plot(t, df['mag_y'], color=C_Y, lw=1, label='By')
    axs[0].plot(t, df['mag_z'], color=C_Z, lw=1, label='Bz')
    axs[0].plot(t, mag_norm, color=C_NORM, lw=1.2, linestyle='--', label='|B| Total')
    axs[0].set_title("Local Magnetic Field (Body Frame) [Tesla]")
    axs[0].legend(loc='upper right', ncol=4)

    axs[1].plot(t, df['mag_dot_x'], color=C_X, lw=1, label='B-dot x')
    axs[1].plot(t, df['mag_dot_y'], color=C_Y, lw=1, label='B-dot y')
    axs[1].plot(t, df['mag_dot_z'], color=C_Z, lw=1, label='B-dot z')
    axs[1].set_title("Magnetic Field Derivative (B-Dot)[Tesla/s]")
    axs[1].legend(loc='upper right', ncol=3)

    rod_cols =[c for c in df.columns if c.startswith('M_')]
    if rod_cols:
        axs[2].plot(t, df[rod_cols], lw=1)
        axs[2].set_title("Hysteresis Rod Magnetizations [Am²]")
        axs[2].legend(rod_cols, loc='upper right', ncol=min(len(rod_cols), 4))
    else:
        axs[2].text(0.5, 0.5, "No Rod Magnetization Data", ha='center', va='center')

    axs[3].plot(t, df['t_rods_x'], color=C_X, lw=1, label='Rod Tx')
    axs[3].plot(t, df['t_rods_y'], color=C_Y, lw=1, label='Rod Ty')
    axs[3].plot(t, df['t_rods_z'], color=C_Z, lw=1, label='Rod Tz')
    axs[3].set_title("Hysteresis Damping Torque [Nm]")
    axs[3].legend(loc='upper right', ncol=3)

    for ax in axs: ax.grid(True, alpha=0.3)
    plt.xlabel("Time[s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(f"{prefix}_04_passive_mag.png", dpi=dpi)
    plt.close()

def plot_orbit_and_environment(t, df, prefix="output", dpi=150):
    fig, axs = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("External Environment & Aerodynamic Forces", fontsize=14, fontweight='bold')

    axs[0].fill_between(t, 0, df['shadow'], alpha=0.3, color='orange', label='Illuminated Fraction (Shadow)')
    axs[0].set_ylabel("Illumination", color='darkorange')
    axs[0].tick_params(axis='y', labelcolor='darkorange')
    axs[0].set_ylim([-0.1, 1.1])

    ax_srp = axs[0].twinx()
    ax_srp.plot(t, df['solar_p'], color='red', lw=1.5, label='Solar Radiation Pressure')
    ax_srp.set_ylabel("SRP [Pa]", color='red')
    ax_srp.tick_params(axis='y', labelcolor='red')
    axs[0].set_title("Illumination & Solar Pressure")

    v_rel_norm = np.linalg.norm(df[['v_rel_x', 'v_rel_y', 'v_rel_z']], axis=1)
    axs[1].semilogy(t, df['rho'], color='green', lw=1.5, label='Density (ρ)')
    axs[1].set_ylabel("Density [kg/m³]", color='green')
    axs[1].tick_params(axis='y', labelcolor='green')

    ax_vrel = axs[1].twinx()
    ax_vrel.plot(t, v_rel_norm, color='purple', lw=1.5, label='Relative Velocity Magnitude')
    ax_vrel.set_ylabel("Velocity [m/s]", color='purple')
    ax_vrel.tick_params(axis='y', labelcolor='purple')
    axs[1].set_title("Atmospheric Density & Relative Velocity")

    if all(c in df.columns for c in['f_face_x', 'f_face_y', 'f_face_z']):
        f_norm = np.linalg.norm(df[['f_face_x', 'f_face_y', 'f_face_z']], axis=1)
        axs[2].plot(t, df['f_face_x'], color=C_X, lw=1, label='Fx')
        axs[2].plot(t, df['f_face_y'], color=C_Y, lw=1, label='Fy')
        axs[2].plot(t, df['f_face_z'], color=C_Z, lw=1, label='Fz')
        axs[2].plot(t, f_norm, color=C_NORM, lw=1.2, linestyle='--', label='|F| Total')
        axs[2].set_title("Total Aerodynamic / Surface Forces (Body Frame) [N]")
        axs[2].legend(loc='upper right', ncol=4)
        axs[2].set_ylabel("Force [N]")

    for ax in axs: ax.grid(True, alpha=0.3)
    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(f"{prefix}_05_environment.png", dpi=dpi)
    plt.close()

def plot_orbit_3d(df, prefix="output", dpi=150):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    plt.suptitle("3D Orbital Path (Color-coded by Sunlight)", fontsize=14, fontweight='bold')

    rx, ry, rz = df['r_x'].values, df['r_y'].values, df['r_z'].values

    u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
    earth_r = 6371000
    x = earth_r * np.cos(u) * np.sin(v)
    y = earth_r * np.sin(u) * np.sin(v)
    z = earth_r * np.cos(v)
    ax.plot_wireframe(x, y, z, color="green", alpha=0.08)

    sc = ax.scatter(rx, ry, rz, c=df['shadow'], cmap='Wistia', s=2, alpha=0.8)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.02, pad=0.1)
    cbar.set_label('Illumination Fraction (1=Sun, 0=Eclipse)')

    interval = max(1, len(df) // 25)
    ax.quiver(rx[::interval], ry[::interval], rz[::interval],
              -rx[::interval], -ry[::interval], -rz[::interval],
              color='red', length=0.8, arrow_length_ratio=0.1,
              alpha=0.6, label='Nadir Direction')

    ax.set_xlabel('ECI X [m]')
    ax.set_ylabel('ECI Y [m]')
    ax.set_zlabel('ECI Z[m]')

    ax.set_aspect('equal')
    ax.legend(loc='upper left')
    plt.tight_layout()

    plt.savefig(f"{prefix}_06_orbit3d.png", dpi=dpi)
    plt.close()

def plot_per_face_effects(t, df, prefix="output", dpi=150):
    fig, axs = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("Per-Face Surface Effects (Body Frame Magnitudes)", fontsize=14, fontweight='bold')

    if 'd_f0_x' not in df.columns:
        axs[0].text(0.5, 0.5, "Per-face data not found in CSV", ha='center', va='center')
        plt.savefig(f"{prefix}_07_per_face.png", dpi=dpi)
        plt.close()
        return

    colors = plt.cm.tab10(np.linspace(0, 1, 6))

    for i, color in enumerate(colors):
        f_cols =[f'd_f{i}_x', f'd_f{i}_y', f'd_f{i}_z']
        if all(col in df.columns for col in f_cols):
            f_mag = np.linalg.norm(df[f_cols], axis=1)
            axs[0].plot(t, f_mag, color=color, lw=1.2, label=f'Face {i}')

        t_cols =[f's_f{i}_x', f's_f{i}_y', f's_f{i}_z']
        if all(col in df.columns for col in t_cols):
            t_mag = np.linalg.norm(df[t_cols], axis=1)
            axs[1].plot(t, t_mag, color=color, lw=1.2, label=f'Face {i}')

    axs[0].set_title("Aerodynamic Force Magnitude per Face")
    axs[0].set_ylabel("Force [N]")
    axs[0].legend(loc='upper right', ncol=6)
    axs[0].grid(True, alpha=0.3)

    axs[1].set_title("SRP Force Magnitude per Face")
    axs[1].set_ylabel("Force [N]")
    axs[1].legend(loc='upper right', ncol=6)
    axs[1].grid(True, alpha=0.3)

    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(f"{prefix}_07_per_face.png", dpi=dpi)
    plt.close()

def plot_position_velocity(t, df, prefix="output", dpi=150):
    fig, axs = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    plt.suptitle("Orbital Position and Velocity (ECI Frame)", fontsize=14, fontweight='bold')

    rx_km, ry_km, rz_km = df['r_x'] / 1000, df['r_y'] / 1000, df['r_z'] / 1000

    if 'r' in df.columns:
        r_norm_km = df['r'] / 1000
    else:
        r_norm_km = np.linalg.norm(df[['r_x', 'r_y', 'r_z']], axis=1) / 1000

    axs[0].plot(t, rx_km, color=C_X, lw=1.2, label='Rx')
    axs[0].plot(t, ry_km, color=C_Y, lw=1.2, label='Ry')
    axs[0].plot(t, rz_km, color=C_Z, lw=1.2, label='Rz')
    axs[0].plot(t, r_norm_km, color=C_NORM, lw=1.5, linestyle='--', label='|R| Total')
    axs[0].axhline(y=6371, color='green', linestyle=':', lw=1.5, label='Earth Radius (~6371 km)')
    axs[0].set_title("ECI Position [km]")
    axs[0].set_ylabel("Position [km]")
    axs[0].legend(loc='upper right', ncol=5)

    vx_km, vy_km, vz_km = df['v_x'] / 1000, df['v_y'] / 1000, df['v_z'] / 1000

    if 'v' in df.columns:
        v_norm_km = df['v'] / 1000
    else:
        v_norm_km = np.linalg.norm(df[['v_x', 'v_y', 'v_z']], axis=1) / 1000

    axs[1].plot(t, vx_km, color=C_X, lw=1.2, label='Vx')
    axs[1].plot(t, vy_km, color=C_Y, lw=1.2, label='Vy')
    axs[1].plot(t, vz_km, color=C_Z, lw=1.2, label='Vz')
    axs[1].plot(t, v_norm_km, color=C_NORM, lw=1.5, linestyle='--', label='|V| Total')
    axs[1].set_title("ECI Velocity [km/s]")
    axs[1].set_ylabel("Velocity [km/s]")
    axs[1].legend(loc='upper right', ncol=4)

    for ax in axs: ax.grid(True, alpha=0.3)
    plt.xlabel("Time [s]")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(f"{prefix}_08_pos_vel.png", dpi=dpi)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="AOS Verification Simulator & Telemetry Plotter")

    # Input definition
    parser.add_argument("input_toml", type=str, help="Path to the TOML configuration file")

    # Simulation Setup
    parser.add_argument("-o", "--output-dir", type=str, default="analysis/verify", help="Output directory (default: 'analysis/verify')")
    parser.add_argument("-d", "--duration", type=str, choices=["2w", "2y"], default="2w", help="Simulation duration (default: 2w)")
    parser.add_argument("-c", "--checkpoint", type=float, default=600.0, help="Checkpoint interval in sec (default: 600.0)")
    parser.add_argument("--force", action="store_true", help="Force re-run of simulation even if CSV output already exists")

    # Plotting Settings
    parser.add_argument("--dpi", type=int, default=150, help="Resolution for export (default: 150)")

    args = parser.parse_args()

    # Parse duration to seconds
    t_end = 2 * 365 * 24 * 3600 if args.duration == "2y" else 2 * 7 * 24 * 3600

    print(f"\n--- PMAOS Verification Mode ---")

    # 1. Run the simulation automatically
    csv_file = prepare_and_run_simulation(args, t_end)

    # 2. Analyze and Plot
    base_name = Path(csv_file).stem
    save_prefix = os.path.join(args.output_dir, base_name)

    try:
        print(f"Loading telemetry from '{csv_file}' for visualization...")
        df = pd.read_csv(csv_file)
        df.columns = df.columns.str.strip()
        t = df['time']

        plot_attitude_kinematics(t, df, prefix=save_prefix, dpi=args.dpi)
        plot_adcs_pointing(t, df, prefix=save_prefix, dpi=args.dpi)
        plot_torque_envelopes(t, df, prefix=save_prefix, dpi=args.dpi)
        plot_passive_magnetic_system(t, df, prefix=save_prefix, dpi=args.dpi)
        plot_orbit_and_environment(t, df, prefix=save_prefix, dpi=args.dpi)
        plot_orbit_3d(df, prefix=save_prefix, dpi=args.dpi)
        plot_per_face_effects(t, df, prefix=save_prefix, dpi=args.dpi)
        plot_position_velocity(t, df, prefix=save_prefix, dpi=args.dpi)

        print(f"\nAll plots successfully saved to directory: '{args.output_dir}/'")

    except Exception as e:
        print(f"Error during visualization: {e}")

if __name__ == "__main__":
    main()
