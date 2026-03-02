import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def visualize_satellite_data(csv_file):
    # Load data
    df = pd.read_csv(csv_file)
    t = df['time']

    # Create a multi-panel figure
    fig, axs = plt.subplots(4, 2, figsize=(15, 20), constrained_layout=True)
    fig.suptitle(f"Satellite Simulation Analysis: {csv_file}", fontsize=16)

    # 1. Angular Velocity (Body Frame)
    axs[0, 0].plot(t, df['w_x'], label='ωx')
    axs[0, 0].plot(t, df['w_y'], label='ωy')
    axs[0, 0].plot(t, df['w_z'], label='ωz')
    axs[0, 0].set_title("Angular Velocity [rad/s]")
    axs[0, 0].legend()
    axs[0, 0].grid(True)

    # 2. Attitude Quaternions
    axs[0, 1].plot(t, df['q_w'], label='qw')
    axs[0, 1].plot(t, df['q_x'], label='qx')
    axs[0, 1].plot(t, df['q_y'], label='qy')
    axs[0, 1].plot(t, df['q_z'], label='qz')
    axs[0, 1].set_title("Attitude Quaternions")
    axs[0, 1].legend()
    axs[0, 1].grid(True)

    # 3. Magnetic Torque Components
    axs[1, 0].plot(t, df['t_mag_x'], label='τ_mag_x')
    axs[1, 0].plot(t, df['t_mag_y'], label='τ_mag_y')
    axs[1, 0].plot(t, df['t_mag_z'], label='τ_mag_z')
    axs[1, 0].set_title("Magnetic Torque [Nm]")
    axs[1, 0].legend()
    axs[1, 0].grid(True)

    # 4. Gravity Gradient Torque
    axs[1, 1].plot(t, df['t_grav_x'], label='τ_grav_x')
    axs[1, 1].plot(t, df['t_grav_y'], label='τ_grav_y')
    axs[1, 1].plot(t, df['t_grav_z'], label='τ_grav_z')
    axs[1, 1].set_title("Gravity Gradient Torque [Nm]")
    axs[1, 1].legend()
    axs[1, 1].grid(True)

    # 5. Environmental Vectors Magnitudes (Verification)
    mag_norm = np.sqrt(df['mag_x']**2 + df['mag_y']**2 + df['mag_z']**2)
    sun_norm = np.sqrt(df['sun_x']**2 + df['sun_y']**2 + df['sun_z']**2)
    axs[2, 0].plot(t, mag_norm, label='|B| (Tesla)', color='blue')
    axs[2, 0].set_ylabel('Magnetic Field [T]', color='blue')
    ax2_0_twin = axs[2, 0].twinx()
    ax2_0_twin.plot(t, sun_norm, label='|Sun| (m)', color='orange')
    ax2_0_twin.set_ylabel('Sun Distance [m]', color='orange')
    axs[2, 0].set_title("Environmental Field Magnitudes")
    axs[2, 0].grid(True)

    # 6. Shadow Factor & SRP Context
    axs[2, 1].fill_between(t, 0, df['shadow'], alpha=0.3, color='gray', label='Shadow Factor')
    axs[2, 1].set_ylim(-0.1, 1.1)
    axs[2, 1].set_title("Eclipse/Sunlight (1=Sun, 0=Umbra)")
    axs[2, 1].grid(True)

    # 7. Rod Magnetizations (Hysteresis state)
    axs[3, 0].plot(t, df['M_1'], label='Rod 1')
    axs[3, 0].plot(t, df['M_2'], label='Rod 2')
    axs[3, 0].plot(t, df['M_3'], label='Rod 3')
    if 'M_4' in df: axs[3, 0].plot(t, df['M_4'], label='Rod 4')
    axs[3, 0].set_title("Rod Magnetizations [Am^2]")
    axs[3, 0].legend()
    axs[3, 0].grid(True)

    # 8. Orbit Path (ECI Projection)
    axs[3, 1].plot(df['r_x'], df['r_y'])
    axs[3, 1].set_aspect('equal')
    axs[3, 1].set_title("Orbital Path (XY Projection) [m]")
    axs[3, 1].grid(True)

    plt.show()

if __name__ == "__main__":
    # Ensure your filename matches the exported CSV
    visualize_satellite_data("output.csv")
