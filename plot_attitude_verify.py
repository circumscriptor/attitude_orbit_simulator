import pandas as pd
import matplotlib.pyplot as plt
import sys


def plot_attitude_verify(filename):
    try:
        df = pd.read_csv(filename)
        df.columns = df.columns.str.strip()
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return

    # Identify available data groups
    has_euler = all(c in df.columns for c in ["roll_deg", "pitch_deg", "yaw_deg"])
    has_rates = all(c in df.columns for c in ["omega_x", "omega_y", "omega_z"])
    has_nadir = "nadir_error_deg" in df.columns

    active_rows = sum([has_euler, has_rates, has_nadir])
    if active_rows == 0:
        print("No recognizable attitude columns found.")
        return

    fig, axes = plt.subplots(active_rows, 1, figsize=(12, 4 * active_rows), sharex=True)
    if active_rows == 1:
        axes = [axes]

    curr = 0
    # 1. Euler Angles
    if has_euler:
        axes[curr].plot(df["time"], df["roll_deg"], label="Roll (x)", alpha=0.7)
        axes[curr].plot(df["time"], df["pitch_deg"], label="Pitch (y)", linewidth=2)
        axes[curr].plot(df["time"], df["yaw_deg"], label="Yaw (z)", alpha=0.7)
        axes[curr].set_ylabel("Angle [deg]")
        axes[curr].set_title("Euler Angles")
        axes[curr].legend(loc="right")
        axes[curr].grid(True, alpha=0.3)
        curr += 1

    # 2. Body Rates
    if has_rates:
        axes[curr].plot(df["time"], df["omega_x"], label="wx")
        axes[curr].plot(df["time"], df["omega_y"], label="wy")
        axes[curr].plot(df["time"], df["omega_z"], label="wz")
        axes[curr].set_ylabel("Rate [rad/s]")
        axes[curr].set_title("Angular Velocity")
        axes[curr].legend(loc="right")
        axes[curr].grid(True, alpha=0.3)
        curr += 1

    # 3. Nadir Error
    if has_nadir:
        axes[curr].plot(df["time"], df["nadir_error_deg"], color="red", linewidth=1.5)
        axes[curr].set_ylabel("Error [deg]")
        axes[curr].set_title("Nadir Pointing Error (Z-axis to Earth Center)")
        axes[curr].set_ylim(bottom=0)
        axes[curr].grid(True, alpha=0.3)

    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_attitude_verify.py <data.csv>")
    else:
        plot_attitude_verify(sys.argv[1])
