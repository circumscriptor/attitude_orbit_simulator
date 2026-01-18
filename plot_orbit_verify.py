import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys


def plot_orbit_verify(filename):
    try:
        df = pd.read_csv(filename)
        df.columns = df.columns.str.strip()
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return

    # Helper to check/calculate magnitude
    if "r_mag" not in df.columns and all(
        c in df.columns for c in ["r_x", "r_y", "r_z"]
    ):
        df["r_mag"] = np.sqrt(df["r_x"] ** 2 + df["r_y"] ** 2 + df["r_z"] ** 2)

    if "v_mag" not in df.columns and all(
        c in df.columns for c in ["v_x", "v_y", "v_z"]
    ):
        df["v_mag"] = np.sqrt(df["v_x"] ** 2 + df["v_y"] ** 2 + df["v_z"] ** 2)

    active_plots = []
    if "r_mag" in df.columns:
        active_plots.append("radius")
    if "v_mag" in df.columns:
        active_plots.append("velocity")

    if not active_plots:
        print("No plottable orbital data found (missing r_mag or v_mag).")
        return

    fig, axes = plt.subplots(
        len(active_plots), 1, figsize=(12, 4 * len(active_plots)), sharex=True
    )
    if len(active_plots) == 1:
        axes = [axes]

    curr = 0
    if "radius" in active_plots:
        axes[curr].plot(df["time"], df["r_mag"] / 1000.0, color="black", label="|r|")
        axes[curr].set_ylabel("Distance [km]")
        axes[curr].set_title("Orbital Radius (ECI)")
        axes[curr].grid(True, alpha=0.3)
        curr += 1

    if "velocity" in active_plots:
        axes[curr].plot(df["time"], df["v_mag"], color="blue", label="|v|")
        axes[curr].set_ylabel("Velocity [m/s]")
        axes[curr].set_title("Orbital Velocity Magnitude")
        axes[curr].grid(True, alpha=0.3)

    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_orbit_verify.py <data.csv>")
    else:
        plot_orbit_verify(sys.argv[1])
