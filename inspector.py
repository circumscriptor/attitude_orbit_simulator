import pandas as pd
import numpy as np
import argparse
import sys

# WGS84 Constants
EARTH_RADIUS_M = 6378137.0


def calculate_magnitude(df, prefix):
    """Calculates magnitude from components if the scalar column doesn't exist."""
    cols = [f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"]

    # If the magnitude column already exists (e.g., 'v' or 'r'), use it
    if prefix in df.columns:
        return df[prefix]

    # Otherwise calculate from components
    if all(c in df.columns for c in cols):
        return np.sqrt(df[cols[0]] ** 2 + df[cols[1]] ** 2 + df[cols[2]] ** 2)

    return None


def analyze_detumbling(df, threshold=0.02):
    """
    Determines if and when the spacecraft stabilized.
    Criteria: Angular velocity magnitude stays below threshold until the end.
    """
    w_mag = calculate_magnitude(df, "w")
    if w_mag is None:
        return "N/A (Missing angular velocity data)"

    # Find indices where omega is ABOVE threshold
    above_threshold = w_mag > threshold

    if not above_threshold.any():
        return "0.00 s (Started stable)"

    last_unstable_idx = above_threshold[::-1].idxmax()

    # If the last unstable point is the very last data point, it never stabilized
    if last_unstable_idx == df.index[-1]:
        return f"Not stabilized (Final rate: {w_mag.iloc[-1]:.4f} rad/s)"

    settling_time = df.loc[last_unstable_idx, "time"]
    return f"{settling_time:.2f} s"


def check_quaternion_norm(df):
    """Checks if quaternions remained normalized (Integration Health)."""
    cols = ["q_w", "q_x", "q_y", "q_z"]
    if not all(c in df.columns for c in cols):
        return "N/A"

    norms = np.sqrt(df["q_w"] ** 2 + df["q_x"] ** 2 + df["q_y"] ** 2 + df["q_z"] ** 2)
    max_error = np.max(np.abs(norms - 1.0))
    return max_error


def inspect(filename, stable_threshold):
    try:
        df = pd.read_csv(filename)
        df.columns = df.columns.str.strip()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)

    t = df["time"]
    dt = t.diff().mean()
    duration = t.iloc[-1] - t.iloc[0]

    print("=" * 60)
    print(f"  SIMULATION REPORT: {filename}")
    print("=" * 60)

    # --- 1. General Stats ---
    print(f"[-] Data Points:      {len(df)}")
    print(f"[-] Duration:         {duration:.2f} s ({duration / 3600:.2f} hours)")
    print(f"[-] Avg Time Step:    {dt:.4f} s")

    # Check for NaNs
    nan_count = df.isna().sum().sum()
    if nan_count > 0:
        print(f"\033[91m[!] DATA CORRUPTION: Found {nan_count} NaN values!\033[0m")
    else:
        print("[-] Data Integrity:   OK (No NaNs)")

    # --- 2. Orbit Analysis ---
    print("-" * 60)
    print("ORBITAL MECHANICS")
    r_mag = calculate_magnitude(df, "r")
    if r_mag is not None:
        altitudes = (r_mag - EARTH_RADIUS_M) / 1000.0  # km
        min_alt = altitudes.min()
        max_alt = altitudes.max()

        print(f"[-] Altitude Range:   {min_alt:.2f} km <-> {max_alt:.2f} km")
        if min_alt < 0:
            print(
                "\033[91m[!] CRITICAL: Spacecraft crashed into Earth (Alt < 0)!\033[0m"
            )
        elif min_alt < 150:
            print(
                "\033[93m[!] WARNING: Altitude dangerously low (< 150 km). Drag will be massive.\033[0m"
            )
    else:
        print("[-] Position data not found.")

    v_mag = calculate_magnitude(df, "v")
    if v_mag is not None:
        print(f"[-] Velocity Avg:     {v_mag.mean():.2f} m/s")

    # --- 3. Attitude & Control ---
    print("-" * 60)
    print("ATTITUDE & CONTROL")

    # Integration Health
    q_err = check_quaternion_norm(df)
    if isinstance(q_err, float):
        status = "OK" if q_err < 1e-4 else "DRIFTING"
        color = "\033[92m" if q_err < 1e-4 else "\033[91m"
        print(f"[-] Quat Norm Error:  {color}{q_err:.2e} ({status})\033[0m")

    # Tumbling
    w_mag = calculate_magnitude(df, "w")
    if w_mag is not None:
        init_tumble = w_mag.iloc[0]
        final_tumble = w_mag.iloc[-1]
        peak_tumble = w_mag.max()

        print(
            f"[-] Initial Tumble:   {np.degrees(init_tumble):.2f} deg/s ({init_tumble:.3f} rad/s)"
        )
        print(f"[-] Peak Tumble:      {np.degrees(peak_tumble):.2f} deg/s")
        print(f"[-] Final Tumble:     {np.degrees(final_tumble):.2f} deg/s")

        # Stabilization Check
        settle_time = analyze_detumbling(df, stable_threshold)
        print(
            f"[-] Settling Time:    \033[1m{settle_time}\033[0m (Threshold: {stable_threshold} rad/s)"
        )
    else:
        print("[-] Angular velocity data not found.")

    # --- 4. Magnetics ---
    mag_cols = [c for c in df.columns if c.startswith("M_")]
    if mag_cols:
        print("-" * 60)
        print(f"MAGNETIC HYSTERESIS ({len(mag_cols)} rods)")
        for col in mag_cols:
            m_max = df[col].max()
            m_min = df[col].min()
            m_abs_avg = df[col].abs().mean()
            print(
                f"[-] {col}: Range [{m_min:.1f}, {m_max:.1f}] A/m | Avg Abs: {m_abs_avg:.1f} A/m"
            )

            if m_max == 0.0 and m_min == 0.0:
                print(
                    f"\033[93m    [!] Warning: Rod {col} is inactive (0.0 magnetization).\033[0m"
                )

    print("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect AOS Simulation Results")
    parser.add_argument(
        "filename", nargs="?", default="output.csv", help="CSV file to inspect"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.02,
        help="Angular velocity threshold for stability [rad/s] (default: 0.02)",
    )

    args = parser.parse_args()
    inspect(args.filename, args.threshold)
