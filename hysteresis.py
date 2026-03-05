import os
import re
import subprocess
import tomllib
import tomli_w
import pandas as pd
import numpy as np
import argparse

# Force matplotlib to not use any Xwindows/Qt backend.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def calculate_loop_area(h, b):
    return 0.5 * np.abs(np.dot(h, np.roll(b, 1)) - np.dot(b, np.roll(h, 1)))

def get_last_cycle_indices(h_data):
    """Finds the last full cycle by detecting direction changes in H."""
    # Find where the direction of H changes (slopes)
    slopes = np.diff(h_data)
    # Filter out tiny numerical noise
    slopes[np.abs(slopes) < 1e-6] = 0
    direction_changes = np.where(np.diff(np.sign(slopes)) != 0)[0]

    # We need at least 3 direction changes to form one full cycle (peak-trough-peak)
    if len(direction_changes) < 3:
        return 0, len(h_data)

    # The last full cycle is between the 3rd to last and the last direction change
    return direction_changes[-3], direction_changes[-1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("param_csv", help="CSV with material parameters")
    parser.add_argument("-c", "--config", default="template.toml")
    parser.add_argument("-b", "--bin", default="./build/pmaos_vh")
    parser.add_argument("-o", "--outdir", default="analysis/hysteresis")
    args = parser.parse_args()

    params_df = pd.read_csv(args.param_csv)
    os.makedirs(args.outdir, exist_ok=True)

    ratings = []

    for i, row in params_df.iterrows():
        mat_clean = re.sub(r'[^a-zA-Z0-9]', '_', str(row['material']))
        mat_name = re.sub(r'_+', '_', mat_clean).strip('_')
        run_id = f"{i:02d}_{mat_name}"

        toml_path = os.path.join(args.outdir, f"{run_id}.toml")
        csv_path = os.path.join(args.outdir, f"{run_id}.csv")
        png_path = os.path.join(args.outdir, f"{run_id}.png")

        # 1. Update TOML
        with open(args.config, "rb") as f:
            config = tomllib.load(f)
        config['satellite']['hysteresis'] = row.drop('material').to_dict()
        with open(toml_path, "wb") as f:
            tomli_w.dump(config, f)

        # 2. Execute
        print(f"[{i+1}/{len(params_df)}] Simulating {row['material']}...")
        subprocess.run([args.bin, "-o", csv_path, toml_path], check=True, capture_output=True)

        # 3. Analyze Data
        df = pd.read_csv(csv_path)
        h, b = df['H_Am'].values, df['B_T'].values

        start_idx, end_idx = get_last_cycle_indices(h)
        h_cycle = h[start_idx:end_idx]
        b_cycle = b[start_idx:end_idx]

        # Calculate stability by comparing the last cycle to the one before it
        # (Assuming the cycle before has the same length)
        cycle_len = end_idx - start_idx
        prev_start = start_idx - cycle_len

        if prev_start >= 0:
            b_prev = b[prev_start:start_idx]
            stability_err = np.mean(np.abs(b_cycle - b_prev)) / (np.max(b) - np.min(b) + 1e-9)
            status = "YES" if stability_err < 0.01 else f"Drifting ({stability_err:.1%})"
        else:
            status = "INCOMPLETE (Need more cycles)"

        # Save plot: Plot the WHOLE path in light grey, and the LAST cycle in bold blue
        plt.figure(figsize=(8, 6))
        plt.plot(h, b, color='lightgray', alpha=0.5, label='Full Path (Drift)')
        plt.plot(h_cycle, b_cycle, color='blue', linewidth=2, label='Last Cycle')
        plt.title(f"{row['material']} | Status: {status}")
        plt.xlabel("H [A/m]")
        plt.ylabel("B [T]")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(png_path)
        plt.close()

if __name__ == "__main__":
    main()
