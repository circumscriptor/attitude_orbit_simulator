import os
import re
import subprocess
import tomllib
import tomli_w
import pandas as pd
import numpy as np
import argparse
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def calculate_loop_area(h, b):
    return 0.5 * np.abs(np.dot(h, np.roll(b, 1)) - np.dot(b, np.roll(h, 1)))

def get_cycle_indices(h_data):
    """Finds indices of full cycles using H-field zero crossings."""
    # Detect where H crosses 0 from negative to positive
    # This marks the start of a new sine period
    indices = np.where((h_data[:-1] < 0) & (h_data[1:] >= 0))[0]
    return indices

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("param_csv", help="CSV with material parameters")
    parser.add_argument("-c", "--config", default="template.toml")
    parser.add_argument("-b", "--bin", default="./build/pmaos_vh")
    parser.add_argument("-o", "--outdir", default="analysis/hysteresis")
    args = parser.parse_args()

    params_df = pd.read_csv(args.param_csv)
    os.makedirs(args.outdir, exist_ok=True)

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

        # 3. Analyze Adaptive Data
        df = pd.read_csv(csv_path)
        h, b = df['H_Am'].values, df['B_T'].values

        cycle_starts = get_cycle_indices(h)
        status = ""
        h_last, b_last = h, b # Fallback

        if len(cycle_starts) >= 3:
            # Indices for the last two full cycles
            s1, s2, s3 = cycle_starts[-3], cycle_starts[-2], cycle_starts[-1]

            h_last, b_last = h[s2:s3], b[s2:s3]
            h_prev, b_prev = h[s1:s2], b[s1:s2]

            # Use interpolation to compare cycles with different step counts (adaptive dt)
            # Normalize cycle progress from 0.0 to 1.0
            interp_last = interp1d(np.linspace(0, 1, len(b_last)), b_last)
            interp_prev = interp1d(np.linspace(0, 1, len(b_prev)), b_prev)

            common_x = np.linspace(0, 1, 1000)
            diff = np.abs(interp_last(common_x) - interp_prev(common_x))

            range_b = np.max(b) - np.min(b) + 1e-9
            stability_err = np.mean(diff) / range_b

            if stability_err > 0.01:
                status = f"\nDrifting ({stability_err:.1%})"
        else:
            status = "\nINCOMPLETE (Need >2 cycles)"

        # 4. Plotting
        plt.figure(figsize=(8, 6))
        plt.plot(h, b, color='lightgray', alpha=0.5, label='Full Path')
        plt.plot(h_last, b_last, color='blue', linewidth=2, label='Last Cycle')
        plt.title(f"{row['material']} {status}")
        plt.xlabel("H [A/m]")
        plt.ylabel("B [T]")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(png_path)
        plt.close()

if __name__ == "__main__":
    main()
