import os
import sys
import toml
import argparse
import subprocess
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import pytz

# Force matplotlib to not use any Xwindows/Qt backend.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# WGS84 Constants
EARTH_RADIUS_M = 6378137.0

def get_baseline_config(t_end, checkpoint_interval):
    """Returns a perfectly stable, baseline 1U CubeSat configuration."""
    return {
        "t_start": 0,
        "t_end": t_end,
        "dt_initial": 0.1,
        "checkpoint_interval": checkpoint_interval,
        "angular_velocity": [0.5, 0.5, 0.5],
        "stepper_function": 1,
        "satellite": {
            "mass": 1.33,
            "uniform": {
                "dimensions":[0.1, 0.1, 0.1],
                "drag_coefficient": 2.2,
                "specular_reflection_coefficient": 0.1,
                "diffuse_reflection_coefficient": 0.2
            },
            "hysteresis": {
                "ms": 750000.0,
                "a": 12.0,
                "k": 1.5,
                "c": 0.02,
                "alpha": 1.0e-5
            },
            "magnet": {
                "remanence": 1.2,
                "relative_permeability": 1.05,
                "orientation": [0.0, 0.0, 1.0],
                "cylindrical": {
                    "length": 0.02,
                    "radius": 0.002
                }
            },
            "rods":[
                {"volume": 2.0e-7, "orientation":[1.0, 0.0, 0.0]},
                {"volume": 2.0e-7, "orientation":[-1.0, 0.0, 0.0]},
                {"volume": 2.0e-7, "orientation": [0.0, 1.0, 0.0]},
                {"volume": 2.0e-7, "orientation":[0.0, -1.0, 0.0]}
            ]
        },
        "orbit": {
            "semi_major_axis": 6818137.0,
            "eccentricity": 0.001,
            "inclination": 1.396263,
            "raan": 0.0,
            "arg_of_periapsis": 0.0,
            "mean_anomaly": 0.0
        },
        "environment": {
            "start_year_decimal": 2026.5,
            "gravity_model_degree": 12,
            "gravity_model_order": 12
        },
        "observer": {
            "exclude_elements": False,
            "exclude_magnitudes": False
        }
    }

def apply_variable_sweep(config, variable, value):
    """Mutates the baseline config with the swept variable."""
    if variable == "mass":
        config["satellite"]["mass"] = float(value)
    elif variable == "rod_volume":
        for rod in config["satellite"]["rods"]:
            rod["volume"] = float(value)
    elif variable == "magnet_remanence":
        config["satellite"]["magnet"]["remanence"] = float(value)
    elif variable == "hysteresis_k":
        config["satellite"]["hysteresis"]["k"] = float(value)
    elif variable == "hysteresis_a":
        config["satellite"]["hysteresis"]["a"] = float(value)
    elif variable == "hysteresis_ms":
        config["satellite"]["hysteresis"]["ms"] = float(value)
    else:
        raise ValueError(f"Unknown sweep variable: {variable}")

def generate_sweep_files(args, t_end):
    """Generates the TOML files for the parameter sweep."""
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    values = np.linspace(args.min, args.max, args.steps)
    files_to_run =[]

    print(f"Sweeping '{args.variable}' from {args.min} to {args.max} in {args.steps} steps.")

    for i, val in enumerate(values):
        filename = f"sweep_{args.variable}_{(i + 1):03d}.toml"
        filepath = os.path.join(args.output_dir, filename)
        csv_path = filepath.replace(".toml", ".csv")

        # Skip generating/running if it already exists
        if os.path.exists(csv_path):
            print(f"Skipping {filename}: Corresponding CSV already exists.")
            continue

        config = get_baseline_config(t_end, args.checkpoint)
        apply_variable_sweep(config, args.variable, val)

        # We inject the current sweep value into a custom metadata block
        # so the analyzer can easily read it back
        config["_sweep_metadata"] = {"variable": args.variable, "value": float(val)}

        with open(filepath, "w") as f:
            toml.dump(config, f)

        files_to_run.append(filepath)

    return files_to_run

def run_simulation(toml_path):
    """Executes the simulator for a single TOML file."""
    csv_path = toml_path.replace(".toml", ".csv")
    cmd =["./build/pmaos_run", "-o", csv_path, toml_path]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return f"[SUCCESS] {toml_path} -> {csv_path}"
        else:
            return f"[FAILED] {toml_path}\nError: {result.stderr}"
    except FileNotFoundError:
        return "[ERROR] './build/pmaos_run' command not found. Ensure it is built."
    except Exception as e:
        return f"[ERROR] {toml_path}: {str(e)}"

# --- ANALYSIS LOGIC ---

def calculate_magnitude(df, prefix):
    cols = [f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"]
    if prefix in df.columns:
        return df[prefix]
    if all(c in df.columns for c in cols):
        return np.sqrt(df[cols[0]]**2 + df[cols[1]]**2 + df[cols[2]]**2)
    return None

def get_settling_time(df, threshold=0.02):
    w_mag = calculate_magnitude(df, "w")
    if w_mag is None:
        return np.nan
    above_threshold = w_mag > threshold
    if not above_threshold.any():
        return 0.0
    last_unstable_idx = above_threshold[::-1].idxmax()
    if last_unstable_idx == df.index[-1]:
        return np.nan
    return df.loc[last_unstable_idx, "time"]

def analyze_and_plot(args):
    """Reads all completed sweep files and generates trend plots."""
    print("\nAnalyzing sweep results and generating plots...")

    csv_files =[f for f in os.listdir(args.output_dir) if f.startswith(f"sweep_{args.variable}_") and f.endswith(".csv")]
    if not csv_files:
        print("No CSV files found for analysis.")
        return

    results =[]

    for csv_file in sorted(csv_files):
        csv_path = os.path.join(args.output_dir, csv_file)
        toml_path = csv_path.replace(".csv", ".toml")

        try:
            with open(toml_path, "r") as f:
                config = toml.load(f)
                sweep_val = config.get("_sweep_metadata", {}).get("value", np.nan)
        except:
            continue

        try:
            df = pd.read_csv(csv_path)
            df.columns = df.columns.str.strip()

            w_mag = calculate_magnitude(df, "w")
            final_tumble = w_mag.iloc[-1] if w_mag is not None else np.nan
            settle_time = get_settling_time(df, args.threshold)

            results.append({
                "Sweep Value": sweep_val,
                "Settling Time (Days)": settle_time / 86400.0 if not np.isnan(settle_time) else np.nan,
                "Final Tumble (rad/s)": final_tumble
            })
        except Exception as e:
            print(f"Error parsing {csv_file}: {e}")

    df_res = pd.DataFrame(results).sort_values(by="Sweep Value")

    if df_res.empty:
        print("No valid data to plot.")
        return

    # Metadata for plotting
    tz = pytz.timezone('Europe/Bratislava')
    current_time = datetime.now(tz).strftime("%Y-%m-%d %H:%M %Z")

    # --- PLOT 1: Swept Variable vs. Settling Time ---
    plt.figure(figsize=(9, 5))

    # We plot stabilized points as solid dots, and unstabilized (NaN) as 'X' at the top of the graph
    stable = df_res.dropna(subset=["Settling Time (Days)"])
    unstable = df_res[df_res["Settling Time (Days)"].isna()]

    plt.plot(stable["Sweep Value"], stable["Settling Time (Days)"], 'b-o', label='Stabilized')

    if not unstable.empty:
        # Plot unstable points at a high arbitrary Y value just to show they failed
        max_y = stable["Settling Time (Days)"].max() if not stable.empty else 14
        plt.scatter(unstable["Sweep Value"], [max_y * 1.1] * len(unstable), color='red', marker='x', s=100, label='Failed to Stabilize')

    plt.title(f"Sensitivity Analysis: Settling Time vs {args.variable}\nGenerated: {current_time}")
    plt.xlabel(f"{args.variable} Value")
    plt.ylabel("Settling Time (Days)")
    plt.grid(True, ls="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, f"plot_{args.variable}_settling.png"), dpi=150)
    plt.close()

    # --- PLOT 2: Swept Variable vs. Final Tumbling Rate ---
    plt.figure(figsize=(9, 5))
    plt.plot(df_res["Sweep Value"], df_res["Final Tumble (rad/s)"], 'g-s')
    plt.axhline(args.threshold, color='red', linestyle='--', label=f'Threshold ({args.threshold} rad/s)')

    plt.title(f"Sensitivity Analysis: Final Tumble Rate vs {args.variable}")
    plt.xlabel(f"{args.variable} Value")
    plt.ylabel("Final Angular Velocity (rad/s)")
    plt.yscale("log")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, f"plot_{args.variable}_final_tumble.png"), dpi=150)
    plt.close()

    print(f"Plots saved to {args.output_dir}/")

def main():
    parser = argparse.ArgumentParser(description="AOS Parameter Sweep Sensitivity Analysis")

    # Sweep Configuration
    parser.add_argument("-v", "--variable", required=True,
                        choices=["mass", "rod_volume", "magnet_remanence", "hysteresis_k", "hysteresis_a", "hysteresis_ms"],
                        help="The variable to sweep")
    parser.add_argument("--min", type=float, required=True, help="Minimum value of the sweep")
    parser.add_argument("--max", type=float, required=True, help="Maximum value of the sweep")
    parser.add_argument("--steps", type=int, default=10, help="Number of intermediate steps (default: 10)")

    # Simulation Execution Configuration
    parser.add_argument("-j", "--threads", type=int, default=8, help="Max parallel jobs (default: 8)")
    parser.add_argument("-o", "--output-dir", type=str, default="sweep_results", help="Output directory (default: 'sweep_results')")
    parser.add_argument("-d", "--duration", type=str, choices=["2w", "2y"], default="2w", help="Simulation duration (default: 2w)")
    parser.add_argument("-c", "--checkpoint", type=float, default=600.0, help="Checkpoint interval in sec (default: 600.0)")
    parser.add_argument("-t", "--threshold", type=float, default=0.02, help="Stability threshold [rad/s] (default: 0.02)")

    args = parser.parse_args()

    # Parse duration to seconds
    t_end = 2 * 365 * 24 * 3600 if args.duration == "2y" else 2 * 7 * 24 * 3600

    # 1. Generate Configuration Files
    toml_files = generate_sweep_files(args, t_end)

    # 2. Run Simulations
    if toml_files:
        print(f"\nStarting {len(toml_files)} simulations...")
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            future_to_file = {executor.submit(run_simulation, path): path for path in toml_files}
            for future in as_completed(future_to_file):
                print(future.result())
    else:
        print("\nAll simulations already complete. Proceeding to analysis.")

    # 3. Analyze and Plot
    analyze_and_plot(args)

    # TODO: Plot using plotly each angular velocity

if __name__ == "__main__":
    main()
