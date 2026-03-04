import os
import sys
import toml
import copy
import random
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import pytz

def randomize_value(base_value, variance_ratio):
    """Applies a uniform random variance to a base value."""
    min_val = base_value * (1.0 - variance_ratio)
    max_val = base_value * (1.0 + variance_ratio)
    return random.uniform(min_val, max_val)

def apply_monte_carlo_variance(config, variance):
    """Mutates the loaded config with randomized values for essential variables."""
    sat = config.get("satellite", {})
    mc_log = {}

    # 1. Randomize Mass
    if "mass" in sat:
        new_mass = round(randomize_value(sat["mass"], variance), 4)
        sat["mass"] = new_mass
        mc_log["mass"] = new_mass

    # 2. Randomize Magnet Properties (Remanence and Dimensions)
    if "magnet" in sat:
        mag = sat["magnet"]
        if "remanence" in mag:
            new_rem = round(randomize_value(mag["remanence"], variance), 4)
            mag["remanence"] = new_rem
            mc_log["magnet_remanence"] = new_rem

        if "cylindrical" in mag:
            cyl = mag["cylindrical"]
            if "length" in cyl:
                cyl["length"] = round(randomize_value(cyl["length"], variance), 5)
            if "radius" in cyl:
                cyl["radius"] = round(randomize_value(cyl["radius"], variance), 6)
        elif "rectangular" in mag:
            rect = mag["rectangular"]
            if "length" in rect:
                rect["length"] = round(randomize_value(rect["length"], variance), 5)
            if "width" in rect:
                rect["width"] = round(randomize_value(rect["width"], variance), 5)
            if "height" in rect:
                rect["height"] = round(randomize_value(rect["height"], variance), 5)

    # 3. Randomize Hysteresis Rod Volumes
    if "rods" in sat:
        mc_log["rods_volume"] = []
        for rod in sat["rods"]:
            # Handle both 'volume' and 'volume_m3' depending on which is used in base TOML
            vol_key = "volume" if "volume" in rod else "volume_m3" if "volume_m3" in rod else None
            if vol_key:
                new_vol = round(randomize_value(rod[vol_key], variance), 11)
                rod[vol_key] = new_vol
                mc_log["rods_volume"].append(new_vol)

    return mc_log

def generate_mc_variants(args, t_end):
    """Generates randomized TOML files based on the input TOML."""
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Load the base TOML
    try:
        with open(args.input_toml, "r") as f:
            base_config = toml.load(f)
    except Exception as e:
        print(f"Error reading base TOML file '{args.input_toml}': {e}")
        sys.exit(1)

    files_to_run =[]

    print(f"Loaded base config: '{args.input_toml}'")
    print(f"Applying +/- {args.variance * 100:.1f}% variance to essential parameters.")

    for i in range(1, args.runs + 1):
        filename = f"{i:04d}.toml"
        filepath = os.path.join(args.output_dir, filename)
        csv_path = filepath.replace(".toml", ".csv")

        # Skip if the corresponding CSV already exists
        if os.path.exists(csv_path):
            print(f"Skipping {filename}: Corresponding CSV '{os.path.basename(csv_path)}' already exists.")
            continue

        # Create a fresh copy of the base config for this MC run
        config = copy.deepcopy(base_config)

        # Override simulation duration and checkpoints
        config["t_end"] = t_end
        config["checkpoint_interval"] = args.checkpoint

        # Apply random spread
        mc_log = apply_monte_carlo_variance(config, args.variance)

        # Append to or create history chain
        if "_mc_history" not in config:
            config["_mc_history"] =[]

        chain_entry = {
            "type": "monte_carlo",
            "base_file": args.input_toml,
            "variance_applied": args.variance,
            "generated_values": mc_log,
            "timestamp": datetime.now(pytz.utc).isoformat()
        }
        config["_mc_history"].append(chain_entry)

        # Save configuration to disk
        with open(filepath, "w") as f:
            toml.dump(config, f)

        files_to_run.append(filepath)

    return files_to_run

def run_simulation(toml_path):
    """Executes the simulator for a single TOML file."""
    csv_path = toml_path.replace(".toml", ".csv")
    cmd = ["./build/pmaos_run", "-o", csv_path, toml_path]

    try:
        # We capture output so threads don't scramble terminal text
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            return f"[SUCCESS] {toml_path} -> {csv_path}"
        else:
            return f"[FAILED] {toml_path}\nError: {result.stderr}"
    except FileNotFoundError:
        return "[ERROR] './build/pmaos_run' command not found. Ensure it is built."
    except Exception as e:
        return f"[ERROR] {toml_path}: {str(e)}"

def main():
    parser = argparse.ArgumentParser(description="Monte Carlo Generator and Runner for AOS Simulation")

    # Positional Argument (The input config file)
    parser.add_argument("input_toml", type=str, help="Path to the base TOML configuration file")

    # Monte Carlo Settings
    parser.add_argument("-n", "--runs", type=int, default=40, help="Number of Monte Carlo iterations (default: 40)")
    # FIXED: Replaced "10%" with "10%%" so argparse doesn't crash formatting the help text
    parser.add_argument("-v", "--variance", type=float, default=0.1, help="Variance range ratio, e.g., 0.1 means +/- 10%% (default: 0.1)")

    # Execution Settings
    parser.add_argument("-j", "--threads", type=int, default=8, help="Maximum number of parallel jobs (default: 8)")
    parser.add_argument("-o", "--output-dir", type=str, default="analysis", help="Output directory (default: 'analysis')")
    parser.add_argument("-d", "--duration", type=str, choices=["2w", "2y"], default="2w", help="Simulation duration: '2w' or '2y' (default: 2w)")
    parser.add_argument("-c", "--checkpoint", type=float, default=600.0, help="Checkpoint interval in seconds (default: 600.0)")

    args = parser.parse_args()

    # Parse duration to seconds
    if args.duration == "2y":
        t_end = 2 * 365 * 24 * 3600
        print(f"Simulation duration set to 2 Years ({t_end} seconds).")
    else:
        t_end = 2 * 7 * 24 * 3600
        print(f"Simulation duration set to 2 Weeks ({t_end} seconds).")

    print(f"Checking for existing/generating up to {args.runs} variations in '{args.output_dir}' directory...")

    toml_files = generate_mc_variants(args, t_end)

    if not toml_files:
        print("\nAll corresponding CSV files already exist. Nothing new to simulate.")
        return

    print(f"\nStarting {len(toml_files)} new simulations using {args.threads} threads...")

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_file = {executor.submit(run_simulation, path): path for path in toml_files}

        for future in as_completed(future_to_file):
            print(future.result())

if __name__ == "__main__":
    main()
