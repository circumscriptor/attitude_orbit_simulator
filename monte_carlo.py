import os
import toml
import math
import copy
import random
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def get_base_config(t_end, checkpoint_interval):
    """Returns the base simulation dictionary structure."""
    return {
        "t_start": 0,
        "t_end": t_end,
        "dt_initial": 0.1,
        "checkpoint_interval": checkpoint_interval,
        "angular_velocity":[0.5, 0.5, 0.5],
        "stepper_function": 1,
        "satellite": {
            "mass": 1.3,
            "hysteresis": {},  # Populated dynamically
            "uniform": {
                "dimensions": [0.1, 0.1, 0.1],
                "drag_coefficient": 2.2,
                "specular_reflection_coefficient": 0.1,
                "diffuse_reflection_coefficient": 0.2
            },
            "magnet": {}, # Populated dynamically
            "rods":[]    # Populated dynamically
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

def generate_mc_variants(args, t_end):
    """Generates randomized TOML files based on Monte Carlo distributions."""
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    files_to_run =[]

    for i in range(1, args.runs + 1):
        filename = f"{i:04d}.toml"
        filepath = os.path.join(args.output_dir, filename)
        csv_path = filepath.replace(".toml", ".csv")

        # Skip if the corresponding CSV already exists
        if os.path.exists(csv_path):
            print(f"Skipping {filename}: Corresponding CSV '{os.path.basename(csv_path)}' already exists.")
            continue

        config = get_base_config(t_end, args.checkpoint)

        # TODO: Make one side a bit longer
        config["satellite"]["uniform"]["dimensions"] = [0.1, 0.1, 0.1]

        # Base mass approx 1.33 kg, ±10% variance
        base_mass = 1.33
        config["satellite"]["mass"] = round(random.uniform(base_mass * 0.9, base_mass * 1.1), 3)

        # 2. Hysteresis Material Parameters (Realistic Soft Magnetic Material e.g., HyMu-80)
        config["satellite"]["hysteresis"] = {
            "ms": round(random.uniform(600000, 800000), 1),
            "a": round(random.uniform(10.0, 30.0), 2),
            "k": round(random.uniform(1.0, 4.0), 2),
            "c": round(random.uniform(0.01, 0.1), 3),
            "alpha": round(random.uniform(1.0e-6, 5.0e-5), 7)
        }

        # 3. Permanent Magnet Properties
        magnet = {
            "remanence": round(random.uniform(1.0, 1.35), 3),
            "relative_permeability": round(random.uniform(1.02, 1.20), 3),
            "orientation": [0.0, 0.0, 1.0]
        }

        if random.choice(["cylindrical", "rectangular"]) == "cylindrical":
            magnet["cylindrical"] = {
                "length": round(random.uniform(0.015, 0.05), 4),
                "radius": round(random.uniform(0.001, 0.003), 4)
            }
        else:
            magnet["rectangular"] = {
                "length": round(random.uniform(0.015, 0.05), 4),
                "width": round(random.uniform(0.002, 0.005), 4),
                "height": round(random.uniform(0.002, 0.005), 4)
            }
        config["satellite"]["magnet"] = magnet

        # 4. Hysteresis Rods (2, 4, or 6)
        num_rods = random.choice([2, 4, 6])
        rods =[]

        for r in range(num_rods):
            if r % 2 == 0:
                orientation =[1.0, 0.0, 0.0]  # X-axis
            else:
                orientation =[0.0, 1.0, 0.0]  # Y-axis

            rods.append({
                "volume": round(random.uniform(5.0e-8, 5.0e-7), 10),
                "orientation": orientation
            })

        config["satellite"]["rods"] = rods

        # Save configuration to disk
        with open(filepath, "w") as f:
            toml.dump(config, f)

        files_to_run.append(filepath)

    return files_to_run

def run_simulation(toml_path):
    """Executes the simulator for a single TOML file."""
    csv_path = toml_path.replace(".toml", ".csv")
    cmd =["./build/pmaos_run", "-o", csv_path, toml_path]

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
    parser.add_argument("-n", "--runs", type=int, default=40, help="Number of Monte Carlo iterations (default: 40)")
    parser.add_argument("-j", "--threads", type=int, default=8, help="Maximum number of parallel jobs (default: 8)")
    parser.add_argument("-o", "--output-dir", type=str, default="analysis", help="Output directory (default: 'analysis')")
    parser.add_argument("-d", "--duration", type=str, choices=["2w", "2y"], default="2w", help="Simulation duration: '2w' for 2 weeks, '2y' for 2 years (default: 2w)")
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
