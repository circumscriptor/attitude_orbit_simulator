import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

# Define shorthand groups for easy CLI use
GROUPS = {
    "pos": ["r_x", "r_y", "r_z", "r_mag"],
    "vel": ["v_x", "v_y", "v_z", "v_mag"],
    "quat": ["q_w", "q_x", "q_y", "q_z"],
    "euler": ["roll_deg", "pitch_deg", "yaw_deg"],
    "omega": ["omega_x", "omega_y", "omega_z"],
    "nadir": ["nadir_error_deg"],
    "rods": ["M_"],  # Special case: matches any column starting with M_
}


def resolve_requested_columns(requested, available):
    """Maps group names or specific column names to actual available columns."""
    to_plot = []

    for item in requested:
        if item in GROUPS:
            # Add all columns in the group that exist in the CSV
            if item == "rods":
                to_plot.extend([c for c in available if c.startswith("M_")])
            else:
                to_plot.extend([c for c in GROUPS[item] if c in available])
        elif item in available:
            to_plot.append(item)
        else:
            print(f"Warning: '{item}' not found in CSV or as a valid group name.")

    return sorted(list(set(to_plot)), key=lambda x: available.index(x))


def group_columns_for_plotting(columns):
    """Groups columns back into logical subplots (e.g., all quats together)."""
    plot_groups = {}
    for col in columns:
        found_group = False
        for g_name, g_cols in GROUPS.items():
            if g_name == "rods" and col.startswith("M_"):
                plot_groups.setdefault("Rod Magnetization", []).append(col)
                found_group = True
                break
            elif col in g_cols:
                plot_groups.setdefault(g_name.upper(), []).append(col)
                found_group = True
                break

        if not found_group:
            plot_groups.setdefault("Other", []).append(col)
    return plot_groups


def main():
    parser = argparse.ArgumentParser(description="AOS Universal Verification Plotter")
    parser.add_argument("filename", help="CSV file to plot")
    parser.add_argument(
        "-c",
        "--cols",
        nargs="+",
        default=["pos", "vel", "quat", "euler", "omega", "nadir", "rods"],
        help="Select groups (pos, vel, quat, euler, omega, nadir, rods) or specific column names.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List all available columns in the CSV and exit.",
    )

    args = parser.parse_args()

    try:
        df = pd.read_csv(args.filename)
        df.columns = df.columns.str.strip()
    except Exception as e:
        print(f"Error: Could not read {args.filename}: {e}")
        sys.exit(1)

    available = list(df.columns)

    if args.list:
        print("Available columns:", ", ".join(available))
        sys.exit(0)

    # 1. Handle automatic magnitude calculation if components exist but mags don't
    if "r_mag" not in available and all(x in available for x in ["r_x", "r_y", "r_z"]):
        df["r_mag"] = np.sqrt(df["r_x"] ** 2 + df["r_y"] ** 2 + df["r_z"] ** 2)
        available.append("r_mag")

    # 2. Resolve what the user wants to see
    selected_cols = resolve_requested_columns(args.cols, available)
    if not selected_cols:
        print("No valid columns selected for plotting.")
        sys.exit(1)

    plot_data = group_columns_for_plotting(selected_cols)

    # 3. Dynamic Plotting
    n_rows = len(plot_data)
    fig, axes = plt.subplots(n_rows, 1, figsize=(12, 3.5 * n_rows), sharex=True)
    if n_rows == 1:
        axes = [axes]

    for ax, (title, cols) in zip(axes, plot_data.items()):
        for c in cols:
            ax.plot(df["time"], df[c], label=c, alpha=0.8)

        ax.set_title(title)
        ax.set_ylabel("Value")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="right", fontsize="small")

    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
