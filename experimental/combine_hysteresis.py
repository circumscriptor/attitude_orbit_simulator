import os
import glob
import argparse
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# The specific materials requested by their expected index in the parameter CSV
TARGET_MATERIALS = {
    17: "Vacodur 50 Mechanical",
    27: "Permalloy / Mumetall",
    29: "Fe50-Ni50",
    32: "Amorphous (Fe-based)",
    33: "Amorphous (Co-based)",
    34: "Nanocrystaline (FINEMENT)"
}

def get_cycle_indices(h_data):
    """Finds indices of full cycles using H-field zero crossings."""
    return np.where((h_data[:-1] < 0) & (h_data[1:] >= 0))[0]

def main():
    parser = argparse.ArgumentParser(description="Plot combined hysteresis loops for selected materials.")
    parser.add_argument("-i", "--indir", default="analysis/hysteresis", help="Directory containing the simulated CSV files")
    parser.add_argument("-o", "--output", default="analysis/hysteresis_materials_bh_loop_comparison.png", help="Output PNG filename")
    args = parser.parse_args()

    if not os.path.exists(args.indir):
        print(f"Error: Directory '{args.indir}' does not exist.")
        return

    plt.figure(figsize=(10, 8))

    # Generate distinct colors for the lines
    colors = plt.cm.tab10(np.linspace(0, 1, len(TARGET_MATERIALS)))

    plotted_count = 0

    for (idx, name), color in zip(TARGET_MATERIALS.items(), colors):
        # Locate the CSV file by matching the prefix (e.g., "17_*.csv")
        search_pattern = os.path.join(args.indir, f"{idx:02d}_*.csv")
        matching_files = glob.glob(search_pattern)

        if not matching_files:
            print(f"Warning: Could not find CSV for index {idx:02d} ({name}) in {args.indir}")
            continue

        csv_path = matching_files[0]
        try:
            df = pd.read_csv(csv_path)
            h = df['H_Am'].values
            b = df['B_T'].values

            # Extract the last full cycle to ensure we plot the stabilized loop
            cycle_starts = get_cycle_indices(h)
            if len(cycle_starts) >= 3:
                s2, s3 = cycle_starts[-2], cycle_starts[-1]
                h_plot, b_plot = h[s2:s3], b[s2:s3]
            else:
                # Fallback to the whole path if not enough cycles were found
                h_plot, b_plot = h, b

            plt.plot(h_plot, b_plot, label=name, color=color, linewidth=1)
            plotted_count += 1
            print(f"Added {name} to plot.")

        except Exception as e:
            print(f"Error reading/plotting {csv_path}: {e}")

    if plotted_count > 0:
        plt.title("Comparison of Hysteresis Loops", fontsize=14, fontweight='bold')
        plt.xlabel("Magnetic Field Strength, H [A/m]", fontsize=12)
        plt.ylabel("Magnetic Flux Density, B [T]", fontsize=12)

        # Add a grid with a bit more detail
        plt.grid(True, which='major', linestyle='-', alpha=0.5)
        plt.grid(True, which='minor', linestyle=':', alpha=0.2)
        plt.minorticks_on()

        plt.axhline(0, color='black', linewidth=1)
        plt.axvline(0, color='black', linewidth=1)

        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()

        plt.savefig(args.output, dpi=300)
        print(f"\nSuccess! Combined plot saved to '{args.output}'")
    else:
        print("\nNo materials were plotted. Check if the input directory contains the required CSV files.")

if __name__ == "__main__":
    main()
