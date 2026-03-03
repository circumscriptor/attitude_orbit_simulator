import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
import toml
from datetime import datetime
import pytz

# Force matplotlib to not use any Xwindows/Qt backend.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# WGS84 Constants
EARTH_RADIUS_M = 6378137.0

def calculate_magnitude(df, prefix):
    cols =[f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"]
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

def process_run(csv_path, toml_path, threshold):
    try:
        df = pd.read_csv(csv_path)
        df.columns = df.columns.str.strip()
    except Exception as e:
        return None, None

    if len(df) < 2 or df.isna().sum().sum() > 0:
        return None, None

    metrics = {"Run": os.path.basename(csv_path).replace(".csv", "")}

    t = df["time"]
    metrics["Duration (Days)"] = (t.iloc[-1] - t.iloc[0]) / 86400.0

    w_mag = calculate_magnitude(df, "w")
    if w_mag is not None:
        metrics["Initial Tumble (rad/s)"] = w_mag.iloc[0]
        metrics["Final Tumble (rad/s)"] = w_mag.iloc[-1]
        metrics["Settling Time (Days)"] = get_settling_time(df, threshold) / 86400.0
    else:
        metrics["Initial Tumble (rad/s)"] = np.nan
        metrics["Final Tumble (rad/s)"] = np.nan
        metrics["Settling Time (Days)"] = np.nan

    r_mag = calculate_magnitude(df, "r")
    if r_mag is not None:
        altitudes = (r_mag - EARTH_RADIUS_M) / 1000.0
        metrics["Min Altitude (km)"] = altitudes.min()
        metrics["Deorbited"] = bool(altitudes.iloc[-1] <= 100.0)
    else:
        metrics["Min Altitude (km)"] = np.nan
        metrics["Deorbited"] = False

    try:
        with open(toml_path, "r") as f:
            config = toml.load(f)
            sat = config.get("satellite", {})
            metrics["Mass (kg)"] = sat.get("mass", np.nan)

            rods = sat.get("rods",[])
            total_vol = sum(r.get("volume", 0.0) for r in rods)
            metrics["Total Rod Volume (m^3)"] = total_vol
            metrics["Num Rods"] = len(rods)

            magnet = sat.get("magnet", {})
            metrics["Magnet Remanence (T)"] = magnet.get("remanence", np.nan)

    except Exception as e:
        print(f"Warning: Could not parse {toml_path}: {e}")

    return metrics, df

def analyze_directory(directory, threshold):
    csv_files = glob.glob(os.path.join(directory, "*.csv"))
    if not csv_files:
        print(f"Error: No CSV files found in {directory}")
        sys.exit(1)

    print(f"Found {len(csv_files)} simulation runs. Processing...")

    results =[]
    time_series_data = {}

    for csv_file in sorted(csv_files):
        toml_file = csv_file.replace(".csv", ".toml")
        if not os.path.exists(toml_file):
            print(f"Skipping {csv_file}: Missing corresponding .toml file.")
            continue

        metrics, df = process_run(csv_file, toml_file, threshold)
        if metrics is not None:
            results.append(metrics)
            if "time" in df.columns and calculate_magnitude(df, "w") is not None:
                step = max(1, len(df) // 1000)
                time_series_data[metrics["Run"]] = {
                    "time_days": df["time"].iloc[::step] / 86400.0,
                    "w_mag": calculate_magnitude(df, "w").iloc[::step]
                }

    df_results = pd.DataFrame(results)
    return df_results, time_series_data

def safe_corr(df, col1, col2):
    if df[col1].std() == 0 or df[col2].std() == 0:
        return np.nan
    return df[col1].corr(df[col2])

def print_summary(df_results, threshold):
    print("\n" + "="*70)
    print(f" MONTE CARLO SIMULATION SUMMARY (Threshold: {threshold} rad/s)")
    print("="*70)

    total_runs = len(df_results)
    deorbited = df_results["Deorbited"].sum()
    stabilized = df_results["Settling Time (Days)"].notna().sum()

    print(f"Total Successful Integrations : {total_runs}")
    print(f"Spacecraft Stabilized         : {stabilized} / {total_runs} ({stabilized/total_runs*100:.1f}%)")
    print(f"Spacecraft Deorbited Early    : {deorbited} / {total_runs} ({deorbited/total_runs*100:.1f}%)")

    print("\n--- DETUMBLING PERFORMANCE ---")
    if stabilized > 0:
        settle_times = df_results["Settling Time (Days)"].dropna()
        best_run = df_results.loc[settle_times.idxmin()]
        worst_run = df_results.loc[settle_times.idxmax()]

        print(f"Average Settling Time : {settle_times.mean():.2f} Days")
        print(f"Median Settling Time  : {settle_times.median():.2f} Days")
        print(f"Fastest Detumble      : {best_run['Settling Time (Days)']:.2f} Days (Run {best_run['Run']})")
        print(f"Slowest Detumble      : {worst_run['Settling Time (Days)']:.2f} Days (Run {worst_run['Run']})")
    else:
        print("No simulations successfully stabilized within the given threshold.")

    print("\n--- PARAMETER CORRELATIONS (w.r.t Settling Time) ---")
    if stabilized > 2:
        df_stab = df_results.dropna(subset=["Settling Time (Days)"])
        corr_mass = safe_corr(df_stab, "Mass (kg)", "Settling Time (Days)")
        corr_vol = safe_corr(df_stab, "Total Rod Volume (m^3)", "Settling Time (Days)")
        corr_mag = safe_corr(df_stab, "Magnet Remanence (T)", "Settling Time (Days)")

        print(f"Mass vs Settling Time        : {corr_mass:+.3f}" if not np.isnan(corr_mass) else "Mass vs Settling Time        : N/A (No Variance)")
        print(f"Hyst. Volume vs Settling Time: {corr_vol:+.3f}" if not np.isnan(corr_vol) else "Hyst. Volume vs Settling Time: N/A (No Variance)")
        print(f"Magnet Remanence vs Settling : {corr_mag:+.3f}" if not np.isnan(corr_mag) else "Magnet Remanence vs Settling : N/A (No Variance)")

    print("="*70)

def plot_results(df_results, time_series_data, threshold, output_dir):
    tz = pytz.timezone('Europe/Bratislava')
    current_time = datetime.now(tz).strftime("%Y-%m-%d %H:%M %Z")

    # --- 1. Enhanced Static Envelope Plot ---
    # Make the plot wider to accommodate the legend on the right side
    fig, ax = plt.subplots(figsize=(14, 7))

    runs = sorted(time_series_data.keys())
    cmap = plt.get_cmap('tab20') # 20 distinct colors

    # Identify best and worst runs for highlighting
    best_run, worst_run = None, None
    if not df_results.empty:
        settle_times = df_results["Settling Time (Days)"].dropna()
        if not settle_times.empty:
            best_run = df_results.loc[settle_times.idxmin(), "Run"]
            worst_run = df_results.loc[settle_times.idxmax(), "Run"]

    # Plot all runs
    for i, run_name in enumerate(runs):
        data = time_series_data[run_name]

        # Highlight specific runs
        if run_name == best_run:
            ax.plot(data["time_days"], data["w_mag"], color='green', linewidth=2.5, zorder=10, label=f"Run {run_name} (Fastest)")
        elif run_name == worst_run:
            ax.plot(data["time_days"], data["w_mag"], color='red', linewidth=2.5, zorder=10, label=f"Run {run_name} (Slowest)")
        else:
            ax.plot(data["time_days"], data["w_mag"], color=cmap(i % 20), alpha=0.6, linewidth=1.0, label=f"Run {run_name}")

    ax.axhline(threshold, color='black', linestyle='--', linewidth=2, label=f'Threshold ({threshold} rad/s)')
    ax.set_title(f"Angular Velocity Decay (All Runs)\nGenerated: {current_time}", fontsize=12)
    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Angular Velocity Magnitude (rad/s)")
    ax.set_yscale("log")
    ax.grid(True, which="both", ls="--", alpha=0.5)

    # Place legend completely outside the plot area. Dynamically adjust columns based on # of runs.
    ncol = max(1, len(runs) // 20 + 1)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', ncol=ncol, fontsize='x-small')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mc_tumble_envelope.png"), dpi=150)
    plt.close()

    # --- 1b. Interactive Plotly Export (Optional) ---
    try:
        import plotly.graph_objects as go
        fig_html = go.Figure()

        for run_name in runs:
            data = time_series_data[run_name]

            # Highlight best and worst in HTML too
            line_opts = {'width': 1, 'color': None}
            name_label = f"Run {run_name}"
            opacity = 0.5

            if run_name == best_run:
                line_opts = {'width': 3, 'color': 'green'}
                name_label += " (Fastest)"
                opacity = 1.0
            elif run_name == worst_run:
                line_opts = {'width': 3, 'color': 'red'}
                name_label += " (Slowest)"
                opacity = 1.0

            fig_html.add_trace(go.Scatter(
                x=data["time_days"], y=data["w_mag"],
                mode='lines', name=name_label,
                line=line_opts, opacity=opacity
            ))

        fig_html.add_hline(y=threshold, line_dash="dash", line_color="black", annotation_text="Threshold")
        fig_html.update_layout(
            title="Interactive Angular Velocity Decay (Hover over lines)",
            xaxis_title="Time (Days)",
            yaxis_title="Angular Velocity Magnitude (rad/s)",
            yaxis_type="log",
            hovermode="closest",
            template="plotly_white"
        )
        html_path = os.path.join(output_dir, "mc_tumble_interactive.html")
        fig_html.write_html(html_path)
        print(f"Generated interactive plot: {html_path}")
    except ImportError:
        pass # Plotly is not installed, fail silently.

    # --- 2. Scatter Plot: Settling Time vs Total Hysteresis Volume ---
    plt.figure(figsize=(8, 5))
    scatter = plt.scatter(
        df_results["Total Rod Volume (m^3)"] * 1e7,
        df_results["Settling Time (Days)"],
        c=df_results["Mass (kg)"],
        cmap='viridis',
        s=100, edgecolor='k', alpha=0.8
    )
    plt.colorbar(scatter, label="Spacecraft Mass (kg)")
    plt.title("Detumbling Performance vs. Hysteresis Material Volume")
    plt.xlabel(r"Total Rod Volume ($\times 10^{-7} m^3$)")
    plt.ylabel("Settling Time (Days)")
    plt.grid(True, ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mc_volume_vs_settling.png"), dpi=150)
    plt.close()

    # --- 3. Histogram of Settling Times ---
    plt.figure(figsize=(8, 5))
    settle_times = df_results["Settling Time (Days)"].dropna()
    if not settle_times.empty:
        plt.hist(settle_times, bins=10, color='skyblue', edgecolor='black')
        plt.axvline(settle_times.mean(), color='red', linestyle='dashed', linewidth=2, label=f"Mean: {settle_times.mean():.1f} d")
        plt.title("Distribution of Settling Times")
        plt.xlabel("Settling Time (Days)")
        plt.ylabel("Frequency")
        plt.legend()
        plt.grid(axis='y', alpha=0.75)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "mc_settling_histogram.png"), dpi=150)
    plt.close()

    print(f"Static plots successfully generated and saved to '{output_dir}/'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze AOS Monte Carlo Results")
    parser.add_argument("-d", "--dir", default="analysis", help="Directory containing Monte Carlo CSV and TOML files")
    parser.add_argument("-t", "--threshold", type=float, default=0.02, help="Angular velocity threshold for stability [rad/s] (default: 0.02)")

    args = parser.parse_args()

    if not os.path.exists(args.dir):
        print(f"Directory {args.dir} does not exist.")
        sys.exit(1)

    df_results, time_series_data = analyze_directory(args.dir, args.threshold)

    if df_results.empty:
        print("No valid data could be extracted.")
        sys.exit(1)

    print_summary(df_results, args.threshold)
    plot_results(df_results, time_series_data, args.threshold, args.dir)
