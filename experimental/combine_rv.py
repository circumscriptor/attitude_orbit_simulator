import os
import glob
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.interpolate import CubicSpline

# Force matplotlib to not use any Xwindows/Qt backend.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import hsv_to_rgb

# Configuration
materials = {
    "vacodur50_rv8": "Vacodur 50 [2.158e-7 m3]",
    "permalloy_rv1": "Permalloy [5.000e-8 m3]",
    "permalloy_rv3": "Permalloy [9.737e-8 m3]",
    "permalloy_rv9": "Permalloy [2.395e-7 m3]",
    "fe50ni50_rv1": "Fe50-Ni50 [5.000e-8 m3]",
    "fe50ni50_rv2": "Fe50-Ni50 [7.368e-8 m3]",
    "amorphous_fe_rv1": "Amorphous (Fe) [5.000e-8]",
    "amorphous_fe_rv5": "Amorphous (Fe) [1.447-7]",
    "amorphous_co_rv2": "Amorphous (Co) [7.368e-8]",
    "amorphous_co_rv4": "Amorphous (Co) [1.211e-7]",
    "finement_rv11": "FINEMENT [2.868e-7]",
    "finement_rv16": "FINEMENT [4.053e-7]"
}

THRESHOLD = 0.02
magnet_remanences = np.linspace(0.7, 1.5, 20)
RESULTS = []

plt.rcParams['lines.markersize'] = 4

# Setup Colors
num_materials = len(materials)
hues = np.linspace(0, 1, num_materials, endpoint=False)
colors_list = [hsv_to_rgb((h, 0.9, 0.9)) for h in hues]
color_map = {label: colors_list[i] for i, label in enumerate(materials.values())}

# Convert RGB tuples (0-1) to Plotly-friendly CSS strings 'rgb(R,G,B)'
plotly_color_map = {
    label: f"rgb({int(c[0]*255)}, {int(c[1]*255)}, {int(c[2]*255)})"
    for label, c in color_map.items()
}

# Dictionary to hold time-series data for Plotly: plotly_data[remanence_index][material_label] = df
plotly_data = {i: {} for i in range(len(magnet_remanences))}

# Data Processing
for folder, label in materials.items():
    files = sorted(glob.glob(f"analysis/{folder}/sweep_magnet_remanence_*.csv"))

    for i, file_path in enumerate(files):
        try:
            df = pd.read_csv(file_path)
            df['w_mag'] = np.sqrt(df['w_x']**2 + df['w_y']**2 + df['w_z']**2)
            df['time_days'] = df['time'] / 86400 # Convert to days for plotting

            # Save full time-series data for Plotly later
            plotly_data[i][label] = df[['time_days', 'w_mag']]

            final_w = df['w_mag'].iloc[-1]
            mask_above = df['w_mag'] > THRESHOLD

            if not mask_above.any():
                stab_time = df['time'].iloc[0]
            else:
                last_idx_above = df.index[mask_above][-1]
                stab_time = df['time'].iloc[last_idx_above + 1] if last_idx_above < len(df)-1 else np.nan

            RESULTS.append({
                "Material": label,
                "Remanence": magnet_remanences[i],
                "Final_W": final_w,
                "Stab_Days": stab_time / 86400
            })
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

res_df = pd.DataFrame(RESULTS)

# =========================================================
# --- FIGURE 1: Final Angular Velocity (Matplotlib) ---
# =========================================================
plt.figure(figsize=(10, 6))
ax1 = plt.gca()

for label in materials.values():
    m_data = res_df[res_df['Material'] == label].dropna(subset=['Final_W'])
    if not m_data.empty:
        x, y = m_data['Remanence'].values, m_data['Final_W'].values
        m_color = color_map[label]

        plt.plot(x, y, zorder=1, linestyle=':', linewidth=1, color=m_color)
        plt.scatter(x, y, zorder=2, color=m_color, label=label)

plt.axhline(y=THRESHOLD, zorder=0, color='black', linestyle='--', alpha=0.6, label=f'Threshold ({THRESHOLD})')

plt.yscale('log')
plt.title("Final Angular Velocity Comparison")
plt.ylabel(r"$|\omega_{final}|$ (rad/s)")
plt.xlabel(r"Magnet Remanence ($m^3$)")
plt.grid(True, which="both", alpha=0.2)
plt.legend(loc='lower right', fontsize='small', ncol=2)

plt.tight_layout()
plt.savefig("magnet_remanence_final_angular_velocity_comparison.png", dpi=300)
print("Saved: magnet_remanence_final_angular_velocity_comparison.png")


# =========================================================
# --- FIGURE 2: Success Histogram (Matplotlib) ---
# =========================================================
plt.figure(figsize=(10, 6))
ax2 = plt.gca()

STAB_LIMIT = 14
success_counts = res_df[res_df['Stab_Days'] <= STAB_LIMIT].groupby('Material').size()
counts = [success_counts.get(label, 0) for label in materials.values()]

material_labels = list(materials.values())
bar_colors = [color_map[label] for label in material_labels]
bars = plt.bar(material_labels, counts, color=bar_colors, alpha=0.8, edgecolor='black')

plt.title(f"Successful Stabilizations (Within {STAB_LIMIT} Days)")
plt.ylabel("Number of Runs (Out of 20)")
plt.ylim(0, 22)

ax2.yaxis.set_major_locator(MultipleLocator(4))

for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
             f'{int(height)}', ha='center', va='bottom')

plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.setp(ax2.get_xticklabels(), rotation=15, ha='right')

plt.tight_layout()
plt.savefig("magnet_remanence_stabilization_success_comparison.png", dpi=300)
print("Saved: magnet_remanence_stabilization_success_comparison.png")


# =========================================================
# --- FIGURE 3: Plotly Interactive Plots (Grouped by Sweep Target) ---
# =========================================================
# try:
#     import plotly.graph_objects as go

#     # Create directory for interactive plots so it doesn't clutter root
#     interactive_dir = "interactive_plots"
#     os.makedirs(interactive_dir, exist_ok=True)

#     print(f"\nGenerating interactive Plotly files in '{interactive_dir}/'...")

#     # Iterate through each sweep target value
#     for i, rem_val in enumerate(magnet_remanences):
#         # Skip if no materials were successfully processed for this sweep step
#         if not plotly_data[i]:
#             continue

#         fig_html = go.Figure()

#         # Add trace for every material available at this sweep target
#         for label in materials.values():
#             if label in plotly_data[i]:
#                 data = plotly_data[i][label]

#                 fig_html.add_trace(go.Scatter(
#                     x=data["time_days"],
#                     y=data["w_mag"],
#                     mode='lines',
#                     name=label,
#                     line=dict(color=plotly_color_map[label]),
#                     opacity=0.8
#                 ))

#         # Add the threshold horizontal line
#         fig_html.add_hline(
#             y=THRESHOLD,
#             line_dash="dash",
#             line_color="red",
#             annotation_text=f"Threshold ({THRESHOLD} rad/s)"
#         )

#         fig_html.update_layout(
#             title=f"Angular Velocity Decay (Magnet Remanence = {rem_val:.4g})",
#             xaxis_title="Time (Days)",
#             yaxis_title="Angular Velocity Magnitude (rad/s)",
#             yaxis_type="log",
#             hovermode="closest",
#             template="plotly_white",
#             legend_title_text="Material"
#         )

#         fig_html.update_traces(hoverlabel_namelength=-1)

#         # Save specific HTML file for this sweep step
#         html_path = os.path.join(interactive_dir, f"plot_remanence_step_{i+1:02d}_{rem_val:.4f}.html")
#         fig_html.write_html(html_path)

#     print(f"Success! {len(magnet_remanences)} interactive HTML plots saved.")

# except ImportError:
#     print("\nNote: 'plotly' is not installed. Skipping interactive HTML plots.")
#     print("To enable interactive plots, run: pip install plotly")

print("\nAnalysis complete.")
