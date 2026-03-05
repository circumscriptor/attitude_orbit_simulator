import pandas as pd
import numpy as np
import glob
from scipy.interpolate import make_interp_spline

# Force matplotlib to not use any Xwindows/Qt backend.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator

# Configuration
materials = {
    "vacodur50": "Vacodur 50",
    "permalloy": "Permalloy",
    "fe50ni50": "Fe50-Ni50",
    "amorphous_fe": "Amorphous (Fe)",
    "amorphous_co": "Amorphous (Co)",
    "finement": "FINEMENT"
}

THRESHOLD = 0.02
rod_volumes = np.linspace(5e-8, 5e-7, 20)
RESULTS = []

# Create a fixed color mapping for each material
prop_cycle = plt.rcParams['axes.prop_cycle']
colors_list = prop_cycle.by_key()['color']
color_map = {label: colors_list[i % len(colors_list)] for i, label in enumerate(materials.values())}

# Data Processing
for folder, label in materials.items():
    files = sorted(glob.glob(f"analysis/{folder}/sweep_rod_volume_*.csv"))

    for i, file_path in enumerate(files):
        try:
            df = pd.read_csv(file_path)
            df['w_mag'] = np.sqrt(df['w_x']**2 + df['w_y']**2 + df['w_z']**2)

            final_w = df['w_mag'].iloc[-1]

            mask_above = df['w_mag'] > THRESHOLD
            if not mask_above.any():
                stab_time = df['time'].iloc[0]
            else:
                last_idx_above = df.index[mask_above][-1]
                stab_time = df['time'].iloc[last_idx_above + 1] if last_idx_above < len(df)-1 else np.nan

            RESULTS.append({
                "Material": label,
                "Volume": rod_volumes[i],
                "Final_W": final_w,
                "Stab_Days": stab_time / 86400
            })
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

res_df = pd.DataFrame(RESULTS)

# --- FIGURE 1: Final Angular Velocity ---
plt.figure(figsize=(10, 6))
ax1 = plt.gca()

for label in materials.values():
    m_data = res_df[res_df['Material'] == label].dropna(subset=['Final_W'])
    x, y = m_data['Volume'].values, m_data['Final_W'].values
    m_color = color_map[label]

    if len(x) > 3:
        x_smooth = np.linspace(x.min(), x.max(), 300)
        spl = make_interp_spline(x, y, k=3)
        plt.plot(x_smooth, spl(x_smooth), label=label, linewidth=1, color=m_color)
    else:
        plt.plot(x, y, label=label, linewidth=1, color=m_color)

plt.axhline(y=THRESHOLD, color='black', linestyle='--', alpha=0.6, label=f'Threshold ({THRESHOLD})')

plt.yscale('log') # Applied as requested
plt.title("Final Angular Velocity Comparison")
plt.ylabel(r"$|\omega_{final}|$ (rad/s)")
plt.xlabel(r"Hysteresis Rod Volume ($m^3$)")
plt.grid(True, which="both", alpha=0.2)
plt.legend(loc='upper right', fontsize='small', ncol=2)

plt.tight_layout()
plt.savefig("hysteresis_materials_final_angular_velocity_comparison.png", dpi=300)

# --- FIGURE 2: Success Histogram ---
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
plt.ylim(0, 22) # Leave space for numbers

# Formatting Y-axis
from matplotlib.ticker import MultipleLocator
ax2.yaxis.set_major_locator(MultipleLocator(4))

# Add numeric labels
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
             f'{int(height)}', ha='center', va='bottom')

plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.setp(ax2.get_xticklabels(), rotation=15, ha='right')

plt.tight_layout()
plt.savefig("hysteresis_materials_stabilization_success_comparison.png", dpi=300)
