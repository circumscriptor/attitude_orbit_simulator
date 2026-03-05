import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator
import glob
from scipy.interpolate import make_interp_spline

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

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 11))

for label in materials.values():
    m_data = res_df[res_df['Material'] == label].dropna(subset=['Final_W'])

    x = m_data['Volume'].values
    y = m_data['Final_W'].values
    m_color = color_map[label]

    if len(x) > 3: # Need points for spline
        x_smooth = np.linspace(x.min(), x.max(), 300)
        spl = make_interp_spline(x, y, k=3)
        y_smooth = spl(x_smooth)
        ax1.plot(x_smooth, y_smooth, label=label, linewidth=1, color=m_color)
    else:
        ax1.plot(x, y, label=label, linewidth=1, color=m_color)

    # ax2.scatter(m_data['Volume'], m_data['Stab_Days'], label=label, s=40, alpha=0.7)

# Add the threshold line to the top axis (ax1)
ax1.axhline(y=THRESHOLD, color='black', linestyle='--', alpha=0.6, label=f'Threshold ({THRESHOLD})')

# Formatting
ax1.set_title("Final Angular Velocity")
ax1.set_ylabel(r"$|\omega_{final}|$ (rad/s)")
ax1.grid(True, alpha=0.2)
ax1.legend(loc='upper right', fontsize='small', ncol=2)

# 1. Filter for runs that stabilized within 14 days
STAB_LIMIT = 14
success_counts = res_df[res_df['Stab_Days'] <= STAB_LIMIT].groupby('Material').size()

# Ensure all materials are present in the series (even if count is 0)
counts = [success_counts.get(label, 0) for label in materials.values()]
material_labels = list(materials.values())

# 2. Plotting (Replace your ax2 code with this)
bar_colors = [color_map[label] for label in material_labels]
bars = ax2.bar(material_labels, counts, color=bar_colors, alpha=0.8, edgecolor='black')

# 3. Formatting
ax2.set_title(f"Successful Stabilizations (Within {STAB_LIMIT} Days)")
ax2.set_ylabel("Number of Runs (Max 20)")
ax2.set_ylim(0, 22) # Extra room for labels
ax2.grid(axis='y', linestyle='--', alpha=0.6)

ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
ax2.yaxis.set_major_locator(MultipleLocator(4))

# Add numeric labels on top of bars
for bar in bars:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
             f'{int(height)}', ha='center', va='bottom', fontweight='bold')

# Rotate x-labels if they overlap
plt.setp(ax2.get_xticklabels(), rotation=15, ha='right')

plt.tight_layout()
plt.savefig("hysteresis_materials_comparison.png", dpi=300)
print("Plot saved as 'hysteresis_materials_comparison.png'")
