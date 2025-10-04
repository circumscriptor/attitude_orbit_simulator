import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# Check if CSV file path is provided as command line argument
if len(sys.argv) < 2:
    print("Usage: python plot_hysteresis.py <your_data.csv>")
    sys.exit(1)

# Get CSV file path from command line argument
csv_file = sys.argv[1]

# Load CSV data
data = pd.read_csv(csv_file)

# Extract columns
H_Am = data['H_Am']
B_T = data['B_T']
M_Am = data['M_Am']

# Set subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Subfigure 1: B_T vs H_Am (hysteresis loop)
ax1.plot(H_Am, B_T, color='blue')
ax1.set_xlabel('Magnetizing Field Amplitude $H_{Am}$ (A/m)')
ax1.set_ylabel('Total Magnetic Flux Density $B_T$ (T)')
ax1.set_title('Hysteresis Curve ($B_T$ vs $H_{Am}$)')
ax1.grid(True)

# Subfigure 2: M_Am vs H_Am (magnetization)
ax2.plot(H_Am, M_Am, color='red')
ax2.set_xlabel('Magnetizing Field Amplitude $H_{Am}$ (A/m)')
ax2.set_ylabel('Magnetization Amplitude $M_{Am}$ (A/m)')
ax2.set_title('Magnetization Curve ($M_{Am}$ vs $H_{Am}$)')
ax2.grid(True)

# Show
plt.tight_layout()
# plt.show()

filename_wo_ext = os.path.splitext(csv_file)[0]
save_path = f"{filename_wo_ext}_curve.png"
plt.savefig(save_path, dpi=300)
