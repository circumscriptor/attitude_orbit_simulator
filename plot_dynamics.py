import sys
import pandas as pd
import matplotlib.pyplot as plt

# Check if CSV file path is provided as command line argument
if len(sys.argv) < 3:
    print("Usage: python plot_dynamics.py <your_data.csv> <number_of_rods>")
    sys.exit(1)

# Get CSV file path from command line argument
csv_file = sys.argv[1]

try:
    number_of_rods = int(sys.argv[2])
except ValueError:
    print("Error: number_of_rods must be an integer.")
    sys.exit(1)

# Load CSV data
data = pd.read_csv(csv_file)

# Plot settings
# plt.style.use('seaborn-darkgrid')
fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

# 1. Plot Quaternion components vs time
axs[0].plot(data['time'], data['q_w'], label='q_w')
axs[0].plot(data['time'], data['q_x'], label='q_x')
axs[0].plot(data['time'], data['q_y'], label='q_y')
axs[0].plot(data['time'], data['q_z'], label='q_z')
axs[0].set_ylabel('Quaternion')
axs[0].set_title('Attitude Quaternion Components')
axs[0].legend()

# 2. Plot angular velocity components vs time
# axs[1].plot(data['time'], data['w'], label='w')
axs[1].plot(data['time'], data['w_x'], label='w_x')
axs[1].plot(data['time'], data['w_y'], label='w_y')
axs[1].plot(data['time'], data['w_z'], label='w_z')
axs[1].set_ylabel('Angular Velocity (rad/s)')
axs[1].set_title('Angular Velocity Components')
axs[1].legend()

# 3. Plot magnetization components vs time
for i in range(1, number_of_rods + 1):
    column_name = f'M_{i}'
    if column_name in data.columns:
        axs[2].plot(data['time'], data[column_name], label=column_name)
    else:
        print(f"Warning: Column '{column_name}' not found in the CSV file.")

axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Magnetization')
axs[2].set_title('Hysteresis Rod Magnetizations')
axs[2].legend()

plt.tight_layout()
plt.show()
