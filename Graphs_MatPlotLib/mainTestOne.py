import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV data from the first file
file_path_1 = '/Users/natanelsolomonov/Desktop/slr-simple/model-runs/status-quo/output/df_model_1_scInt_ne0.5.csv'
print(file_path_1)
try:
    data_1 = pd.read_csv(file_path_1)
except FileNotFoundError:
    print(f"File not found: {file_path_1}")
    exit()

# Load the CSV data from the second file
file_path_2 = '/Users/natanelsolomonov/Desktop/slr-simple/model-runs/status-quo-2/output/df_model_1_scInt_ne0.5.csv'
print(file_path_2)
try:
    data_2 = pd.read_csv(file_path_2)
except FileNotFoundError:
    print(f"File not found: {file_path_2}")
    exit()

# Load the CSV data from the third file
file_path_3 = '/Users/natanelsolomonov/Desktop/slr-simple/model-runs/status-quo-3/output/df_model_1_scInt_ne0.5.csv'
print(file_path_3)
try:
    data_3 = pd.read_csv(file_path_3)
except FileNotFoundError:
    print(f"File not found: {file_path_3}")
    exit()

# Load the CSV data from the fourth file
file_path_4 = '/Users/natanelsolomonov/Desktop/slr-simple/model-runs/Previous data/df_model_1_scInt_ne0.5_RLModel.csv'
print(file_path_4)
try:
    data_4 = pd.read_csv(file_path_4)
except FileNotFoundError:
    print(f"File not found: {file_path_4}")
    exit()

# Extract the necessary columns from the first file
time_1 = data_1['time']
n_unoccupied_1 = data_1['n_unoccupied']

# Extract the necessary columns from the second file
time_2 = data_2['time']
n_unoccupied_2 = data_2['n_unoccupied']

# Extract the necessary columns from the third file
time_3 = data_3['time']
n_unoccupied_3 = data_3['n_unoccupied']

# Extract the necessary columns from the fourth file
time_4 = data_4['time']
n_unoccupied_4 = data_4['n_unoccupied']

# Define a function to calculate the moving average
def moving_average(x, w=10):
    return np.convolve(x, np.ones(w), 'valid') / w

# Apply moving average smoothing
time_1_smooth = time_1[len(time_1) - len(moving_average(n_unoccupied_1)):]
n_unoccupied_1_smooth = moving_average(n_unoccupied_1)

time_2_smooth = time_2[len(time_2) - len(moving_average(n_unoccupied_2)):]
n_unoccupied_2_smooth = moving_average(n_unoccupied_2)

time_3_smooth = time_3[len(time_3) - len(moving_average(n_unoccupied_3)):]
n_unoccupied_3_smooth = moving_average(n_unoccupied_3)

time_4_smooth = time_4[len(time_4) - len(moving_average(n_unoccupied_4)):]
n_unoccupied_4_smooth = moving_average(n_unoccupied_4)

# Plot the smooth data from the first file
# plt.plot(time_1_smooth, n_unoccupied_1_smooth, marker='', linestyle='-', color='b', label="Test for 50% Exposure (status-quo)")
plt.plot(time_1, n_unoccupied_1, marker='', linestyle='-', color='b', label="Test for 10% Exposure (status-quo)")

# Plot the smooth data from the second file
# plt.plot(time_2_smooth, n_unoccupied_2_smooth, marker='', linestyle='--', color='r', label="Test for 30% Exposure (status-quo-2)")
plt.plot(time_2, n_unoccupied_2, marker='', linestyle='--', color='r', label="Test for 30% Exposure (status-quo-2)")

# Plot the smooth data from the third file
# plt.plot(time_3_smooth, n_unoccupied_3_smooth, marker='', linestyle='-', color='y', label="Test for 10% Exposure (status-quo-3)")
plt.plot(time_3, n_unoccupied_3, marker='', linestyle='-', color='y', label="Test for 50% Exposure (status-quo-3)")

# Plot the smooth data from the fourth file
# plt.plot(time_4_smooth, n_unoccupied_4_smooth, marker='', linestyle='-', color='purple', label="Test for RL Model Exposure (Previous data)")
plt.plot(time_4, n_unoccupied_4, marker='', linestyle='-', color='purple', label="Test for RL Model Exposure (Previous data)")

# Add titles and labels
plt.title('Number of Unoccupied Buildings Over Time')
plt.xlabel('Time (Years)')
plt.ylabel('Number of Unoccupied Buildings')

# Add legend
plt.legend()

# Show the plot
plt.grid(True)
plt.show()