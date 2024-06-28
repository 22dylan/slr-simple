import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Base path and file name pattern for the output files
base_path = '/Users/natanelsolomonov/Desktop/slr-simple/model-runs/status-quo/output/'
file_pattern = 'df_model_{}_scInt_ne0.5.csv'

# Initialize lists to store data for averaging
all_times = []
all_unoccupied = []

# Load the CSV data from the "status-quo" folder
file_path_1 = base_path + file_pattern.format(1)
print(file_path_1)
try:
    data_1 = pd.read_csv(file_path_1)
except FileNotFoundError:
    print(f"File not found: {file_path_1}")
    exit()

# Extract the necessary columns from the first file
time_1 = data_1['time']
n_unoccupied_1 = data_1['n_unoccupied']

# Plot the data from the first file
plt.plot(time_1, n_unoccupied_1, marker='', linestyle='-', color='b', linewidth=0.01, label="Utility Theory Test (status-quo)")

# Add data to lists for averaging
all_times.append(time_1)
all_unoccupied.append(n_unoccupied_1)

# Loop through the 100 output files
for i in range(1, 101):
    file_path_4 = base_path + file_pattern.format(i)
    try:
        data_4 = pd.read_csv(file_path_4)
    except FileNotFoundError:
        print(f"File not found: {file_path_4}")
        continue

    # Extract the necessary columns from each file
    time_4 = data_4['time']
    n_unoccupied_4 = data_4['n_unoccupied']

    # Plot the data from each file with a very thin line
    if i == 1:
        plt.plot(time_4, n_unoccupied_4, marker='', linestyle='-', color='g', linewidth=0.5, label="Data from all 100 runs")
    else:
        plt.plot(time_4, n_unoccupied_4, marker='', linestyle='-', color='g', linewidth=0.2)

    # Add data to lists for averaging
    all_times.append(time_4)
    all_unoccupied.append(n_unoccupied_4)

# Calculate the average number of unoccupied buildings at each time point
all_times = np.array(all_times)
all_unoccupied = np.array(all_unoccupied)

# Ensure all_times have the same length by selecting the unique time points from the first file
unique_times = all_times[0]
average_unoccupied = np.mean(all_unoccupied, axis=0)

# Plot the average line
plt.plot(unique_times, average_unoccupied, marker='', linestyle='-', color='blue', linewidth=2, label="Average (status-quo)")

# Load the CSV data from the fourth file
file_path_4 = '/Users/natanelsolomonov/Desktop/slr-simple/model-runs/Previous data/df_model_1_scInt_ne0.5_RLModel.csv'
print(file_path_4)
try:
    data_4 = pd.read_csv(file_path_4)
except FileNotFoundError:
    print(f"File not found: {file_path_4}")
    exit()
# Extract the necessary columns from the fourth file
time_4 = data_4['time']
n_unoccupied_4 = data_4['n_unoccupied']

#Plot Previous Data
plt.plot(time_4, n_unoccupied_4, marker='', linestyle='-', color='y', linewidth=1, label="Previous Data")

# Add titles and labels
plt.title('Number of Unoccupied Buildings Over Time')
plt.xlabel('Time (Years)')
plt.ylabel('Number of Unoccupied Buildings')

# Set limits to zoom into the graph
plt.xlim(30, 70)  # Adjust the range as needed
plt.ylim(700, 2500)  # Adjust the range as needed

# Add legend
plt.legend()

# Show the plot
plt.grid(True)
plt.show()