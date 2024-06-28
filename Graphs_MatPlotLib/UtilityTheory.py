import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV data from the "status-quo" folder
file_path_1 = '/Users/natanelsolomonov/Desktop/slr-simple/model-runs/Status-quo-2/output/df_model_1_scInt_ne0.5.csv'
print(file_path_1)
try:
    data_1 = pd.read_csv(file_path_1)
except FileNotFoundError:
    print(f"File not found: {file_path_1}")
    exit()

# Load the CSV data from the "Previous data" folder
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

# Extract the necessary columns from the fourth file
time_4 = data_4['time']
n_unoccupied_4 = data_4['n_unoccupied']

# Plot the data from the first file
plt.plot(time_1, n_unoccupied_1, marker='', linestyle='-', color='b', label="Utility Theory Test (status-quo)")

# Plot the data from the fourth file
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