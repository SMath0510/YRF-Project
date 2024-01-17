import matplotlib.pyplot as plt
import numpy as np
# Assuming 'file7' is already open and read
# Replace 'file7' with the actual file object or filename

# Lists to store the extracted data
frame_no = []
num_paths = []

# Read data from the file
with open('PerFrame_HistCount_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        frame_no.append(int(values[0]))
        num_paths.append(int(values[4]))

# Calculate the number of bins based on the bin width
bin_width = 100
num_bins = int((max(num_paths) - min(num_paths)) / bin_width)

# Create a histogram with specified bin width
mean_value = np.mean(num_paths)
variance_value = np.var(num_paths)

# Print or use the results

plt.hist(num_paths, bins=num_bins, range=(min(num_paths), max(num_paths)), color='blue', edgecolor='black', alpha=0.7)

# Set labels and title
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title(f'Number of Paths\nMean: {mean_value}, Variance: {variance_value}')

# Show the plot
plt.savefig("Fig5.png")
plt.show()
