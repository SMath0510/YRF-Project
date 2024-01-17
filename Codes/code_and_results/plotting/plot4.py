import matplotlib.pyplot as plt
import numpy as np
# Assuming 'file7' is already open and read
# Replace 'file7' with the actual file object or filename

# Lists to store the extracted data
frame_no = []
hbonds = []

# Read data from the file
with open('../hamsa_maam_outputs/Network_Stats_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        frame_no.append(int(values[0]))
        hbonds.append(int(values[1]))

# Calculate the number of bins based on the bin width
bin_width = 4
num_bins = int((max(hbonds) - min(hbonds)) / bin_width)
plt.hist(hbonds, bins=num_bins, range=(min(hbonds), max(hbonds)), color='blue', edgecolor='black', alpha=0.7, label="Normal")

# Lists to store the extracted data
frame_no = []
hbonds = []

# Read data from the file
with open('../shaun_outputs/kdtree/Network_Stats_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        frame_no.append(int(values[0]))
        hbonds.append(int(values[1]))

# Calculate the number of bins based on the bin width
bin_width = 4
num_bins = int((max(hbonds) - min(hbonds)) / bin_width)
plt.hist(hbonds, bins=num_bins, range=(min(hbonds), max(hbonds)), color='green', edgecolor='black', alpha=0.7, label="KDTree")

# Lists to store the extracted data
frame_no = []
hbonds = []

# Read data from the file
with open('../shaun_outputs/threading/Network_Stats_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        frame_no.append(int(values[0]))
        hbonds.append(int(values[1]))

# Calculate the number of bins based on the bin width
bin_width = 4
num_bins = int((max(hbonds) - min(hbonds)) / bin_width)
plt.hist(hbonds, bins=num_bins, range=(min(hbonds), max(hbonds)), color='red', edgecolor='black', alpha=0.7, label="Threading")

# Set labels and title
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title(f'Number of HBonds')
plt.legend()

# Show the plot
plt.savefig("Fig4.png")
plt.show()
