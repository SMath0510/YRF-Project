import matplotlib.pyplot as plt

# Assuming 'file7' is already open and read
# Replace 'file7' with the actual file object or filename

# Lists to store the extracted data
path_length = []
line_path = []

# Read data from the file
with open('Pathlength_histogram_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        path_length.append(int(values[0]))
        line_path.append(int(values[1]))

# Plotting the histogram for hbweight histogram
plt.plot(path_length, line_path, 'go-', label='line-path', alpha=0.7)
print(path_length, line_path)

# Setting labels and title
plt.xlabel('pathlength')
plt.ylabel('linepath')
plt.title('pathlength-histogram')
plt.legend()

# Show the plot
plt.savefig("Fig3.png")
plt.show()
