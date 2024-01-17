import matplotlib.pyplot as plt

# Assuming 'file7' is already open and read
# Replace 'file7' with the actual file object or filename

# Lists to store the extracted data
weight = []
whist = []

# Read data from the file
with open('../hamsa_maam_outputs/hbweight_histogram_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        weight.append(float(values[0]))
        whist.append(int(values[1]))

# Plotting the histogram for hbweight histogram
plt.bar(weight, whist, width=0.05, label='hb-weight : Normal', alpha=0.7, color = "blue")


weight = []
whist = []

# Read data from the file
with open('../shaun_outputs/kdtree/hbweight_histogram_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        weight.append(float(values[0]))
        whist.append(int(values[1]))

# Plotting the histogram for hbweight histogram
plt.bar(weight, whist, width=0.05, label='hb-weight : kdtree', alpha=0.7, color = "red")

weight = []
whist = []

# Read data from the file
with open('../shaun_outputs/threading/hbweight_histogram_.dat', 'r') as file7:
    heading = True
    for line in file7:
        if heading:
            heading = False
            continue
        values = line.split()
        weight.append(float(values[0]))
        whist.append(int(values[1]))

# Plotting the histogram for hbweight histogram
plt.bar(weight, whist, width=0.05, label='hb-weight : threading', alpha=0.7, color = "green")

# Setting labels and title
plt.xlabel('weight')
plt.ylabel('hb-weight')
plt.title('hb-weight histogram')
plt.legend()

# Show the plot
plt.savefig("Fig2.png")
plt.show()