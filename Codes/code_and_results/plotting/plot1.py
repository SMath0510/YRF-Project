import matplotlib.pyplot as plt

# Lists to store the extracted data
closeness_frame_numbers = []
closeness_values = []
betweenness_frame_numbers = []
betweenness_values = []

# Read data from the file
with open('../hamsa_maam_outputs/Centrality_histogram_.dat', 'r') as file:
    heading = True
    for line in file:
        if heading:
            heading = False
            continue
        values = line.split()
        closeness_frame_numbers.append(float(values[0]))
        closeness_values.append(int(values[1]))
        betweenness_frame_numbers.append(float(values[2]))
        betweenness_values.append(int(values[3]))

# Plotting the histogram for Closeness Centrality
# plt.bar(closeness_frame_numbers, closeness_values, width=0.01, label='Closeness Centrality', alpha=0.7)
plt.plot(closeness_frame_numbers, closeness_values, 'o-r',  label='Closeness Centrality : Normal', alpha=0.7)

# Plotting the histogram for Betweenness Centrality
# plt.bar(betweenness_frame_numbers, betweenness_values, width=0.005, label='Betweenness Centrality', alpha=0.7)
plt.plot(betweenness_frame_numbers, betweenness_values, 'o-b', label='Betweenness Centrality : Normal', alpha=0.7)


# Lists to store the extracted data
closeness_frame_numbers = []
closeness_values = []
betweenness_frame_numbers = []
betweenness_values = []

# Read data from the file
with open('../shaun_outputs/kdtree/Centrality_histogram_.dat', 'r') as file:
    heading = True
    for line in file:
        if heading:
            heading = False
            continue
        values = line.split()
        closeness_frame_numbers.append(float(values[0]))
        closeness_values.append(int(values[1]))
        betweenness_frame_numbers.append(float(values[2]))
        betweenness_values.append(int(values[3]))

# Plotting the histogram for Closeness Centrality
# plt.bar(closeness_frame_numbers, closeness_values, width=0.01, label='Closeness Centrality', alpha=0.7)
plt.plot(closeness_frame_numbers, closeness_values, 'o-c',  label='Closeness Centrality : kdTree', alpha=0.7)

# Plotting the histogram for Betweenness Centrality
# plt.bar(betweenness_frame_numbers, betweenness_values, width=0.005, label='Betweenness Centrality', alpha=0.7)
plt.plot(betweenness_frame_numbers, betweenness_values, 'o-y', label='Betweenness Centrality : kdTree', alpha=0.7)


# Lists to store the extracted data
closeness_frame_numbers = []
closeness_values = []
betweenness_frame_numbers = []
betweenness_values = []

# Read data from the file
with open('../shaun_outputs/threading/Centrality_histogram_.dat', 'r') as file:
    heading = True
    for line in file:
        if heading:
            heading = False
            continue
        values = line.split()
        closeness_frame_numbers.append(float(values[0]))
        closeness_values.append(int(values[1]))
        betweenness_frame_numbers.append(float(values[2]))
        betweenness_values.append(int(values[3]))

# Plotting the histogram for Closeness Centrality
# plt.bar(closeness_frame_numbers, closeness_values, width=0.01, label='Closeness Centrality', alpha=0.7)
plt.plot(closeness_frame_numbers, closeness_values, 'o-g',  label='Closeness Centrality : Threading', alpha=0.7)

# Plotting the histogram for Betweenness Centrality
# plt.bar(betweenness_frame_numbers, betweenness_values, width=0.005, label='Betweenness Centrality', alpha=0.7)
plt.plot(betweenness_frame_numbers, betweenness_values, 'o-m', label='Betweenness Centrality : Threading', alpha=0.7)

# Setting labels and title
plt.xlabel('Frame Number')
plt.ylabel('Centrality Values')
plt.title('Centrality Values vs Frame Number')
plt.legend()

# Show the plot
plt.savefig("Fig1.png")
plt.show()
