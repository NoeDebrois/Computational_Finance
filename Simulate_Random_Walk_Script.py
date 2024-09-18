import numpy as np
import matplotlib.pyplot as plt
from collections import deque

random_x_y = np.random.standard_normal(size=(2, 1000))
random = np.apply_along_axis(np.cumsum, 1, random_x_y)

# Create a fixed-length deque to store the data points
data_points = deque(maxlen=random.shape[1]) 

# Create an empty plot
fig, ax = plt.subplots()
line = ax.plot([])

#plt.figure(figsize=(8,8))
#plt.title("Random walk")
for k in range(random.shape[1]):
    x, y = random[0,k], random[1,k]
    data_points.append((x, y))
    # Update the plot with the new data points
    x_values = [x for x, y in data_points]
    y_values = [y for x, y in data_points]
    plt.scatter(x_values, y_values, color='red')
    # pause the plot for 0.01s before next point is shown
    plt.pause(0.01)
plt.show()