import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Generate a sample NumPy array
data = np.array([[1, 2, 3, 4],
                 [5, 6, 7, 8],
                 [9, 10, 11, 12]])

# Create a Seaborn clustermap with hexagonal cells
sns.clustermap(data, cmap="YlGnBu", linewidths=.5)

# Show the plot
plt.show()
