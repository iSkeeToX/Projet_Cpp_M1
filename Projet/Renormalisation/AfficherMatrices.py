import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize

# Load matrices from text files
matrix_c = np.loadtxt("Matrice_c.txt")
matrix_c = (matrix_c + matrix_c.T) / 2

matrix_d = np.loadtxt("Matrice_d.txt")
matrix_d = (matrix_d + matrix_d.T) / 2

# Find the common color scale range
min_value = min(np.min(matrix_c), np.min(matrix_d))
max_value = max(np.max(matrix_c), np.max(matrix_d))

# Set Seaborn style
sns.set(style="whitegrid")

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
cmap = sns.diverging_palette(20, 230, as_cmap=True)
cmap = 'viridis'
# Plot matrix d with Seaborn color palette
sns.heatmap(matrix_d, cmap=cmap, vmin=min_value, vmax=max_value, ax=axs[0], cbar=False, square=True, linewidths=0.5, cbar_kws={"shrink": .5})
axs[0].set_title("Densité second voisins $< d >_J$")

# Plot matrix c with Seaborn color palette
sns.heatmap(matrix_c, cmap=cmap, vmin=min_value, vmax=max_value, ax=axs[1], cbar=False, square=True, linewidths=0.5, cbar_kws={"shrink": .5}, yticklabels=[])
axs[1].set_title("Densité premier voisins $< c >_{J'}$")

# Add colorbar manually
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = fig.colorbar(axs[1].collections[0], cax=cbar_ax)
cbar.set_label('Common Color Scale')

plt.savefig("Densité")
plt.subplots_adjust(wspace=0.07)
plt.show()

# Load matrices from text files
InteractionMap = np.loadtxt("Ancienne_Interaction_Map.txt")
InteractionMap = (InteractionMap + InteractionMap.T) / 2

NewInteractionMap = np.loadtxt("Nouvelle_Interaction_Map.txt")
NewInteractionMap = (NewInteractionMap + NewInteractionMap.T) / 2
NewInteractionMap = NewInteractionMap/100

min_value = min(np.min(InteractionMap), np.min(NewInteractionMap))
max_value = max(np.max(InteractionMap), np.max(NewInteractionMap))

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
cmap = sns.diverging_palette(230, 20, as_cmap=True)
# Plot matrix d with Seaborn color palette
sns.heatmap(InteractionMap, cmap=cmap, vmin=min_value, vmax=max_value, ax=axs[0], cbar=False, square=True, center=0, linewidths=0.5, cbar_kws={"shrink": .5})
axs[0].set_title("Carte d'interaction J")

# Plot matrix c with Seaborn color palette
sns.heatmap(NewInteractionMap, cmap=cmap, vmin=min_value, vmax=max_value, ax=axs[1], cbar=False, square=True, center=0, linewidths=0.5, cbar_kws={"shrink": .5}, yticklabels=[])
axs[1].set_title("Carte d'interaction J'/100")

# Add colorbar manually
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = fig.colorbar(axs[1].collections[0], cax=cbar_ax)
cbar.set_label('$J_{ab} [kT]')

plt.savefig("InteractionMap")
plt.subplots_adjust(wspace=0.07)
plt.show()