import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


A = np.loadtxt("FluctFitness.dat")
IntMap = np.loadtxt("IntMap.dat")


fig, axs = plt.subplots(1, 2, figsize=(12, 6))


mask_A = (A <0)
axs[0].hist(A[mask_A], bins=30, alpha=0.7, color="k")
axs[0].set_title(f'Mean = {A.mean():.2f}, Stdev = {A.std():.2f}')
axs[0].set_xlabel("Fitness")


B = np.zeros((6,6))
k = 0
for j in range(0,6):
    for i in range(j,6):
        B[i,j] = IntMap[k]
        k += 1
mask = np.triu(np.ones_like(B), k=1)

ax = sns.heatmap(B, mask=mask, vmin=-10, vmax=10, cbar=True, square=True, annot=True, linewidths=0.5, xticklabels=[], yticklabels=[], ax=axs[1])
ax.set_title('Interaction Map')

cbar = ax.collections[0].colorbar
cbar.set_label('Energy [kT]')


plt.tight_layout()
plt.show()
