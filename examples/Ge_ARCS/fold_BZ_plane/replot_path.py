import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

plt.rcParams.update({"font.size": 14})

fig, ax = plt.subplots(figsize=(7.5, 3.5))


data = np.load("folding_progess/folded_path_plot_20.npy")


# Plotting
norm = SymLogNorm(linthresh=7e-4, vmin=7e-4, vmax=5e-2)
im = plt.pcolormesh(data.T, norm=norm)

# im = ax.imshow(np.log10(np.flipud(new_data_reshaped.T)), vmin=vmin, vmax=vmax, aspect='auto')
cbar = plt.colorbar(im, pad=0.02)


# Customize ticks and labels
xpos = [40, 78, 114, 148, 180, 210, 238, 264, 288, 310, 330, 346, 358, 366]
for pos in xpos:
    ax.axvline(x=pos, color="black", linewidth=1)
plt.xlim([0, 370])
plt.xticks("")

# Set labels and title
plt.ylabel("Energy (meV)")
ypos = np.arange(0, 80, 20)
yticks = np.arange(0, 40, 10)
fmt = lambda x: "{:.0f}".format(x)
plt.yticks(ypos, [fmt(i) for i in yticks])
plt.ylim([0, 80])

plt.savefig("replot_folded_path_plot_20.png")
