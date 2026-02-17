import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from mantid.simpleapi import *

# === USER CONFIGURATION ===

slices_folder = 'slices/fully_folded/'
data_path = 'folding_progess/folded_path_plot_37.npy'
pathSQE_module_path = 'pathSQE_input.py'

# Energy axis config
E_ylim = (0, 20)
E_tick_step = 5

# Color scaling
norm = Normalize(vmin=0, vmax=2e-2) # LogNorm(vmin=2e-4, vmax=2e-2)
colormap = 'plasma'



# === LOAD INPUT FILE ===
from importlib.util import spec_from_file_location, module_from_spec
spec = spec_from_file_location("pathSQE_input", pathSQE_module_path)
pathSQE_input = module_from_spec(spec)
spec.loader.exec_module(pathSQE_input)
params = pathSQE_input.define_pathSQE_params()

# Parse energy binning info
E_min, E_step, E_max = [float(x) for x in params['E bins'].split(',')]
E_bins = np.arange(E_min, E_max + E_step, E_step)

# === LOAD DATA (for Y-dimension only) ===
data = np.load(data_path)
print('Data shape:', data.shape)

# === FIGURE SETUP ===
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(7.5, 3.5))

# === X-AXIS SEGMENT LOADING ===
segment_list = params['user defined Qpoints']['path']

dsl_fold = []
cumulative_x = 0
for start, end in segment_list:
    segment_name = f"{start}2{end}"
    filename = os.path.join(slices_folder, f"{segment_name}_folded.nxs")

    if not os.path.isfile(filename):
        print(f"WARNING: Missing file {filename}")
        continue

    ws_name = f"tmp_{segment_name}"
    LoadMD(Filename=filename, OutputWorkspace=ws_name)
    signal_array = np.squeeze(mtd[ws_name].getSignalArray())
    x_len = signal_array.shape[0]

    dsl_fold.append({
        'seg_start_name': start,
        'seg_end_name': end,
        'x_start': cumulative_x,
        'x_end': cumulative_x + x_len
    })

    cumulative_x += x_len

print("Segment info:", dsl_fold)

# === PLOTTING ===
X = np.arange(data.shape[0] + 1)
Y = E_bins[:data.shape[1] + 1]
im = ax.pcolormesh(X, Y, data.T, norm=norm, cmap=colormap)
cbar = plt.colorbar(im, pad=0.02)

# === Y-AXIS CUSTOMIZATION ===
yticks = np.arange(E_ylim[0], E_ylim[1] + E_tick_step, E_tick_step)
ax.set_yticks(yticks)
ax.set_yticklabels([f"{y:.0f}" for y in yticks])
ax.set_ylabel("Energy (meV)")
ax.set_ylim(E_ylim)

# === X-TICK LABELING ===
known_mathtext_symbols = {
    'Gamma', 'Delta', 'Sigma', 'Pi', 'Lambda', 'Omega', 'Theta',
    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta',
    'theta', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'pi', 'rho',
    'sigma', 'tau', 'phi', 'chi', 'psi', 'omega'
}

def safe_label(pt):
    if pt in known_mathtext_symbols:
        return rf'$\{pt}$'
    elif len(pt) == 1:
        return rf'${pt}$'
    else:
        return pt

def set_seg_xlabels(i, dsl_fold):
    ticks = []
    labels = []

    if i == len(dsl_fold) - 1:
        prev_end = dsl_fold[i - 1]['seg_end_name']
        curr_start = dsl_fold[i]['seg_start_name']
        curr_end = dsl_fold[i]['seg_end_name']
        pt_labels = [safe_label(pt) for pt in [prev_end, curr_start, curr_end]]
        if curr_start == prev_end:
            ticks = [dsl_fold[i]['x_start'], dsl_fold[i]['x_end']]
            labels = [pt_labels[1], pt_labels[2]]
        else:
            ticks = [dsl_fold[i]['x_start'], dsl_fold[i]['x_end']]
            labels = [f'{pt_labels[0]},{pt_labels[1]}', pt_labels[2]]

    elif i == 0:
        curr_start = dsl_fold[i]['seg_start_name']
        pt_label = safe_label(curr_start)
        ticks = [dsl_fold[i]['x_start']]
        labels = [pt_label]

    else:
        prev_end = dsl_fold[i - 1]['seg_end_name']
        curr_start = dsl_fold[i]['seg_start_name']
        pt_labels = [safe_label(pt) for pt in [prev_end, curr_start]]
        if curr_start == prev_end:
            ticks = [dsl_fold[i]['x_start']]
            labels = [pt_labels[1]]
        else:
            ticks = [dsl_fold[i]['x_start']]
            labels = [f'{pt_labels[0]},{pt_labels[1]}']

    return ticks, labels

# Draw vertical lines
for seg in dsl_fold:
    ax.axvline(x=seg['x_start'], color='black', linewidth=1)
ax.axvline(x=dsl_fold[-1]['x_end'], color='black', linewidth=1)

# Collect ticks + labels
xtick_positions = []
xtick_labels = []

for i in range(len(dsl_fold)):
    ticks, labels = set_seg_xlabels(i, dsl_fold)
    xtick_positions.extend(ticks)
    xtick_labels.extend(labels)

ax.set_xticks(xtick_positions)
ax.set_xticklabels(xtick_labels)
ax.set_xlim([0, dsl_fold[-1]['x_end']])
ax.tick_params(axis='x', labelrotation=0)

# === FINALIZE ===
plt.tight_layout()
plt.savefig('replot_folded_path.png', dpi=300)
print("Saved: replot_folded_path.png")
