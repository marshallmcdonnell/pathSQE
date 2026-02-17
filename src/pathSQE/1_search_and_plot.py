import os
import re
import numpy as np
from mantid.simpleapi import *
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm



# Set global font sizes for consistency
plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10})



def build_search_regex(search_dict):
    """Build regex to match Mantid-generated filenames using flexible keys."""
    pattern = r"BZ"
    BZ = search_dict.get("BZ")
    pattern += re.escape(str(BZ)) if BZ is not None else r"\[.*?\]"
    pattern += "_"

    seg_name = search_dict.get("seg_name")
    pattern += re.escape(seg_name) if seg_name is not None else r".+?"
    pattern += "_"

    pt1 = search_dict.get("pt1")
    pt2 = search_dict.get("pt2")
    pattern += re.escape(str(pt1)) if pt1 is not None else r"\[.*?\]"
    pattern += "to"
    pattern += re.escape(str(pt2)) if pt2 is not None else r"\[.*?\]"

    pattern += "_ID"
    sliceID = search_dict.get("sliceID")
    pattern += str(sliceID) if sliceID is not None else r"\d+"
    pattern += "_data\.nxs$"

    return pattern

def search_files(folder, search_dict):
    pattern = build_search_regex(search_dict)
    regex = re.compile(pattern)
    return [os.path.join(folder, f) for f in os.listdir(folder) if regex.match(f)]

def load_and_normalize(file_path, ws_prefix="slice"):
    """Loads both data and norm for a given file and returns normalized workspace name."""
    norm_path = file_path.replace("_data.nxs", "_norm.nxs")
    base_name = os.path.splitext(os.path.basename(file_path))[0].replace("_data", "")
    
    ws_data = f"{ws_prefix}_data_{base_name}"
    ws_norm = f"{ws_prefix}_norm_{base_name}"
    ws_plot = f"{ws_prefix}_plot_{base_name}"

    LoadMD(Filename=file_path, OutputWorkspace=ws_data, LoadHistory=False)
    LoadMD(Filename=norm_path, OutputWorkspace=ws_norm, LoadHistory=False)

    # Divide to normalize
    ws_plot_obj = mtd[ws_data] / mtd[ws_norm]
    
    return ws_plot_obj

def batch_search_and_load(folder, search_terms):
    """Batch search and load slices using a list of flexible search dictionaries."""
    loaded_ws = []

    for i, term in enumerate(search_terms):
        results = search_files(folder, term)

        # Build a nice string summary of provided keys
        summary_parts = []
        for key in ["BZ", "seg_name", "pt1", "pt2", "sliceID"]:
            if key in term:
                summary_parts.append(f"{key}={term[key]}")
        summary = ", ".join(summary_parts) if summary_parts else "no keys specified"

        if len(results) == 0:
            print(f"[{i}] No match found for {summary}")
        elif len(results) > 2:
            print(f"[{i}] {len(results)} matches found for {summary}:")
            for r in results:
                print("   ", os.path.basename(r))
        elif len(results) == 2:
            print(f"[{i}] {len(results)} matches found for {summary}:")
            for r in results:
                print("   ", os.path.basename(r))
            plot_ws = load_and_normalize(results[0], ws_prefix=f"seg{i}")
            loaded_ws.append(plot_ws)
        else:
            print(f"[{i}] Found 1 match for {summary}:")
            print("   ", os.path.basename(results[0]))
            plot_ws = load_and_normalize(results[0], ws_prefix=f"seg{i}")
            loaded_ws.append(plot_ws)

    return loaded_ws


def generate_search_terms(user_qpoints, BZ_offset, seg_name=None, sliceID=None):
    """
    Generate search terms based on user-defined Q-path and BZ offset.
    """
    search_terms = []

    for label1, label2 in user_qpoints["path"]:
        pt1 = np.array(user_qpoints["point_coords"][label1])
        pt2 = np.array(user_qpoints["point_coords"][label2])

        entry = {
            "pt1": BZ_offset + pt1,
            "pt2": BZ_offset + pt2
        }
        if seg_name is not None:
            entry["seg_name"] = seg_name
        if sliceID is not None:
            entry["sliceID"] = sliceID

        search_terms.append(entry)

    return search_terms


################################# Searching for specific slice IDs #########################################
'''
# Specific with IDs
search_terms = [
    {"sliceID": 43},
    {"sliceID": 43},
    {"sliceID": 43}
]

folder = "/SNS/ARCS/IPTS-13861/shared/Aiden/comprehensive/pathSQE_testing_Ge/paperData2_5K_40meV/all_slices"
plot_workspaces = batch_search_and_load(folder, search_terms)
'''


################################# Searching for same path in many BZs #########################################
# path search results in given BZ
user_defined_Qpoints = {
    "path": [('Gamma', 'X'), ('X', 'U'), ('K', 'Gamma'), ('Gamma', 'L'), ('L', 'W'), ('W', 'X')],
    "1d_points": ['K', 'L'],
    "point_coords": {
        'Gamma': [0.0, 0.0, 0.0], 'X': [0.0, 1.0, 0.0], 'L': [0.5, 0.5, 0.5],
        'W': [0.5, 1.0, 0.0], 'K': [0.75, 0.75, 0.0], 'U': [0.25, 1.0, 0.25]
    }
}
folder = "/SNS/ARCS/IPTS-13861/shared/Aiden/comprehensive/pathSQE_testing_Ge/paperData2_5K_40meV/all_slices/"


BZ_list = np.array([
    [4., 2., 0.],    [2, 4, 0] ,    [4, 0, 0],    [2, 2, 0],
    [4, -2, 0],    [3, 3, 1],    [3, -1, 1],    [2, 0, 0],
    [3, 3, -1],    [3, 1, -1],    [3, 1, 1],    [5, -1, -1],
    [5, 1, -1],    [1, 5, -1],    [3, -1, -1],    [1, 3, -1],
    [1, 3, 1],    [0, 2, 0],    [4, 4, 0],    [5, -1, 1],
    [1, 5, 1]
])

for BZ_offset in BZ_list:
    search_terms = generate_search_terms(user_defined_Qpoints, BZ_offset)
    plot_workspaces = batch_search_and_load(folder, search_terms)

    print(len(plot_workspaces))
    plot_data = np.vstack([np.squeeze(data.getSignalArray()) for data in plot_workspaces])
    print(plot_data.shape)

    plt.figure(figsize=(4.25, 2.5))

    # Filter out NaNs and zeros
    valid_vals = plot_data[np.isfinite(plot_data) & (plot_data != 0)]

    # Compute percentiles
    vmin = np.percentile(valid_vals, 5)
    vmax = np.percentile(valid_vals, 98)

    # Create the plot with dynamic vmin and vmax
    cb = plt.pcolormesh(
        plot_data.T,
        shading='auto',
        cmap='viridis',
        norm=SymLogNorm(linthresh=5e-4, vmin=vmin, vmax=vmax)
    )

    # Set X-axis ticks and labels
    x_tick_positions = [0, 40, 50, 80, 100, 120, 140]
    x_tick_labels = [r"$\Gamma$", "X", "U,K", r"$\Gamma$", "L", "W", "X"]
    plt.xticks(ticks=x_tick_positions, labels=x_tick_labels, fontsize=12)

    # Set Y-axis ticks and labels
    y_tick_positions = [0, 20, 40, 60, 80]
    y_tick_labels = [0, 10, 20, 30, 40]
    plt.yticks(ticks=y_tick_positions, labels=y_tick_labels, fontsize=12)

    # Add colorbar
    cb = plt.colorbar(cb, shrink=0.8)
    cb.set_label('Intensity (a.u.)', fontsize=8)
    cb.ax.tick_params(labelsize=6)
    
    # Save figure
    plt.tight_layout()
    plt.savefig("temp_searching/Ge_path_BZ{}.png".format(BZ_offset), dpi=600, bbox_inches="tight")
    np.save('temp_searching/Ge_path_BZ{}.npy'.format(BZ_offset), plot_data)
    plt.close()