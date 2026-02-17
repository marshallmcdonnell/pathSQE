import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize, SymLogNorm
from matplotlib.cm import ScalarMappable, viridis
from matplotlib.backends.backend_pdf import PdfPages
from mantid.simpleapi import *

from . import simulations
from . import helper


def update_BZ_folding_plot(BZ_full_array, BZ_offset, pathSQE_params):
    """
    Updates and saves the BZ Coverage scatter plot during iterative folding.

    Parameters:
        BZ_full_list (list): List of (BZ_center, Fractional_Coverage) tuples.
        BZ_offset (array): The BZ center that was just processed.
        pathSQE_params (dict): Dictionary containing 'output directory' path.
        min_fraction (float): Minimum fraction threshold line.
    """
    # Find the index of the last processed BZ (BZ_offset) in the sorted array
    min_fraction = pathSQE_params["BZ to process"]
    BZ_centers = np.vstack(BZ_full_array[:, 0])  # Extract BZ center coordinates
    match_idx = np.where(np.all(BZ_centers == BZ_offset, axis=1))[0]

    if match_idx.size == 0:
        print(
            "Warning: BZ_offset not found in BZ_full_array. It will not be highlighted."
        )
        match_idx = None
    else:
        match_idx = match_idx[0]  # Convert to integer index

    # --- Scatter Plot: Coverage vs. Rank ---
    plt.figure(figsize=(6, 4))
    plt.scatter(
        range(len(BZ_full_array)),
        BZ_full_array[:, 1],
        color="tab:blue",
        s=10,
        label="BZ Coverage",
    )
    plt.axhline(
        y=min_fraction,
        color="tab:grey",
        linestyle="--",
        label=f"Threshold ({min_fraction:.2f})",
    )

    # Mark the last processed BZ with a red star
    if match_idx is not None:
        plt.scatter(
            match_idx,
            BZ_full_array[match_idx, 1],
            color="red",
            marker="*",
            s=30,
            label="Last Processed BZ",
        )

    # Labels and title
    plt.xlabel("BZ Rank")
    plt.ylabel("Fractional Coverage")
    plt.title("BZ Coverage vs. Rank")
    plt.legend()
    plt.grid(True)

    # Save the updated plot
    plt.savefig(
        os.path.join(pathSQE_params["output directory"], "BZ_coverage_iteration.png")
    )
    plt.close()


# List of known math symbols that are safe with \ in mathtext
known_mathtext_symbols = {
    "Gamma",
    "Delta",
    "Sigma",
    "Pi",
    "Lambda",
    "Omega",
    "Theta",
    "alpha",
    "beta",
    "gamma",
    "delta",
    "epsilon",
    "zeta",
    "eta",
    "theta",
    "kappa",
    "lambda",
    "mu",
    "nu",
    "xi",
    "pi",
    "rho",
    "sigma",
    "tau",
    "phi",
    "chi",
    "psi",
    "omega",
}


def safe_label(pt):
    if pt in known_mathtext_symbols:
        return rf"$\{pt}$"  # LaTeX mathtext symbol
    elif len(pt) == 1:
        return rf"${pt}$"  # Single character, keep in mathtext
    else:
        return pt  # Plain string


def set_seg_xlabels(seg_axis, i, dsl_fold):
    if i == len(dsl_fold) - 1:
        prev_end = dsl_fold[i - 1]["seg_end_name"]
        curr_start = dsl_fold[i]["seg_start_name"]
        curr_end = dsl_fold[i]["seg_end_name"]

        pt_labels = [safe_label(pt) for pt in [prev_end, curr_start, curr_end]]

        if curr_start == prev_end:
            seg_axis.set_xticklabels([pt_labels[1], pt_labels[2]])
        else:
            seg_axis.set_xticklabels([f"{pt_labels[0]},{pt_labels[1]}", pt_labels[2]])

    elif i == 0:
        curr_start = dsl_fold[i]["seg_start_name"]
        pt_label = safe_label(curr_start)
        seg_axis.set_xticklabels([pt_label, ""])

    else:
        prev_end = dsl_fold[i - 1]["seg_end_name"]
        curr_start = dsl_fold[i]["seg_start_name"]
        pt_labels = [safe_label(pt) for pt in [prev_end, curr_start]]

        if curr_start == prev_end:
            seg_axis.set_xticklabels([pt_labels[1], ""])
        else:
            seg_axis.set_xticklabels([f"{pt_labels[0]},{pt_labels[1]}", ""])


def plot_along_path_foldedBZ(foldedSegNames, dsl_fold, Ei, vmi, vma, cma="jet"):

    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["savefig.dpi"] = 150
    plt.rcParams.update({"font.size": 28})

    num_segments = len(foldedSegNames)
    ratios = [dsl_fold[i]["inv_angstrom_ratio"] for i in range(num_segments)]

    # Add space for the colorbar
    fig, axes = plt.subplots(
        1,
        num_segments,
        gridspec_kw={"width_ratios": ratios},
        figsize=(12, 5),
        subplot_kw={"projection": "mantid"},
    )  # constrained_layout=True

    colormesh_pars = {
        "norm": SymLogNorm(linthresh=vmi, vmin=vmi, vmax=vma),
        "cmap": cma,
    }

    for i in range(num_segments):
        x_start = dsl_fold[i]["qdim0_range"][0]
        x_end = dsl_fold[i]["qdim0_range"][2]

        if i == 0:
            ax1 = axes[i]
            ax1.pcolormesh(mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars)
            ax1.set_ylabel("E (meV)")
            ax1.set_xlabel("")
            # ax1.set_ylim(0., Ei)
            ax1.set_xticks([x_start, x_end])
            set_seg_xlabels(ax1, i, dsl_fold)
            ax1.tick_params(direction="in")
        elif i == num_segments - 1:
            axL = axes[i]
            mappable = axL.pcolormesh(
                mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars
            )
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel("")
            axL.set_xticks([x_start, x_end])
            set_seg_xlabels(axL, i, dsl_fold)
            axL.tick_params(direction="in")
        else:
            ax = axes[i]
            ax.pcolormesh(mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars)
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel("")
            ax.set_xticks([x_start, x_end])
            set_seg_xlabels(ax, i, dsl_fold)
            ax.tick_params(direction="in")

    plt.subplots_adjust(wspace=0, hspace=0)

    # Create the colorbar
    cb_ax = fig.add_axes([0.91, 0.124, 0.02, 0.754])
    cbar = fig.colorbar(mappable, orientation="vertical", cax=cb_ax)
    cbar.ax.tick_params(labelsize=15)

    plt.rcParams.update({"font.size": 10})

    return fig


def generate_BZ_Report_morePages(
    pathSQE_params, BZ_offset, all_sliceInfo_in_BZ, dsl_fold
):
    pdf_filename = os.path.join(
        pathSQE_params["output directory"], "Reports/BZ_Report_{}.pdf".format(BZ_offset)
    )
    with PdfPages(pdf_filename) as pdf:
        # Create a new figure for the final folded path in given BZ
        fig_path_foldedBZ, ax_path_foldedBZ = plt.subplots(
            subplot_kw={"projection": "mantid"}, figsize=(8, 6)
        )

        # plot the final folded path
        pathNames = [ds["Name"] for ds in dsl_fold]
        sortedWorkspaceNames = pathNames.copy()
        sortedWorkspaceNames = [name + ".nxs" for name in sortedWorkspaceNames]

        # Adaptive colorbars
        pathMin, pathMax = helper.adaptive_colorbars(
            sortedWorkspaceNames, pathSQE_params, percentile_min=10, percentile_max=97
        )

        fig_path_foldedBZ = plot_along_path_foldedBZ(
            foldedSegNames=sortedWorkspaceNames,
            dsl_fold=dsl_fold,
            Ei=pathSQE_params["T and Ei conditions"][0][1],
            vmi=pathMin,
            vma=pathMax,
            cma="viridis",
        )

        # Add BZ_offset information at the top of the final folded path figure
        fig_path_foldedBZ.suptitle(
            f"Folded data along path in BZ {BZ_offset}", fontsize=18
        )

        # Save the final folded path figure to the PDF
        pdf.savefig(fig_path_foldedBZ)
        plt.close(fig_path_foldedBZ)

        # now make page for each path seg showing all of the individual slices that were combined
        for i in range(len(all_sliceInfo_in_BZ)):  # range(len(path['path'])):
            path_seg = [dsl_fold[i]["seg_start_name"], dsl_fold[i]["seg_end_name"]]

            # Calculate the number of pages needed for this path segment
            num_pages = int(np.ceil(len(all_sliceInfo_in_BZ[i]) / 12))

            for page_num in range(num_pages):
                # Calculate the start and end index for slices on this page
                start_idx = page_num * 12
                end_idx = min((page_num + 1) * 12, len(all_sliceInfo_in_BZ[i]))

                # Create a new figure for each path_seg
                fig, axes = plt.subplots(
                    nrows=3,
                    ncols=4,
                    subplot_kw={"projection": "mantid"},
                    figsize=(12, 9),
                )
                plt.rcParams["figure.dpi"] = 150
                plt.rcParams["savefig.dpi"] = 150
                plt.set_cmap("viridis")

                for j, slice_info in enumerate(
                    all_sliceInfo_in_BZ[i][start_idx:end_idx]
                ):
                    slice_name = slice_info[0]
                    slice_ID = slice_info[1]
                    slice_eval = slice_info[2]

                    # Adjust the subplot position based on the number of rows and columns
                    row_index = j // 4
                    col_index = j % 4

                    step = float(pathSQE_params["E bins"].split(",")[1])
                    ind_aboveElasticLine = int(np.ceil(3 / step))
                    slice = mtd[slice_name].getSignalArray()[
                        :, 0, 0, ind_aboveElasticLine:
                    ]
                    # Mask NaN and zero values
                    mask = ~np.isnan(slice) & (slice > 0)

                    if np.any(mask != 0):
                        segMin = np.percentile(slice[mask], 10)
                        segMax = np.percentile(slice[mask], 97)
                        axes[row_index, col_index].pcolormesh(
                            mtd[slice_name], vmin=segMin, vmax=segMax, rasterized=True
                        )
                        axes[row_index, col_index].set_title("ID {}".format(slice_ID))
                        if not slice_eval:
                            axes[row_index, col_index].set_title(
                                "Flagged Slice (ID {})".format(slice_ID), color="red"
                            )

                    else:
                        axes[row_index, col_index].text(
                            0.5,
                            0.5,
                            "No data",
                            ha="center",
                            va="center",
                            fontsize=10,
                            color="red",
                        )

                    # x ticks and labels
                    point1_str, point2_str = slice_name.replace("pathSQE_", "").split(
                        "_to_"
                    )
                    point1 = (
                        np.array([float(x) for x in point1_str.split(",")]) + BZ_offset
                    )
                    point2 = (
                        np.array([float(x) for x in point2_str.split(",")]) + BZ_offset
                    )
                    # print(slice_name, point1_str, point2_str, point1, point2)
                    x_min_exp, x_max_exp = axes[row_index, col_index].get_xlim()
                    # print(x_min_exp, x_max_exp)
                    axes[row_index, col_index].set_xticks([x_min_exp, x_max_exp])
                    axes[row_index, col_index].set_xticklabels(
                        [f"{tick}" for tick in [point1, point2]]
                    )

                # Hide unused subplots (case 3: non-existent slices on the last page)
                total_plots = len(all_sliceInfo_in_BZ[i][start_idx:end_idx])
                n_rows, n_cols = axes.shape

                for idx in range(total_plots, n_rows * n_cols):
                    row_index = idx // n_cols
                    col_index = idx % n_cols
                    axes[row_index, col_index].set_visible(False)

                # Add BZ_offset and path_seg information at the top of each page
                fig.suptitle(
                    f"BZ: {BZ_offset}, Path Segment: {path_seg[0]} to {path_seg[1]} - Page {page_num + 1}",
                    y=0.98,
                    fontsize=14,
                )

                # Adjust layout to prevent overlapping labels
                fig.tight_layout(h_pad=0.3, w_pad=0.3)

                # Save the current figure to the PDF
                pdf.savefig(fig)
                plt.close("all")


def IqE_to_SED(SED_ws, T, Ebins):
    SED_ws_shape = SED_ws.getSignalArray().shape
    k = 8.617e-2
    # Split the string to extract min, step, and max values
    min_val, step_val, max_val = map(float, Ebins.split(","))

    # Generate the vector using numpy's arange function
    E_vals = np.arange(min_val + step_val / 2, max_val + step_val / 2, step_val)
    E_terms = np.exp(E_vals / (k * T))
    factors = E_vals / (1 + (1 / (E_terms - 1)))
    reshaped_factors = np.tile(
        factors.reshape(1, 1, 1, len(factors)),
        (SED_ws_shape[0], SED_ws_shape[1], SED_ws_shape[2], 1),
    )
    SED_ws.setSignalArray(SED_ws.getSignalArray() * reshaped_factors)

    return SED_ws


def plot_SED_along_path_foldedBZ(foldedSegNames, dsl_fold, pathSQE_params):
    T = pathSQE_params["T and Ei conditions"][0][0]
    Ebins = pathSQE_params["E bins"]
    pathSQE_params["T and Ei conditions"][0][1]
    cma = "viridis"

    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["savefig.dpi"] = 150
    plt.rcParams.update({"font.size": 28})

    fig = plt.figure(figsize=(12, 5))
    colormesh_pars = {}
    colormesh_pars["cmap"] = cma
    vmi = 0
    vma = 0

    num_segments = len(foldedSegNames)
    ratios = []
    for i in range(num_segments):
        ratios.append(dsl_fold[i]["inv_angstrom_ratio"])

    gs = GridSpec(1, num_segments, width_ratios=ratios, wspace=0)

    for i in range(num_segments):
        x_start = dsl_fold[i]["qdim0_range"][0]
        x_end = dsl_fold[i]["qdim0_range"][2]

        SED_ws = mtd[foldedSegNames[i]].clone()
        SED_ws = IqE_to_SED(SED_ws=SED_ws, T=T, Ebins=Ebins)
        SaveMD(
            SED_ws,
            Filename=os.path.join(
                pathSQE_params["output directory"],
                "slices/fully_folded",
                "{}2{}_folded_SED.nxs".format(
                    dsl_fold[i]["seg_start_name"], dsl_fold[i]["seg_end_name"]
                ),
            ),
        )

        if i == 0:
            slice = SED_ws.getSignalArray()
            mask = ~np.isnan(slice) & (slice > 0)
            vmi = np.percentile(slice[mask], 10)
            vma = np.percentile(slice[mask], 97)
            colormesh_pars["norm"] = SymLogNorm(
                linthresh=vmi, vmin=vmi, vmax=vma
            )  # (linthresh=9e-6,vmin=9e-6,vmax=2e-3)

            ax1 = plt.subplot(gs[i], projection="mantid")
            ax1.pcolormesh(SED_ws, **colormesh_pars)
            ax1.set_ylabel("E (meV)")
            ax1.set_xlabel("")
            # ax1.set_ylim(0.,Ei)
            ax1.set_xticks([x_start, x_end])
            set_seg_xlabels(ax1, i, dsl_fold)
            ax1.tick_params(direction="in")
        elif i == num_segments - 1:
            axL = plt.subplot(gs[i], sharey=ax1, projection="mantid")
            cb = axL.pcolormesh(SED_ws, **colormesh_pars)
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel("")
            axL.set_xticks([x_start, x_end])
            set_seg_xlabels(axL, i, dsl_fold)
            axL.tick_params(direction="in")
        else:
            ax = plt.subplot(gs[i], sharey=ax1, projection="mantid")
            ax.pcolormesh(SED_ws, **colormesh_pars)
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel("")
            ax.set_xticks([x_start, x_end])
            set_seg_xlabels(ax, i, dsl_fold)
            ax.tick_params(direction="in")

    fig.colorbar(cb)
    plt.rcParams.update({"font.size": 10})

    return fig


def generate_BZ_Report_1DSymPoints(
    pathSQE_params, symPoint, symPt_slice_arrays, symPt_slice_names_allTemp, u, v
):
    pdf_filename = os.path.join(
        pathSQE_params["output directory"],
        "Reports/{}point_Report.pdf".format(symPoint),
    )
    with PdfPages(pdf_filename) as pdf:
        fig_scatter, ax_scatter = plt.subplots(figsize=(8, 6))

        # Convert list to NumPy array for sorting
        symPt_slice_array = np.array(symPt_slice_arrays)

        # Compute out-of-plane unit normal
        u = u.astype(float) / np.linalg.norm(u)
        v = v.astype(float) / np.linalg.norm(v)
        normal_vec = np.cross(u, v)
        normal_vec /= np.linalg.norm(normal_vec)  # Normalize

        # Transform (H, K, L) points into (u, v) plane coordinates
        x = np.dot(symPt_slice_array[:, :3], u)  # Project onto u
        y = np.dot(symPt_slice_array[:, :3], v)  # Project onto v
        z = np.dot(symPt_slice_array[:, :3], normal_vec)  # Out-of-plane coordinate

        # Sort by out-of-plane coordinate
        sorted_indices = np.argsort(z)
        x, y, z = x[sorted_indices], y[sorted_indices], z[sorted_indices]

        # Normalize z values between 0 and 1
        norm = Normalize(vmin=np.min(z), vmax=np.max(z))
        norm_z = norm(z)

        # Define marker size based on z values (larger for closer points)
        marker_size = 500 * (1 - norm_z) ** 2 + 30

        # Define colors based on z values
        colors = viridis(norm_z)

        # Plot scatter plot
        ax_scatter.scatter(x, y, c=colors, s=marker_size, alpha=1)

        # Add color bar
        sm = ScalarMappable(norm=norm, cmap=viridis)
        sm.set_array([])
        plt.colorbar(sm, ax=ax_scatter, label="Out-of-Plane Coord")

        # Set grid lines aligned with u and v
        x_half_integer_locs = np.arange(
            np.floor(np.min(x)) - 0.5, np.ceil(np.max(x)) + 1.5, 1
        )
        y_half_integer_locs = np.arange(
            np.floor(np.min(y)) - 0.5, np.ceil(np.max(y)) + 1.5, 1
        )
        ax_scatter.xaxis.set_major_locator(plt.FixedLocator(x_half_integer_locs))
        ax_scatter.yaxis.set_major_locator(plt.FixedLocator(y_half_integer_locs))
        ax_scatter.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax_scatter.yaxis.set_minor_locator(plt.MultipleLocator(1))
        ax_scatter.xaxis.set_minor_formatter(
            plt.FuncFormatter(lambda x, _: f"{int(x)}")
        )
        ax_scatter.yaxis.set_minor_formatter(
            plt.FuncFormatter(lambda y, _: f"{int(y)}")
        )
        ax_scatter.grid(
            which="major", linestyle="-", linewidth=0.5, color="gray", alpha=0.3
        )

        # Set axis labels based on u and v vectors
        ax_scatter.set_xlabel(f"u = {u}")
        ax_scatter.set_ylabel(f"v = {v}")
        fig_scatter.suptitle(
            "{} {} points with data".format(len(symPt_slice_arrays), symPoint),
            y=0.95,
            fontsize=18,
        )

        # Save figure to PDF
        pdf.savefig(fig_scatter)
        plt.close(fig_scatter)

        # Get e range vector for plotting
        spectrum_labels = pathSQE_params["T and Ei conditions"]
        start, step, end = map(float, pathSQE_params["E bins"].split(","))
        evec_tot = np.arange(start + 0.5 * step, end - 0.5 * step + step, step)

        # Calculate number of pages needed
        num_pages = int(np.ceil(len(symPt_slice_arrays) / 12))

        for page_num in range(num_pages):
            start_idx = page_num * 12
            end_idx = min((page_num + 1) * 12, len(symPt_slice_arrays))
            symPts_subset = symPt_slice_arrays[start_idx:end_idx]

            # Create figure for each page
            fig, axes = plt.subplots(
                nrows=3, ncols=4, subplot_kw={"projection": "mantid"}, figsize=(12, 9)
            )
            plt.rcParams["figure.dpi"] = 150
            plt.rcParams["savefig.dpi"] = 150

            for j in range(len(symPts_subset)):
                row_index = j // 4
                col_index = j % 4

                # Get colormap
                cmap = plt.get_cmap("coolwarm")

                for k in range(len(symPt_slice_names_allTemp)):
                    slice_names_subset = symPt_slice_names_allTemp[k][start_idx:end_idx]
                    slice_name = slice_names_subset[j]

                    color = None
                    if len(symPt_slice_names_allTemp) > 1:
                        norm_k = k / (len(symPt_slice_names_allTemp) - 1)
                        color = cmap(norm_k)
                    else:
                        color = "tab:green"

                    if k == 0:
                        # First temperature: Normalize using percentile
                        signal_tot = mtd[slice_name].getSignalArray()[0, 0, 0, :]
                        nonNaN_indices = np.where(~np.isnan(signal_tot))[0]
                        signal = signal_tot[nonNaN_indices]

                        vert_offset = 1.1 * np.percentile(
                            signal, 98
                        )  # Use 95th percentile instead of max
                        peak_scaling = np.percentile(
                            signal, 98
                        )  # Scale sim by 95th percentile

                        signal += vert_offset
                        error = np.sqrt(
                            mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :]
                        )
                        error = error[nonNaN_indices]
                        e_vec = evec_tot[nonNaN_indices]

                        axes[row_index, col_index].errorbar(
                            x=e_vec, y=signal, yerr=error, color=color
                        )
                        axes[row_index, col_index].axhline(
                            y=vert_offset, color="black", alpha=0.15
                        )

                        # **Only perform simulation if pathSQE_params['perform simulations'] is True**
                        if pathSQE_params["perform simulations"]:
                            # print('T? ', pathSQE_params['T and Ei conditions'][k][0])
                            # print(symPts_subset[j])
                            sim_SQE_output = simulations.sim_SQE(
                                pathSQE_params,
                                [
                                    np.linalg.inv(
                                        pathSQE_params[
                                            "primitive to mantid transformation matrix"
                                        ]
                                    )
                                    @ np.array(symPts_subset[j])
                                ],
                                Temperature=pathSQE_params["T and Ei conditions"][k][0],
                            )
                            unique_freq, summed_intensities, freq_range, spectrum = (
                                simulations.SQE_to_1d_spectrum(
                                    pathSQE_params, sim_SQE_output
                                )
                            )

                            sim_mask = (unique_freq > np.min(e_vec)) & (
                                unique_freq < np.max(e_vec)
                            )
                            unique_freq = unique_freq[sim_mask]
                            summed_intensities = summed_intensities[sim_mask]

                            sim_mask = (freq_range > np.min(e_vec)) & (
                                freq_range < np.max(e_vec)
                            )
                            freq_range = freq_range[sim_mask]
                            spectrum = spectrum[sim_mask]

                            summed_intensities *= peak_scaling / np.max(spectrum)
                            spectrum *= peak_scaling / np.max(spectrum)

                            axes[row_index, col_index].scatter(
                                unique_freq, summed_intensities, color="black", s=5
                            )
                            axes[row_index, col_index].plot(
                                freq_range, spectrum, color="black"
                            )

                            # if symPoint == 'R':
                            #    np.save(os.path.join(pathSQE_params['output directory'],'nxs_files/',slice_name,'_simSpectrum.npy'), np.vstack((freq_range, spectrum)))
                            #    np.save(os.path.join(pathSQE_params['output directory'],'nxs_files/',slice_name,'_simScatter.npy'), np.vstack((unique_freq, summed_intensities)))

                    else:
                        # Other temperatures
                        signal_tot2 = mtd[slice_name].getSignalArray()[0, 0, 0, :]
                        nonNaN_indices2 = np.where(~np.isnan(signal_tot2))[0]
                        signal2 = signal_tot2[nonNaN_indices2]
                        signal2 += (k + 1) * vert_offset

                        error2 = np.sqrt(
                            mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :]
                        )
                        error2 = error2[nonNaN_indices2]
                        e_vec2 = evec_tot[nonNaN_indices2]

                        axes[row_index, col_index].errorbar(
                            x=e_vec2, y=signal2, yerr=error2, color=color
                        )
                        axes[row_index, col_index].axhline(
                            y=(k + 1) * vert_offset, color="black", alpha=0.15
                        )

                        if k == len(symPt_slice_names_allTemp) - 1:
                            axes[row_index, col_index].set_ylim(
                                [0, (k + 3) * vert_offset]
                            )

                    axes[row_index, col_index].set_xlim([start, end])
                    axes[row_index, col_index].set_xlabel("Energy (meV)")
                    axes[row_index, col_index].set_title("{}".format(symPts_subset[j]))

                # **Add a legend to the top-left subplot**
                if (
                    j == 0 and len(symPt_slice_names_allTemp) > 1
                ):  # Top-left plot on each page
                    legend_handles = [
                        plt.Line2D(
                            [0], [0], color=cmap(i / (len(spectrum_labels) - 1)), lw=3
                        )
                        for i in range(len(spectrum_labels))
                    ]
                    axes[row_index, col_index].legend(
                        legend_handles, spectrum_labels, loc="upper right", fontsize=8
                    )

            # Adjust layout and save the PDF page
            fig.tight_layout(h_pad=0.3, w_pad=0.3)
            pdf.savefig(fig)
            plt.close("all")


def plot_simple_1d_pt(pathSQE_params, symPoint, symPt_init, slice_names):
    # Get e range vector for plotting
    spectrum_labels = pathSQE_params["T and Ei conditions"]
    start, step, end = map(float, pathSQE_params["E bins"].split(","))
    evec_tot = np.arange(start + 0.5 * step, end - 0.5 * step + step, step)

    # Create figure for each page
    fig, axes = plt.subplots(subplot_kw={"projection": "mantid"})
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["savefig.dpi"] = 150

    # Get colormap
    cmap = plt.get_cmap("coolwarm")

    for k in range(len(pathSQE_params["T and Ei conditions"])):
        slice_name = slice_names[k]

        color = None
        if len(pathSQE_params["T and Ei conditions"]) > 1:
            norm_k = k / (len(pathSQE_params["T and Ei conditions"]) - 1)
            color = cmap(norm_k)
        else:
            color = "tab:green"

        if k == 0:
            # First temperature: Normalize using percentile
            signal_tot = mtd[slice_name].getSignalArray()[0, 0, 0, :]
            nonNaN_indices = np.where(~np.isnan(signal_tot))[0]
            signal = signal_tot[nonNaN_indices]

            vert_offset = 1.1 * np.percentile(
                signal, 98
            )  # Use 95th percentile instead of max
            peak_scaling = np.percentile(signal, 98)  # Scale sim by 95th percentile

            signal += vert_offset
            error = np.sqrt(mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :])
            error = error[nonNaN_indices]
            e_vec = evec_tot[nonNaN_indices]

            axes.errorbar(x=e_vec, y=signal, yerr=error, color=color)
            axes.axhline(y=vert_offset, color="black", alpha=0.15)

            # **Only perform simulation if pathSQE_params['perform simulations'] is True**
            if pathSQE_params["perform simulations"]:
                # print('T? ', pathSQE_params['T and Ei conditions'][k][0])
                # print(symPts_subset[j])
                sim_SQE_output = simulations.sim_SQE(
                    pathSQE_params,
                    [
                        np.linalg.inv(
                            pathSQE_params["primitive to mantid transformation matrix"]
                        )
                        @ symPt_init
                    ],
                    Temperature=pathSQE_params["T and Ei conditions"][k][0],
                )
                unique_freq, summed_intensities, freq_range, spectrum = (
                    simulations.SQE_to_1d_spectrum(pathSQE_params, sim_SQE_output)
                )

                sim_mask = (unique_freq > np.min(e_vec)) & (unique_freq < np.max(e_vec))
                unique_freq = unique_freq[sim_mask]
                summed_intensities = summed_intensities[sim_mask]

                sim_mask = (freq_range > np.min(e_vec)) & (freq_range < np.max(e_vec))
                freq_range = freq_range[sim_mask]
                spectrum = spectrum[sim_mask]

                summed_intensities *= peak_scaling / np.max(spectrum)
                spectrum *= peak_scaling / np.max(spectrum)

                axes.scatter(unique_freq, summed_intensities, color="black", s=5)
                axes.plot(freq_range, spectrum, color="black")

                if pathSQE_params["save individual slices"]:
                    np.save(
                        os.path.join(
                            pathSQE_params["output directory"],
                            "slices/all_individual/",
                            slice_name + "_simSpectrum.npy",
                        ),
                        np.vstack((freq_range, spectrum)),
                    )
                    np.save(
                        os.path.join(
                            pathSQE_params["output directory"],
                            "slices/all_individual/",
                            slice_name + "_simScatter.npy",
                        ),
                        np.vstack((unique_freq, summed_intensities)),
                    )

        else:
            # Other temperatures
            signal_tot2 = mtd[slice_name].getSignalArray()[0, 0, 0, :]
            nonNaN_indices2 = np.where(~np.isnan(signal_tot2))[0]
            signal2 = signal_tot2[nonNaN_indices2]
            signal2 += (k + 1) * vert_offset

            error2 = np.sqrt(mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :])
            error2 = error2[nonNaN_indices2]
            e_vec2 = evec_tot[nonNaN_indices2]

            axes.errorbar(x=e_vec2, y=signal2, yerr=error2, color=color)
            axes.axhline(y=(k + 1) * vert_offset, color="black", alpha=0.15)

            if k == len(pathSQE_params["T and Ei conditions"]) - 1:
                axes.set_ylim([0, (k + 3) * vert_offset])

        axes.set_xlim([start, end])
        axes.set_xlabel("Energy (meV)")
        axes.set_title("{}".format(symPt_init))

    # **Add a legend if applicable**
    if len(pathSQE_params["T and Ei conditions"]) > 1:
        legend_handles = [
            plt.Line2D([0], [0], color=cmap(i / (len(spectrum_labels) - 1)), lw=3)
            for i in range(len(spectrum_labels))
        ]
        axes.legend(legend_handles, spectrum_labels, loc="upper right", fontsize=8)

    # Adjust layout and save the PDF page
    fig.savefig(
        os.path.join(
            pathSQE_params["output directory"],
            "{}point_{}.png".format(symPoint, symPt_init),
        )
    )
    plt.close("all")


def plot_along_path_foldedBZ_sim(sim_folded_data, dsl_fold, Ei, vmi, vma, cma="jet"):

    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["savefig.dpi"] = 150
    plt.rcParams.update({"font.size": 28})

    num_segments = len(dsl_fold)
    ratios = [dsl_fold[i]["inv_angstrom_ratio"] for i in range(num_segments)]

    # Extract binning info and generate the correct energy values
    E_min, E_step, E_max = map(float, dsl_fold[0]["Dimension3Binning"].split(","))
    num_y_bins = sim_folded_data[0].shape[1]  # Assuming (x_bins, y_bins)
    E_values = np.linspace(
        E_min, E_max, num_y_bins + 1
    )  # Use num_y_bins + 1 for bin edges

    # Add space for the colorbar
    fig, axes = plt.subplots(
        1,
        num_segments,
        gridspec_kw={"width_ratios": ratios},
        figsize=(12, 5),
        subplot_kw={"projection": "mantid"},
    )  # constrained_layout=True

    colormesh_pars = {
        "norm": SymLogNorm(linthresh=vmi, vmin=vmi, vmax=vma),
        "cmap": cma,
    }

    for i in range(num_segments):
        x_start = 0
        x_end = sim_folded_data[i].shape[0]

        if i == 0:
            ax1 = axes[i]
            ax1.pcolormesh(
                np.arange(sim_folded_data[i].shape[0] + 1),
                E_values,
                sim_folded_data[i].T,
                rasterized=True,
                **colormesh_pars,
            )
            ax1.set_ylabel("E (meV)")
            ax1.set_xlabel("")
            # ax1.set_ylim(0., Ei)
            ax1.set_xticks([x_start, x_end])
            set_seg_xlabels(ax1, i, dsl_fold)
            ax1.tick_params(direction="in")
        elif i == num_segments - 1:
            axL = axes[i]
            mappable = axL.pcolormesh(
                np.arange(sim_folded_data[i].shape[0] + 1),
                E_values,
                sim_folded_data[i].T,
                rasterized=True,
                **colormesh_pars,
            )
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel("")
            axL.set_xticks([x_start, x_end])
            set_seg_xlabels(axL, i, dsl_fold)
            axL.tick_params(direction="in")
        else:
            ax = axes[i]
            ax.pcolormesh(
                np.arange(sim_folded_data[i].shape[0] + 1),
                E_values,
                sim_folded_data[i].T,
                rasterized=True,
                **colormesh_pars,
            )
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel("")
            ax.set_xticks([x_start, x_end])
            set_seg_xlabels(ax, i, dsl_fold)
            ax.tick_params(direction="in")

    plt.subplots_adjust(wspace=0, hspace=0)

    # Create the colorbar
    cb_ax = fig.add_axes([0.91, 0.124, 0.02, 0.754])
    cbar = fig.colorbar(mappable, orientation="vertical", cax=cb_ax)
    cbar.ax.tick_params(labelsize=15)

    plt.rcParams.update({"font.size": 10})

    return fig


def generate_BZ_Report_with_Sims(
    pathSQE_params,
    BZ_offset,
    all_sliceInfo_in_BZ,
    dsl_fold,
    folded_results,
    all_sims_in_BZ,
):
    pdf_filename = os.path.join(
        pathSQE_params["output directory"], f"Reports/BZ_Report_{BZ_offset}.pdf"
    )

    with PdfPages(pdf_filename) as pdf:
        # Create a new figure for the final folded path in given BZ
        fig_path_foldedBZ, ax_path_foldedBZ = plt.subplots(
            subplot_kw={"projection": "mantid"}, figsize=(8, 6)
        )

        pathNames = [ds["Name"] for ds in dsl_fold]
        sortedWorkspaceNames = [name + ".nxs" for name in pathNames]

        # Adaptive colorbars
        pathMin, pathMax = helper.adaptive_colorbars(
            sortedWorkspaceNames, pathSQE_params, percentile_min=10, percentile_max=97
        )

        # Plot experimental folded data
        fig_path_foldedBZ = plot_along_path_foldedBZ(
            foldedSegNames=sortedWorkspaceNames,
            dsl_fold=dsl_fold,
            Ei=pathSQE_params["T and Ei conditions"][0][1],
            vmi=pathMin,
            vma=pathMax,
            cma="viridis",
        )

        # Add BZ_offset information at the top of the final folded path figure
        fig_path_foldedBZ.suptitle(
            f"Folded experimental data along path in BZ {BZ_offset}", fontsize=18
        )

        # Save the final folded path figure to the PDF
        pdf.savefig(fig_path_foldedBZ)
        plt.close(fig_path_foldedBZ)

        # Create a new figure for the final folded path in given BZ
        fig_path_foldedBZ2, ax_path_foldedBZ2 = plt.subplots(
            subplot_kw={"projection": "mantid"}, figsize=(8, 6)
        )

        # Plot simulated folded data
        sim_folded_data = [
            result[2] for result in folded_results
        ]  # Extract final folded sim data
        simMax = 0
        simMin = 1e4
        # Find min/max percentiles across all slices
        for path_sim in sim_folded_data:
            mask = ~np.isnan(path_sim) & (path_sim > 0)
            if np.any(mask):
                segMin = np.percentile(path_sim[mask], 60)
                segMax = np.max(path_sim[mask])
                simMax = max(pathMax, segMax)
                simMin = min(pathMin, segMin)

        fig_path_foldedBZ2 = plot_along_path_foldedBZ_sim(
            sim_folded_data,
            dsl_fold,
            pathSQE_params["T and Ei conditions"][0][1],
            simMin,
            simMax,
            "viridis",
        )

        # Add BZ_offset information at the top of the final folded path figure
        fig_path_foldedBZ2.suptitle(
            f"Folded simulated data along path in BZ {BZ_offset}", fontsize=18
        )

        # Save the final folded path figure to the PDF
        pdf.savefig(fig_path_foldedBZ2)
        plt.close(fig_path_foldedBZ2)

        # Loop over path segments and create slice comparison pages
        for i in range(len(all_sliceInfo_in_BZ)):
            path_seg = [dsl_fold[i]["seg_start_name"], dsl_fold[i]["seg_end_name"]]
            num_pages = int(
                np.ceil(len(all_sliceInfo_in_BZ[i]) / 6)
            )  # 6 pairs (12 subplots) per page

            for page_num in range(num_pages):
                start_idx = page_num * 6
                end_idx = min((page_num + 1) * 6, len(all_sliceInfo_in_BZ[i]))

                fig, axes = plt.subplots(
                    nrows=3,
                    ncols=4,
                    subplot_kw={"projection": "mantid"},
                    figsize=(12, 9),
                )
                plt.rcParams["figure.dpi"] = 150
                plt.rcParams["savefig.dpi"] = 150
                plt.set_cmap("viridis")

                for j, slice_idx in enumerate(range(start_idx, end_idx)):
                    slice_info = all_sliceInfo_in_BZ[i][slice_idx]
                    slice_name = slice_info[0]
                    slice_ID = slice_info[1]
                    slice_eval = slice_info[2]

                    if not slice_eval:
                        print(
                            "bad slice inside plotting ",
                            slice_name,
                            slice_ID,
                            slice_eval,
                        )

                    sim_data = all_sims_in_BZ[i][
                        slice_idx
                    ]  # Extract final folded simulation data
                    row_index = j // 2  # 3 rows
                    col_index = (j % 2) * 2  # Exp goes in 0, 2; Sim goes in 1, 3

                    # --- Experimental Slice ---
                    step = float(pathSQE_params["E bins"].split(",")[1])
                    ind_aboveElasticLine = int(np.ceil(3 / step))
                    slice = mtd[slice_name].getSignalArray()[
                        :, 0, 0, ind_aboveElasticLine:
                    ]
                    mask = ~np.isnan(slice) & (slice > 0)

                    if np.any(mask):
                        segMin = np.percentile(slice[mask], 10)
                        segMax = np.percentile(slice[mask], 97)
                        axes[row_index, col_index].pcolormesh(
                            mtd[slice_name], vmin=segMin, vmax=segMax, rasterized=True
                        )
                        axes[row_index, col_index].set_title("ID {}".format(slice_ID))
                        if not slice_eval:
                            axes[row_index, col_index].set_title(
                                "Flagged Slice (ID {})".format(slice_ID), color="red"
                            )

                        # x ticks and labels
                        point1_str, point2_str = slice_name.replace(
                            "pathSQE_", ""
                        ).split("_to_")
                        point1 = (
                            np.array([float(x) for x in point1_str.split(",")])
                            + BZ_offset
                        )
                        point2 = (
                            np.array([float(x) for x in point2_str.split(",")])
                            + BZ_offset
                        )
                        # print(slice_name, point1_str, point2_str, point1, point2)
                        x_min_exp, x_max_exp = axes[row_index, col_index].get_xlim()
                        # print(x_min_exp, x_max_exp)
                        axes[row_index, col_index].set_xticks([x_min_exp, x_max_exp])
                        axes[row_index, col_index].set_xticklabels(
                            [f"{tick}" for tick in [point1, point2]]
                        )

                        # y ticks and labels
                        exp_yticks = axes[row_index, col_index].get_yticks()
                        exp_yticklabels = [
                            tick.get_text()
                            for tick in axes[row_index, col_index].get_yticklabels()
                        ]
                        y_min_real, y_max_real = axes[row_index, col_index].get_ylim()
                        _, num_y_bins = sim_data.shape
                        sim_yticks = np.interp(
                            exp_yticks, [y_min_real, y_max_real], [0, num_y_bins - 1]
                        )
                        sim_yticks = np.round(sim_yticks).astype(int)
                        sim_yticks = np.clip(sim_yticks, 0, num_y_bins - 1)
                        axes[row_index, col_index + 1].set_yticks(sim_yticks)
                        axes[row_index, col_index + 1].set_yticklabels(
                            [
                                f"{float(tick.strip().replace('−', '-')):.2g}"
                                for tick in exp_yticklabels
                            ]
                        )

                        # Copy axis labels
                        axes[row_index, col_index + 1].set_xlabel(
                            axes[row_index, col_index].get_xlabel()
                        )
                        axes[row_index, col_index + 1].set_ylabel(
                            axes[row_index, col_index].get_ylabel()
                        )

                        # --- Simulated Slice (Right of Exp) ---
                        mask = ~np.isnan(sim_data) & (sim_data > 0)

                        if np.any(mask):
                            segMin = np.percentile(sim_data[mask], 60)
                            segMax = np.percentile(sim_data[mask], 99)

                        if sim_data is not None:
                            axes[row_index, col_index + 1].pcolormesh(
                                sim_data.T, vmin=segMin, vmax=segMax, rasterized=True
                            )
                            x_min_sim, x_max_sim = axes[
                                row_index, col_index + 1
                            ].get_xlim()
                            # print(x_min_sim, x_max_sim)
                            axes[row_index, col_index + 1].set_xticks(
                                [x_min_sim, x_max_sim]
                            )
                            axes[row_index, col_index + 1].set_xticklabels(
                                [f"{tick}" for tick in [point1, point2]]
                            )
                    else:
                        axes[row_index, col_index].text(
                            0.5,
                            0.5,
                            "No data",
                            ha="center",
                            va="center",
                            fontsize=10,
                            color="red",
                        )
                        axes[row_index, col_index + 1].text(
                            0.5,
                            0.5,
                            "No data",
                            ha="center",
                            va="center",
                            fontsize=10,
                            color="red",
                        )

                fig.suptitle(
                    f"BZ: {BZ_offset}, Path Segment: {path_seg[0]} to {path_seg[1]} - Page {page_num + 1}",
                    fontsize=14,
                )
                fig.tight_layout(h_pad=0.3, w_pad=0.3)

                pdf.savefig(fig)
                plt.close(fig)


def set_seg_xlabels_rerun(seg_axis, full_path, i):
    if i == len(full_path) - 1:
        prev_end = full_path[i - 1][1]
        curr_start = full_path[i][0]
        curr_end = full_path[i][1]

        pt_names = [prev_end, curr_start, curr_end]
        labels = [pt if len(pt) == 1 else "\\" + pt for pt in pt_names]

        if curr_start == prev_end:
            seg_axis.set_xticklabels(
                [r"${}$".format(labels[1]), r"${}$".format(labels[2])]
            )
        else:
            seg_axis.set_xticklabels(
                [r"${},{}$".format(labels[0], labels[1]), r"${}$".format(labels[2])]
            )

    elif i == 0:
        curr_start = full_path[i][0]

        pt_names = [curr_start]
        labels = [pt if len(pt) == 1 else "\\" + pt for pt in pt_names]
        seg_axis.set_xticklabels([r"${}$".format(labels[0]), ""])

    else:
        prev_end = full_path[i - 1][1]
        curr_start = full_path[i][0]

        pt_names = [prev_end, curr_start]
        labels = [pt if len(pt) == 1 else "\\" + pt for pt in pt_names]

        if curr_start == prev_end:
            seg_axis.set_xticklabels([r"${}$".format(labels[1]), ""])
        else:
            seg_axis.set_xticklabels([r"${},{}$".format(labels[0], labels[1]), ""])


def plot_along_path_rerun(
    foldedSegNames, user_defined_Qpoints, vmi, vma, cma="viridis"
):

    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["savefig.dpi"] = 150
    plt.rcParams.update({"font.size": 28})

    full_path = user_defined_Qpoints["path"]
    num_segments = len(full_path)

    ratios = []
    for k in range(len(full_path)):
        pt1 = np.array(user_defined_Qpoints["point_coords"][full_path[k][0]])
        pt2 = np.array(user_defined_Qpoints["point_coords"][full_path[k][1]])
        ratios.append(np.linalg.norm(pt2 - pt1))

    # Add space for the colorbar
    fig, axes = plt.subplots(
        1,
        num_segments,
        gridspec_kw={"width_ratios": ratios},
        figsize=(12, 5),
        subplot_kw={"projection": "mantid"},
    )  # constrained_layout=True

    colormesh_pars = {
        "norm": SymLogNorm(linthresh=vmi, vmin=vmi, vmax=vma),
        "cmap": cma,
    }

    for i in range(num_segments):
        x_start = mtd[foldedSegNames[i]].getDimension(0).getMinimum()
        x_end = mtd[foldedSegNames[i]].getDimension(0).getMaximum()

        if i == 0:
            ax1 = axes[i]
            ax1.pcolormesh(mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars)
            ax1.set_ylabel("E (meV)")
            ax1.set_xlabel("")
            # ax1.set_ylim(0., Ei)
            ax1.set_xticks([x_start, x_end])
            set_seg_xlabels_rerun(ax1, full_path, i)
            ax1.tick_params(direction="in")
        elif i == num_segments - 1:
            axL = axes[i]
            mappable = axL.pcolormesh(
                mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars
            )
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel("")
            axL.set_xticks([x_start, x_end])
            set_seg_xlabels_rerun(axL, full_path, i)
            axL.tick_params(direction="in")
        else:
            ax = axes[i]
            ax.pcolormesh(mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars)
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel("")
            ax.set_xticks([x_start, x_end])
            set_seg_xlabels_rerun(ax, full_path, i)
            ax.tick_params(direction="in")

    plt.subplots_adjust(wspace=0, hspace=0)

    # Create the colorbar
    cb_ax = fig.add_axes([0.91, 0.124, 0.02, 0.754])
    cbar = fig.colorbar(mappable, orientation="vertical", cax=cb_ax)
    cbar.ax.tick_params(labelsize=15)

    plt.rcParams.update({"font.size": 10})

    return fig
