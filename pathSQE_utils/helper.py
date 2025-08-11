import numpy as np
import os
from mantid.simpleapi import *
import time
import re
from seekpath.util import atoms_num_dict
from mantid.geometry import SpaceGroupFactory
from . import core

def simple_read_poscar(fname):
    """Read a POSCAR file."""
    with open(fname) as f:
        lines = [l.partition("!")[0] for l in f.readlines()]

    alat = float(lines[1])
    v1 = [float(_) * alat for _ in lines[2].split()]
    v2 = [float(_) * alat for _ in lines[3].split()]
    v3 = [float(_) * alat for _ in lines[4].split()]
    cell = [v1, v2, v3]

    species = lines[5].split()
    num_atoms = [int(_) for _ in lines[6].split()]

    next_line = lines[7]
    if next_line.strip().lower() != "direct":
        raise ValueError("This simple routine can only deal with 'direct' POSCARs")
    # Note: to support also cartesian, remember to multiply the coordinates by alat

    positions = []
    atomic_numbers = []
    cnt = 8
    for el, num in zip(species, num_atoms):
        atom_num = atoms_num_dict[el.capitalize()]
        for _ in range(num):
            atomic_numbers.append(atom_num)
            positions.append([float(_) for _ in lines[cnt].split()])
            cnt += 1

    return (cell, positions, atomic_numbers)



def conv_to_desc_string(input_item):
    # works for np arrays, lists, and tuples
    string = ''
    for item in input_item:
        string = string + str(item) +','
    string = string[:-1]
    
    return string



def make_slice_desc(pathSQE_params, q_dims_and_bins, point1, point2, path_seg):
    diff = point2 - point1

    # slice description with desired dims and bins and unique name
    slice_desc={'QDimension0':conv_to_desc_string(q_dims_and_bins[0]),
                'QDimension1':conv_to_desc_string(q_dims_and_bins[1]),
                'QDimension2':conv_to_desc_string(q_dims_and_bins[2]),
                'Dimension0Name':'QDimension0',
                'Dimension0Binning':conv_to_desc_string(q_dims_and_bins[3]),
                'Dimension1Name':'QDimension1',
                'Dimension1Binning':conv_to_desc_string(q_dims_and_bins[4]),
                'Dimension2Name':'QDimension2',
                'Dimension2Binning':conv_to_desc_string(q_dims_and_bins[5]),
                'Dimension3Name':'DeltaE',
                'Dimension3Binning':pathSQE_params['E bins'],
                'SymmetryOperations':'x,y,z',
                'Name':'pathSQE_'+conv_to_desc_string(point1)+'_to_'+conv_to_desc_string(point2),
                'qdim0_range':q_dims_and_bins[3],
                'seg_start_name':path_seg[0],
                'seg_end_name':path_seg[1],
                'inv_angstrom_ratio':np.linalg.norm(diff)/0.5,
                'good_slice':True}

    #slice_desc['Name'] = slice_desc['Name'].replace(',', '_').replace('.', '').replace('-','m')
    
    return slice_desc



def find_matching_spacegroup(pathSQE_params):
    query_spacegroup = pathSQE_params['space group']
    print(f'Searching for spacegroup {query_spacegroup}')
    spacegroup_list = SpaceGroupFactory.getAllSpaceGroupSymbols()
    query_spacegroup = ''.join(query_spacegroup.split())  # Remove spacing from the query_spacegroup

    for spacegroup in spacegroup_list:
        stripped_spacegroup = ''.join(spacegroup.split())  # Remove spacing from the spacegroup in the list
        if query_spacegroup == stripped_spacegroup:
            print(f"Matching spacegroup found: {spacegroup}")
            print(SpaceGroupFactory.createSpaceGroup(spacegroup))
            print(SpaceGroupFactory.createSpaceGroup(spacegroup).getPointGroup())
            return spacegroup

    raise ValueError("No matching spacegroup found.")



def get_next_sliceID(all_slices_dir):
    """
    Searches the given directory for files matching the slice pattern and 
    returns the next available unique slice ID (max existing + 1).
    """
    slice_pattern = re.compile(r'_ID(\d+)_data\.nxs$')
    existing_ids = []

    if not os.path.isdir(all_slices_dir):
        os.makedirs(all_slices_dir)
        return 0

    for fname in os.listdir(all_slices_dir):
        match = slice_pattern.search(fname)
        if match:
            existing_ids.append(int(match.group(1)))

    if existing_ids:
        return max(existing_ids) + 1
    else:
        return 0



def project_to_uv_plane(BZ_list, u, v):
    """
    Projects BZ centers onto the (u, v) plane.
    
    Parameters:
        BZ_list: List of tuples, where each tuple contains (BZ center, coverage)
        u: First basis vector defining the plane
        v: Second basis vector defining the plane
        
    Returns:
        np.array of projected (u, v) coordinates.
    """
    u = np.array(u) / np.linalg.norm(u)  # Normalize u
    v = np.array(v) / np.linalg.norm(v)  # Normalize v

    uv_matrix = np.vstack([u, v]).T  # Form transformation matrix
    uv_inv = np.linalg.pinv(uv_matrix)  # Compute pseudo-inverse for projection

    BZ_uv = np.array([uv_inv @ bz[0] for bz in BZ_list])  # Project each BZ center
    return BZ_uv



def get_existing_BZ_counts(output_folder_path):
    """
    Retrieve BZ arrays that already have processed .nxs files and count occurrences.

    Parameters:
        output_folder_path (str): Path to the folder containing processed BZ files.

    Returns:
        dict: Mapping of tuple(BZ center) -> count of occurrences.
    """
    BZ_counts = {}
    pattern = re.compile(r'_BZ\[(.*?)\]_final\.nxs$')  # Regex to capture the BZ array

    for filename in os.listdir(output_folder_path):
        match = pattern.search(filename)
        if match:
            BZ_str = match.group(1)  # Extract array content inside brackets
            try:
                BZ_array = tuple(float(x) for x in BZ_str.split())  # Convert to tuple for dict keys
                BZ_counts[BZ_array] = BZ_counts.get(BZ_array, 0) + 1
            except ValueError:
                print(f"Warning: Could not parse BZ array from filename '{filename}'")

    return BZ_counts



def resume_analysis(BZ_list, output_folder_path, num_paths, all_slices_dir):
    """
    Resume analysis by removing already processed BZ arrays from BZ_list,
    and delete any files in all_slices_dir not matching fully processed BZs.

    Parameters:
        BZ_list (np.ndarray): List of BZ centers to process.
        output_folder_path (str): Path to the folder containing processed BZ files.
        num_paths (int): Number of paths each BZ should have to be considered complete.
        all_slices_dir (str): Directory where individual slice files are saved.

    Returns:
        new_BZ_list (np.ndarray): Filtered list of BZs still needing processing.
        removed_BZ (np.ndarray): List of BZs that were fully processed.
    """
    existing_BZ_counts = get_existing_BZ_counts(output_folder_path)

    # Determine which BZs are done vs. still needed
    mask = np.array([
        existing_BZ_counts.get(tuple(bz), 0) < num_paths
        for bz in BZ_list
    ])
    
    new_BZ_list = BZ_list[mask]
    removed_BZ = BZ_list[~mask]

    print(f"{len(removed_BZ)} BZ arrays fully processed and removed:")
    print(removed_BZ.tolist())

    # Create valid tag strings like 'BZ[ 0  0 -6]'
    removed_tags = {f'BZ{bz}' for bz in removed_BZ}

    # Delete slice files not associated with a fully processed BZ
    if os.path.exists(all_slices_dir):
        for fname in os.listdir(all_slices_dir):
            full_path = os.path.join(all_slices_dir, fname)
            if fname.endswith('.nxs'):
                keep = any(tag in fname for tag in removed_tags)
                if not keep:
                    os.remove(full_path)

    return new_BZ_list, removed_BZ



def count_total_bzfold_slices(pathSQE_params, BZ_list, mtd_spacegroup):
    """
    Calculates the total number of slices expected before processing starts.

    Parameters:
    - BZ_list: list of BZ offsets
    - user_defined_Qpoints: dictionary containing 'path' and 'point_coords'
    - mtd_spacegroup: Mantid space group object used for path generation
    - pathSQE_params: dictionary containing parameters including 'primitive to mantid transformation matrix'

    Returns:
    - total_slices: int, total number of slices to be processed
    """

    num_BZ = len(BZ_list)  # Number of BZs
    unique_slices_per_BZ = 0  # Count slices for one BZ

    # Only iterate through path segments once (since they are the same for all BZs)
    for path_seg in pathSQE_params['user defined Qpoints']['path']:
        pt1_init = np.array(pathSQE_params['user defined Qpoints']['point_coords'][path_seg[0]])
        pt2_init = np.array(pathSQE_params['user defined Qpoints']['point_coords'][path_seg[1]])

        # Generate unique paths for this segment
        pt1_array, _ = core.generate_unique_paths(mtd_spacegroup, pt1_init, pt2_init, pathSQE_params)
        print("\n{} {} symmetrically equivalent path segments".format(pt1_array.shape[0], path_seg))

        # Add the number of unique paths generated (which corresponds to slices)
        unique_slices_per_BZ += pt1_array.shape[0]

    # Multiply by the number of BZs to get the total slices
    num_total_slices = unique_slices_per_BZ * num_BZ

    return num_total_slices



def count_total_symPt_slices(pathSQE_params, BZ_list, mtd_spacegroup, num_datasets):
    """
    Calculates the total number of slices expected before processing starts.

    Parameters:
    - pathSQE_params: dictionary containing parameters, including 'user defined Qpoints'
    - BZ_list: list of BZ offsets
    - mtd_spacegroup: Mantid space group object used for symmetry generation
    - num_datasets: int, number of datasets being processed

    Returns:
    - total_slices: int, estimated total number of slices to be processed
    """

    num_BZ = len(BZ_list)  # Number of BZs
    unique_slices_per_BZ = 0  # Count unique slices in a single BZ

    # Iterate over all defined 1D symmetry points
    for symPoint in pathSQE_params['user defined Qpoints']['1d_points']:
        symPt_init = np.array(pathSQE_params['user defined Qpoints']['point_coords'][symPoint])

        # Find symmetry-equivalent points in a single BZ
        symPts_array = core.generate_unique_symPts(mtd_spacegroup, symPoint, symPt_init)

        # Find unique symmetry points across all BZs
        allSymPts_array = core.generate_unique_symPts_inAllBZ(BZ_list, symPts_array)

        # Count total symmetry-equivalent slices
        unique_slices_per_BZ += allSymPts_array.shape[0]

    # Multiply by the number of datasets since each dataset requires slicing
    num_total_slices = unique_slices_per_BZ * num_datasets

    return num_total_slices



def should_print_update(num_processed_slices, num_total_slices, last_update_step, update_percent=5):
    """
    Determines whether to print a progress update based on percentage completion.

    Ensures that the first slice always triggers a print.

    Parameters:
    - num_processed_slices: int, number of slices processed so far.
    - num_total_slices: int, total number of slices.
    - last_update_step: int, last progress checkpoint (percentage completed).
    - update_percent: int, interval at which to print updates (default 5%).

    Returns:
    - should_print: bool, whether to print an update.
    - new_update_step: int, updated progress checkpoint.
    """
    progress = (num_processed_slices / num_total_slices) * 100
    
    # Always print after first slice
    if num_processed_slices == 1:
        return True, last_update_step

    # Print at each update_percent step
    if progress >= last_update_step + update_percent:
        return True, (last_update_step + update_percent)
    
    return False, last_update_step


# Convert to readable format (hh:mm:ss)
def format_time(seconds):
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{int(h):02d}:{int(m):02d}:{int(s):02d}"



def print_progress(num_processed_slices, num_total_slices, fold_timer, processing_type):
    """
    Prints the estimated time to completion based on the rolling average processing speed.

    Parameters:
    - num_processed_slices: int, number of slices processed so far
    - num_total_slices: int, total number of slices to be processed
    - fold_timer: float, timestamp when processing started (time.time())
    """

    elapsed_time = time.time() - fold_timer  # Time since start
    avg_time_per_slice = elapsed_time / max(num_processed_slices, 1)  # Avoid division by zero
    remaining_slices = num_total_slices - num_processed_slices
    estimated_remaining_time = remaining_slices * avg_time_per_slice

    print(f"\n\nProcessed: {num_processed_slices}/{num_total_slices} "
          f"({(num_processed_slices / num_total_slices) * 100:.2f}%) | "
          f"Elapsed: {format_time(elapsed_time)} | "
          f"{processing_type} estimated remaining time: {format_time(estimated_remaining_time)}\n\n")



def make_slice_desc_1DSymPoints(pathSQE_params, symPoint, point):
    # slice description with desired dims and bins and unique name
    slice_desc={'QDimension0':'1,0,0',
                'QDimension1':'0,1,0',
                'QDimension2':'0,0,1',
                'Dimension0Name':'QDimension0',
                'Dimension0Binning':conv_to_desc_string([point[0]-pathSQE_params['qdim0 step size'], point[0]+pathSQE_params['qdim0 step size']]),
                'Dimension1Name':'QDimension1',
                'Dimension1Binning':conv_to_desc_string([point[1]-pathSQE_params['qdim1 integration range'], point[1]+pathSQE_params['qdim1 integration range']]),
                'Dimension2Name':'QDimension2',
                'Dimension2Binning':conv_to_desc_string([point[2]-pathSQE_params['qdim2 integration range'], point[2]+pathSQE_params['qdim2 integration range']]),
                'Dimension3Name':'DeltaE',
                'Dimension3Binning':pathSQE_params['E bins'],
                'Name':'pathSQE_{}point_{}'.format(symPoint, point)}

    #slice_desc['Name'] = slice_desc['Name'].replace(',', '_').replace('.', '').replace('-','m')
    
    return slice_desc



def adaptive_colorbars(slice_ws_names, pathSQE_params, percentile_min=10, percentile_max=95, elastic_frac_threshold=0.05):
    """
    Computes adaptive colorbar vmin and vmax based on percentiles of signal arrays from sliced workspaces.
    
    Parameters:
        slice_ws_names (list of str): Names of the Mantid slice workspaces.
        pathSQE_params (dict): Dictionary containing keys:
            - 'E bins': str of comma-separated energy binning (e.g., "-10,0.25,50")
            - 'T and Ei conditions': tuple or list, where index 1 is Ei (incident energy)
        percentile_min (float): Lower percentile for colorbar scaling (e.g., 10).
        percentile_max (float): Upper percentile for colorbar scaling (e.g., 95).
        elastic_frac_threshold (float): Fraction of Ei used to define cutoff above the elastic line.
    
    Returns:
        tuple: (vmin, vmax) values for consistent colorbar scaling.
    """

    # Extract energy bin step and incident energy
    step = float(pathSQE_params['E bins'].split(',')[1])
    Ei = pathSQE_params['T and Ei conditions'][0][1]

    # Use a fraction of Ei to define index above elastic line
    elastic_meV = elastic_frac_threshold * Ei
    ind_above_elastic = int(np.ceil(elastic_meV / step))

    pathMax = -1
    pathMin = 1e2

    for seg in slice_ws_names:
        try:
            full_signal = np.squeeze(mtd[seg].getSignalArray())
            slice_signal = full_signal[:, ind_above_elastic:]

            mask = ~np.isnan(slice_signal) & (slice_signal > 0)
            if np.any(mask):
                segMin = np.percentile(slice_signal[mask], percentile_min)
                segMax = np.percentile(slice_signal[mask], percentile_max)

                if segMax > pathMax:
                    pathMax = segMax
                if segMin < pathMin:
                    pathMin = segMin
            else:
                print(f"Warning: No valid signal values in segment {seg}, skipping.")
        except Exception as e:
            print(f"Error processing segment {seg}: {e}")

    return pathMin, pathMax

