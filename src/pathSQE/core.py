import numpy as np
import os
import matplotlib.pyplot as plt
from mantid.simpleapi import *
from mantid.geometry import SpaceGroupFactory
from . import slice_utils_07142023
from . import helper


def find_BZ_with_data(mde_data, pathSQE_params, hkl_range=(-10, 10, -10, 10, -10, 10)): 
    """
    Identify and rank Brillouin Zones (BZs) by fractional data coverage.

    Parameters:
        mde_data: Mantid workspace containing the data.
        pathSQE_params (dict): Dictionary containing parameters, including the transformation matrix.
        hkl_range (tuple): Defines BZ grid in primitive basis.

    Returns:
        list: Ranked list of BZ centers sorted by fractional coverage (descending).
    """
    E_bins = pathSQE_params['E bins'].split(',')
    min_fraction = pathSQE_params['BZ to process']

    # Transform to current basis
    prim2current = pathSQE_params['primitive to mantid transformation matrix'].copy().T

    # Generate all integer BZ centers in primitive basis
    Q_primitive = np.array([[h, k, l] for h in range(hkl_range[0], hkl_range[1] + 1)
                                      for k in range(hkl_range[2], hkl_range[3] + 1)
                                      for l in range(hkl_range[4], hkl_range[5] + 1)])

    #print(prim2current[0], prim2current[1], prim2current[2])

    slice_desc = {
        'QDimension0': helper.conv_to_desc_string(prim2current[0]),  # make it current2prim columns instead of simple HKL aligned basis
        'QDimension1': helper.conv_to_desc_string(prim2current[1]),
        'QDimension2': helper.conv_to_desc_string(prim2current[2]),
        'Dimension0Name': 'QDimension0',
        'Dimension0Binning': f'{hkl_range[0] - 0.5},0.5,{hkl_range[1] + 0.5}',  # Step size = 0.5 for octants
        'Dimension1Name': 'QDimension1',
        'Dimension1Binning': f'{hkl_range[2] - 0.5},0.5,{hkl_range[3] + 0.5}',
        'Dimension2Name': 'QDimension2',
        'Dimension2Binning': f'{hkl_range[4] - 0.5},0.5,{hkl_range[5] + 0.5}',
        'Dimension3Name': 'DeltaE',
        'Dimension3Binning': f'{E_bins[0]},{0.2*(float(E_bins[2])-float(E_bins[0]))},{E_bins[2]}',
        'Name': 'pathSQE_BZwithData'
    }

    slice_utils_07142023.make_slice(mde_data, slice_desc)
    data_slice = mtd[slice_desc['Name']].getSignalArray()
    #print('cov slice shape', data_slice.shape)
    #np.save(os.path.join(pathSQE_params['output directory'],'BZ_coverage_slice.npy'), data_slice)

    # Step 3: Compute BZ fractional coverage
    BZ_list = []
    bin_edges_h = np.arange(hkl_range[0] - 0.5, hkl_range[1] + 1.0, 0.5)
    bin_edges_k = np.arange(hkl_range[2] - 0.5, hkl_range[3] + 1.0, 0.5)
    bin_edges_l = np.arange(hkl_range[4] - 0.5, hkl_range[5] + 1.0, 0.5)

    BZ_full_list = []  # Store all BZ centers with their fractional coverages

    for q_center in Q_primitive:
        Q_current= q_center @ prim2current # used to be .T here
        h_idx = np.digitize(q_center[0] + np.array([-0.25, 0.25]), bin_edges_h) - 1
        k_idx = np.digitize(q_center[1] + np.array([-0.25, 0.25]), bin_edges_k) - 1
        l_idx = np.digitize(q_center[2] + np.array([-0.25, 0.25]), bin_edges_l) - 1

        num_bins = 0
        num_valid_bins = 0

        for hi in h_idx:
            for ki in k_idx:
                for li in l_idx:
                    octant_data = data_slice[hi, ki, li, :]
                    num_bins += np.size(octant_data)
                    num_valid_bins += np.sum(np.logical_and(~np.isnan(octant_data), octant_data != 0))

        # Compute the fractional coverage
        fraction_coverage = num_valid_bins / num_bins if num_bins > 0 else 0

        # Only include BZs with coverage > 0
        if fraction_coverage > 0:
            BZ_full_list.append((Q_current, fraction_coverage))

        # Keep only those meeting the threshold for later analysis
        if fraction_coverage >= min_fraction:
            BZ_list.append((Q_current, fraction_coverage))

    # Sort by coverage descending, then by norm of BZ (column 0) descending
    BZ_full_array = np.array(BZ_full_list, dtype=object)
    coverage = np.array([entry[1] for entry in BZ_full_array])
    bz_norms = np.linalg.norm([entry[0] for entry in BZ_full_array], axis=1)
    sort_idx = np.lexsort(( -bz_norms, -coverage ))  # Primary: coverage, Secondary: BZ norm
    BZ_full_array = BZ_full_array[sort_idx]

    # Save to file
    np.savetxt(os.path.join(pathSQE_params['output directory'],"BZ_coverage_sorted.txt"), 
               np.column_stack([BZ_full_array[:, 0], BZ_full_array[:, 1]]), 
               fmt='%s', header="BZ_center   Fraction_Coverage")

    # --- Scatter Plot: Coverage vs. Rank ---
    plt.figure(figsize=(6, 4))
    plt.scatter(range(len(BZ_full_array)), BZ_full_array[:, 1], color='tab:blue', s=10, label="BZ Coverage")
    plt.axhline(y=min_fraction, color='tab:grey', linestyle='--', label=f"Threshold ({min_fraction:.2f})")
    plt.xlabel("BZ Rank")
    plt.ylabel("Fractional Coverage")
    plt.title("BZ Coverage vs. Rank")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(pathSQE_params['output directory'],'BZ_coverage.png'))
    plt.close()

    # --- 2D Scatter Plot: BZs with Data in u,v plane ---
    # Project BZ Centers to (u, v) Plane
    BZ_uv_all = helper.project_to_uv_plane(BZ_full_list, pathSQE_params['u_vec'], pathSQE_params['v_vec'])  # All BZs with data
    BZ_uv_threshold = helper.project_to_uv_plane(BZ_list, pathSQE_params['u_vec'], pathSQE_params['v_vec'])  # BZs above coverage threshold
    np.save(os.path.join(pathSQE_params['output directory'], 'BZ_uv_all.npy'), BZ_uv_all)
    np.save(os.path.join(pathSQE_params['output directory'], 'BZ_uv_threshold.npy'), BZ_uv_threshold)

    plt.figure(figsize=(6, 6))
    plt.scatter(BZ_uv_all[:, 0], BZ_uv_all[:, 1], color='black', s=10, label="BZ with Data")
    plt.scatter(BZ_uv_threshold[:, 0], BZ_uv_threshold[:, 1], color='red', s=20, label=f"BZ > {min_fraction:.2f} Coverage")

    plt.xlabel("{}".format(pathSQE_params['u_vec']))
    plt.ylabel("{}".format(pathSQE_params['v_vec']))
    plt.title("BZ Centers with Data (Projected to u-v Plane)")
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.savefig(os.path.join(pathSQE_params['output directory'], 'BZ_uv_scatter.png'))
    plt.close()

    # Step 4: Sort and return the ranked BZ list
    BZ_list.sort(key=lambda x: (x[1], np.linalg.norm(x[0])), reverse=True)  # Sort descending by coverage, secondary by BZ norm

    # Extract separate lists for BZ centers and their corresponding fractional coverages
    BZ_centers = [bz[0] for bz in BZ_list]
    fractional_coverages = [bz[1] for bz in BZ_list]

    return BZ_centers, fractional_coverages, BZ_full_array 





def find_qdim_1and2(qdim0, u, v):
    """
    Find qdim1 and qdim2 given qdim0, ensuring qdim1 lies in the plane defined by u and v.

    Parameters:
        qdim0 (array): The difference vector between input Q points, defining the path direction.
        u (array): A vector defining the first basis direction of the horizontal scattering plane.
        v (array): A vector defining the second basis direction of the horizontal scattering plane.

    Returns:
        qdim1 (array): A vector perpendicular to qdim0 and lying in the plane defined by (u, v).
        qdim2 (array): A vector perpendicular to both qdim0 and qdim1.
    """
    # Normalize input vectors
    qdim0 = qdim0.astype(float) / np.linalg.norm(qdim0)
    u = u.astype(float) / np.linalg.norm(u)
    v = v.astype(float) / np.linalg.norm(v)

    # Compute the normal to the scattering plane
    z_scattering = np.cross(u, v)
    z_scattering /= np.linalg.norm(z_scattering)

    # Attempt to construct qdim1 as the cross product of qdim0 and z_scattering
    qdim1 = np.cross(z_scattering, qdim0)

    if np.linalg.norm(qdim1) < 1e-10:  # qdim0 is collinear with z_scattering
        # Find a new qdim1 as a linear combination of u and v, ensuring perpendicularity to qdim0
        qdim1 = u - np.dot(u, qdim0) * qdim0  # Start with u and remove its component along qdim0

        # If still degenerate, use v instead
        if np.linalg.norm(qdim1) < 1e-10:
            qdim1 = v - np.dot(v, qdim0) * qdim0

        # Normalize qdim1
        if np.linalg.norm(qdim1) < 1e-10:
            raise ValueError("Failed to construct qdim1. Check input vectors for collinearity issues.")
    
    qdim1 /= np.linalg.norm(qdim1)

    # Compute qdim2 as the cross product (ensuring it is perpendicular to both qdim0 and qdim1)
    qdim2 = np.cross(qdim0, qdim1)

    # Normalize qdim2 safely
    if np.linalg.norm(qdim2) > 1e-10:
        qdim2 /= np.linalg.norm(qdim2)
    else:
        raise ValueError("Failed to construct qdim2. Check for collinearity issues.")

    return qdim1, qdim2



def find_qbins(pathSQE_params, qdim0, qdim1, qdim2, pt1, pt2, BZ_offset):
    # we start with bins of 0 so hopefully can catch it if code isn't working
    qdim1_range = 0
    qdim2_range = 0 

    # find the new path endpoint coords given change of basis from cart to qdim0,1,2
    CoB_mat = np.linalg.inv(np.array([qdim0,qdim1,qdim2]).T) # back to og, edited from Bi version with np.array([qdim0,qdim1,qdim2])
    path_start_transf = CoB_mat @ (pt1 + BZ_offset)
    path_end_transf = CoB_mat @ (pt2 + BZ_offset)
    print('Slicing ', pt1+BZ_offset, ' to ', pt2+BZ_offset)
    #print('find_qbins start and end: ', path_start_transf, path_end_transf)
    
    # Compute absolute step size based on fractional step // commenting bc not sure this is what we want..
    #qdim0_step_size = np.abs((path_end_transf[0] - path_start_transf[0]) * pathSQE_params['qdim0 step size'])

    # Construct the array with the computed step size
    qdim0_range = np.array([path_start_transf[0], pathSQE_params['qdim0 step size'], path_end_transf[0]])
    
    # might need to check bin/path endpoints...
    # old check_bin0_and_path_endpoints function
    
    # fixed integration range for qdim1 and qdim2
    qdim1_int_range = pathSQE_params['qdim1 integration range']
    qdim1_range = np.array([path_start_transf[1]-qdim1_int_range, path_end_transf[1]+qdim1_int_range])
    
    qdim2_int_range = pathSQE_params['qdim2 integration range']
    qdim2_range = np.array([path_start_transf[2]-qdim2_int_range, path_end_transf[2]+qdim2_int_range])
    
    # if qdim1 or qdim2 bin range gets too big raise error because code above may be wrong...
    #if (np.absolute(qdim1_range[1]-qdim1_range[0]) > 0.2) or (np.absolute(qdim2_range[1]-qdim2_range[0]) > 0.2):
    #    raise ValueError("Might be something weird with bins for qdim1 or dqim2")
    
    return qdim0_range, qdim1_range, qdim2_range



def choose_dims_and_bins(pathSQE_params, point1, point2, u, v, BZ_offset=np.array([0,0,0])):
    diff = point2 - point1

    # calculating default QDimension0 axis direction and scaling
    qdim0 = diff / np.max(np.absolute(diff))

    # determine the directions of qdim1 and qdim2
    qdim1, qdim2 = find_qdim_1and2(qdim0, u, v)
    #print('qdims ', qdim0, qdim1, qdim2)
    qdim0_range, qdim1_range, qdim2_range = find_qbins(pathSQE_params, qdim0, qdim1, qdim2, point1, point2, BZ_offset)

    # return list contains various types (mainly np arrays) that are converted to properly formatted strings
    # as a part of the slice description generation process
    q_dims_and_bins = [qdim0, qdim1, qdim2, qdim0_range, qdim1_range, qdim2_range]
    return q_dims_and_bins




def generate_unique_paths(mtd_spacegroup, pt1, pt2, pathSQE_params):
    prim2mantid = pathSQE_params['primitive to mantid transformation matrix']
    hex_rhomb = pathSQE_params['rhomb in hex notation']

    # Get the symmetry operations for the spacegroup
    spacegroup = SpaceGroupFactory.createSpaceGroup(mtd_spacegroup)
    pg = spacegroup.getPointGroup()
    #print(pg)
    
    rotations = [np.array([so.transformHKL((1, 0, 0)), so.transformHKL((0, 1, 0)), so.transformHKL((0, 0, 1))])
                 for so in pg.getSymmetryOperations()]

    if hex_rhomb:
        new_pt1s = rotations @ pt1 @ prim2mantid.T
        new_pt2s = rotations @ pt2 @ prim2mantid.T
    else:
        new_pt1s = rotations @ pt1
        new_pt2s = rotations @ pt2       

    # Combine new_pt1s and new_pt2s into a single array for comparison
    combined_array = np.column_stack((new_pt1s, new_pt2s))
    combined_array = np.round(combined_array, 2)

    # Find unique rows in the combined array
    unique_combined_array = np.unique(combined_array, axis=0)

    # Split the unique_combined_array back into separate arrays for pt1 and pt2
    unique_pt1_array = unique_combined_array[:, :3]
    unique_pt2_array = unique_combined_array[:, 3:]

    return unique_pt1_array, unique_pt2_array



def combine_data_within_bz(data, norm, mtd_data, mtd_norm):
    data.setSignalArray(data.getSignalArray() + mtd_data.getSignalArray())
    norm.setSignalArray(norm.getSignalArray() + mtd_norm.getSignalArray())
    data.setErrorSquaredArray(data.getErrorSquaredArray() + mtd_data.getErrorSquaredArray())
    return data, norm



def fold_BZs(BZ_fold_list, base_BZ, pathSQE_params, user_defined_Qpoints, rerun=False):
    """
    Fold multiple Brillouin zones' MD data files into a single dataset.

    Parameters:
        BZ_fold_list (list): List of BZ indices to combine (determined externally).
        BZ_offset (int): A reference BZ index used as the base for folding.
        pathSQE_params (dict): Contains 'output directory' as a key for file paths.
        user_defined_Qpoints (dict): Contains 'path', a list of label pairs used to generate segment names.
    """
    print('\nCombining all {} BZ(s)...'.format(len(BZ_fold_list)))

    data_and_norm_dir = os.path.join(pathSQE_params['output directory'], 'slices/intraBZ_folded/DataAndNorms/')
    all_files = os.listdir(data_and_norm_dir)

    # Use base zone to accumulate data into and others to combine
    base_BZ_files = [f for f in all_files if str(base_BZ) in f]

    # FIXED: ensure safe string comparison without NumPy ambiguity
    BZ_list_files = [
        f for f in all_files
        if any(str(bz) in f for bz in map(str, BZ_fold_list) if str(bz) != str(base_BZ))
    ]

    # Generate segment names like 'G2X', 'X2U', etc.
    seg_names = ['{}2{}'.format(p[0], p[1]) for p in user_defined_Qpoints['path']]

    for path_seg in seg_names:
        print(f"  Folding segment: {path_seg}")

        # Files from the reference BZ (base_BZ)
        pathSeg_files_baseBZ = [f for f in base_BZ_files if path_seg in f]
        first_data_file = next((f for f in pathSeg_files_baseBZ if 'data' in f), None)
        first_norm_file = next((f for f in pathSeg_files_baseBZ if 'norm' in f), None)

        if first_data_file is None or first_norm_file is None:
            print(f"    Skipping {path_seg}: Missing data or norm files in reference BZ {base_BZ}")
            continue

        data = LoadMD(Filename=os.path.join(data_and_norm_dir, first_data_file))
        norm = LoadMD(Filename=os.path.join(data_and_norm_dir, first_norm_file))

        # Files from other BZs to fold in
        pathSeg_files_otherBZs = [f for f in BZ_list_files if path_seg in f]
        data_files = [f for f in pathSeg_files_otherBZs if 'data' in f]
        norm_files = [f for f in pathSeg_files_otherBZs if 'norm' in f]

        if len(data_files) != len(norm_files):
            print(f"    Warning: Mismatched data/norm file count for segment {path_seg}")

        for dfile, nfile in zip(data_files, norm_files):
            data_temp = LoadMD(Filename=os.path.join(data_and_norm_dir, dfile), LoadHistory=False)
            norm_temp = LoadMD(Filename=os.path.join(data_and_norm_dir, nfile), LoadHistory=False)

            data.setSignalArray(data.getSignalArray() + data_temp.getSignalArray())
            norm.setSignalArray(norm.getSignalArray() + norm_temp.getSignalArray())
            data.setErrorSquaredArray(data.getErrorSquaredArray() + data_temp.getErrorSquaredArray())

        # Normalize and save the folded result
        folded_ws = data / norm
        if rerun:
            output_dir = os.path.join(pathSQE_params['output directory'], 'Reruns')
        else:
            output_dir = os.path.join(pathSQE_params['output directory'], 'slices/fully_folded')
        SaveMD(folded_ws, Filename=output_dir+f'/{path_seg}_folded.nxs')
        SaveMD(data, Filename=output_dir+f'/{path_seg}_folded_data.nxs')
        SaveMD(norm, Filename=output_dir+f'/{path_seg}_folded_norm.nxs')
        LoadMD(Filename=output_dir+f'/{path_seg}_folded.nxs', OutputWorkspace=path_seg, LoadHistory=False)

        print(f"    Saved folded {path_seg} segment to: {output_dir}")

    return seg_names



def generate_unique_symPts(mtd_spacegroup, symPoint, symPt_init):
    spacegroup = SpaceGroupFactory.createSpaceGroup(mtd_spacegroup)
    pg = spacegroup.getPointGroup()
    rotations = [np.array([so.transformHKL((1, 0, 0)), so.transformHKL((0, 1, 0)), so.transformHKL((0, 0, 1))])
                 for so in pg.getSymmetryOperations()]

    new_symPts = np.dot(rotations, symPt_init)

    # Find unique rows
    unique_symPts_array = np.unique(new_symPts, axis=0)

    return unique_symPts_array



def generate_unique_symPts_inAllBZ(BZ_list, symPts_array):
    unique_points_set = set()
    
    # Iterate through each BZ center, which is a NumPy array
    for bz_center in BZ_list:
        # Translate symPts_array for the current BZ center
        translated_points = symPts_array + bz_center
        
        # Add each translated point to the set (as a tuple to ensure hashability)
        for point in translated_points:
            unique_points_set.add(tuple(point))
    
    # Convert the set back to an Nx3 numpy array
    allSymPts_array = np.array(list(unique_points_set))
    
    return allSymPts_array
