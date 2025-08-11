import numpy as np
import os
import re

from . import Resolution as Res


def load_previous_simulation_progress(output_folder, user_defined_Qpoints, removed_BZ):
    """
    Loads and combines previously saved simulation data and norm files, but only for BZs that are fully processed.

    Parameters:
        output_folder (str): Path to the folder containing saved simulation results.
        user_defined_Qpoints (dict): Dictionary containing the path information.
        removed_BZ (np.ndarray): Nx3 array of BZ centers that were fully processed.

    Returns:
        tuple: (datas, norms) where:
            - datas: List of summed data arrays, one per path segment.
            - norms: List of summed norm arrays, one per path segment.
            - Returns ([], []) if no valid data is found.
    """
    datas = []
    norms = []

    # Convert removed_BZ to string format for easier matching in filenames
    removed_BZ_strs = {str(bz.tolist()) for bz in removed_BZ}
    print('BZ strings for match:', removed_BZ_strs)

    # Initialize empty lists for each segment (but return [] if nothing is found)
    num_segments = len(user_defined_Qpoints['path'])
    found_any_files = False  # Track if any valid files were found
    datas = [None] * num_segments
    norms = [None] * num_segments

    # Search for all BZ-related files in the output folder
    for filename in os.listdir(output_folder):
        match = re.search(r'BZ(\[.*?\])_foldedSim_data.npz', filename)  # Extract BZ array from filename
        if match:
            print('Match found:', filename)
            BZ_str = match.group(1)  # Extract the array portion inside brackets
            BZ_array = np.array([float(x) for x in BZ_str.strip("[]").split()])  # Convert to array

            # Check if this BZ is fully processed
            if str(BZ_array.tolist()) in removed_BZ_strs:
                data_path = os.path.join(output_folder, filename)
                norm_path = data_path.replace("_foldedSim_data.npz", "_foldedSim_norm.npz")

                print('Checking files:', data_path, norm_path)

                if os.path.exists(data_path) and os.path.exists(norm_path):
                    print('Files exist - loading data')
                    found_any_files = True  # At least one file was found

                    with np.load(data_path, allow_pickle=True) as data_npz:
                        data_arrays = [data_npz[key] for key in data_npz]  # Extract all arrays

                    with np.load(norm_path, allow_pickle=True) as norm_npz:
                        norm_arrays = [norm_npz[key] for key in norm_npz]  # Extract all arrays

                    # Add data to corresponding segment index
                    for seg_idx in range(num_segments):
                        if datas[seg_idx] is None:
                            datas[seg_idx] = data_arrays[seg_idx]
                            norms[seg_idx] = norm_arrays[seg_idx]
                        else:
                            datas[seg_idx] += data_arrays[seg_idx]
                            norms[seg_idx] += norm_arrays[seg_idx]

    # If no valid files were found, return empty lists
    if not found_any_files:
        print('No previous simulation data found.')
        return [], []

    print('Final data and norm lists populated')
    return datas, norms



def determine_forces_file():
    """Determine whether FORCE_CONSTANTS or FORCE_SETS is present in the directory."""
    if os.path.exists("FORCE_CONSTANTS"):
        forceconstants_file = "FORCE_CONSTANTS"
        forcesets_file = None
    elif os.path.exists("FORCE_SETS"):
        forceconstants_file = None
        forcesets_file = "FORCE_SETS"
    else:
        raise FileNotFoundError("Error: Neither FORCE_CONSTANTS nor FORCE_SETS found in the current directory.")

    return forceconstants_file, forcesets_file



def get_coherent_scattering_lengths(sample_formula, filename="pathSQE_utils/ScattLengths.txt"):
    """
    Given a sample's chemical formula, find the corresponding coherent scattering lengths.

    Parameters:
        sample_formula (str): Sample composition as a space-separated string (e.g., "Fe Si").
        filename (str): Name of the scattering length data file.

    Returns:
        dict: Coherent scattering lengths for each element in the sample.
    """
    # Check if the scattering lengths file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Error: {filename} not found.")

    # Read the scattering lengths file and store values in a dictionary
    coh_scatter_length = {}

    with open(filename, "r") as file:
        lines = file.readlines()

        for line in lines:
            # Skip headers or malformed lines
            if line.startswith("atom_name") or line.strip() == "":
                continue
            
            parts = line.split()
            if len(parts) < 3:
                continue  # Ensure enough columns exist
            
            element = parts[0]  # First column is the element name
            coh_b = parts[1]  # Second column is the coherent scattering length

            # Some entries may have complex values; extract real part if needed
            try:
                coh_b = float(coh_b.split("-")[0])  # Extract real part (ignoring imaginary)
            except ValueError:
                continue  # Skip if conversion fails (e.g., for weird formats)

            coh_scatter_length[element] = coh_b

    # Extract scattering lengths for the sample elements
    sample_elements = sample_formula.split()
    sample_scattering = {elem: coh_scatter_length[elem] for elem in sample_elements if elem in coh_scatter_length}

    return sample_scattering



def construct_sim_Qpts(pt1, pt2, q_diff, step_size, prim2mantid, BZ_offset=np.array([0,0,0])):
    """
    Constructs simulated Q points along the direction of q_diff.

    Parameters:
    - pt1: ndarray, starting point of the Q trajectory
    - pt2: ndarray, ending point of the Q trajectory
    - BZ_offset: ndarray, offset to apply to both start and end points
    - q_diff: ndarray, direction vector along which Q points are generated
    - step_size: float, step size along the q_diff direction

    Returns:
    - Qpoints: list of ndarrays, generated Q points
    """
    from scipy.linalg import norm

    q_start = pt1 + BZ_offset
    q_end = pt2 + BZ_offset  # Explicitly set the correct endpoint

    # Compute the projection of q_end onto q_start + q_diff * t
    total_distance = np.max(np.abs(q_end - q_start))
    step_size = step_size / 4 # for finer sim                            # EDIT USED TO BE /2
    num_steps = int(np.round(total_distance / step_size))

    # Generate Q points
    Qpoints = [q_start + i * step_size * q_diff for i in np.arange(0.5, num_steps)]
    #print('len Q pts ', len(Qpoints), Qpoints[0], Qpoints[-1])

    # Avoid DW singularity at Q = [0,0,0]
    Qpoints = np.array(Qpoints)
    #print('Q pt array shape ', Qpoints.shape)
    zero_index = np.where(norm(Qpoints, axis=1) == 0)[0]
    if zero_index.size > 0:
        Qpoints[zero_index[0]] = np.array([1e-6, 1e-6, 1e-6])

    Qpoints = Qpoints @ np.linalg.inv(prim2mantid)

    return Qpoints



def run(phonon, Qpoints, temperature, atomic_form_factor_func=None, scattering_lengths=None):
    from phonopy import load

    # Transformation to the Q-points in reciprocal primitive basis vectors
    Q_prim = np.dot(Qpoints, phonon.primitive_matrix)
    # Q_prim must be passed to the phonopy dynamical structure factor code.
    phonon.run_dynamic_structure_factor(
        Q_prim,
        temperature,
        atomic_form_factor_func=atomic_form_factor_func,
        scattering_lengths=scattering_lengths,
        freq_min=8e-2)
    dsf = phonon.dynamic_structure_factor
    q_cartesian = np.dot(dsf.qpoints,
                         np.linalg.inv(phonon.primitive.get_cell()).T)
    distances = np.sqrt((q_cartesian ** 2).sum(axis=1))

    SandE = np.array([dsf.frequencies,dsf.dynamic_structure_factors])
    return SandE



def sim_SQE(pathSQE_params, Qpoint, Temperature):
    from phonopy import load

    #### START OF USER INPUTS #####
    #print('Q pt ', Qpoint)   # needs to be a list of arrays so like [array([-1.5, -1.5,  0.5])] when printed
    ## Q inputs ##
    primitive_cell = [[1,0,0],[0,1,0],[0,0,1]]
    supercell = pathSQE_params['supercell dimensions']
    forceconstants_file, forcesets_file = determine_forces_file()

    # Neutron coherent scattering length can be found at https://www.ncnr.nist.gov/resources/n-lengths/
    coh_scatter_length = get_coherent_scattering_lengths(pathSQE_params['sample'])
    THz2meV = 4.1357

    #### END OF USER INPUTS ####
    phonon = load(supercell_matrix=supercell,
                  primitive_matrix=primitive_cell,
                  unitcell_filename="POSCAR",
                  force_sets_filename=forcesets_file,
                  force_constants_filename=forceconstants_file
                  )

    #print(Qpoints)
    # Mesh sampling phonon calculation is needed for Debye-Waller factor.
    # This must be done with is_mesh_symmetry=False and with_eigenvectors=True.
    mesh = [11, 11, 11]
    phonon.run_mesh(mesh,
                    is_mesh_symmetry=False, # symmetry must be off
                    with_eigenvectors=True) # eigenvectors must be true
    temperature = Temperature

    # For INS, scattering length has to be given.
    # The following values is obtained at (Coh b)
    # https://www.nist.gov/ncnr/neutron-scattering-lengths-list
    output = run(phonon,
                 Qpoint,
                 temperature,
                 scattering_lengths=coh_scatter_length)
    
    ## output has shape as (2,len(Qpoints),branches), the [0,:,:] is for frequency and the [1,:,:] is for SQE; The frequency is in THz unit
    for i in range(len(Qpoint)):
        output[0,i,:] *= THz2meV

    return output



def SQE_to_1d_spectrum(pathSQE_params, output):
    frequencies = output[0, 0, :]
    intensities = output[1, 0, :]

    # define spectral binning info
    E_min = float(pathSQE_params['E bins'].split(',')[0])
    E_max = float(pathSQE_params['E bins'].split(',')[2])
    
    e_resolution = 1.177 * pathSQE_params['energy blurring sigma'] * float(pathSQE_params['E bins'].split(',')[1]) # gauss sigma to lorentz hwhm

    # Tolerance for considering frequencies equal (0.1%)
    tolerance = 0.001

    # Initialize arrays to store unique frequencies and their corresponding summed intensities
    unique_freq = []
    summed_intensities = []

    # Iterate over the frequency array
    for freq in np.unique(frequencies):
        # Exclude frequencies less than 0.1 to avoid singularity at gamma
        if freq < 0.1:
            continue
        
        # Check if the current frequency is sufficiently separated from existing unique frequencies
        if all(abs(freq - existing_freq) > tolerance * existing_freq for existing_freq in unique_freq):
            # Find indices where the current frequency occurs in the original array
            indices = np.where(np.abs(frequencies - freq) / freq <= tolerance)[0]
            # Sum corresponding intensities for the current frequency
            total_intensity = np.sum(intensities[indices])
            # Append the unique frequency and its corresponding summed intensity
            unique_freq.append(freq)
            summed_intensities.append(total_intensity)


    # Convert lists to numpy arrays
    unique_freq = np.array(unique_freq)
    summed_intensities = np.array(summed_intensities)
    
    # Generate frequency range
    freq_range = np.arange(E_min, E_max, 0.02)
    
    # Initialize spectrum array
    spectrum = np.zeros(freq_range.shape)
    ind_spec = []

    # Iterate over each peak
    for freq, intensity in zip(unique_freq, summed_intensities):
        # Add the Gaussian function to the spectrum
        #gaussian_peak = intensity * np.exp(-((freq_range - freq) / (2 * e_resolution)) ** 2)
        lorentzian_peak = intensity / (1 + ((freq_range - freq) / e_resolution) ** 2)
        ind_spec.append(lorentzian_peak)
        spectrum += lorentzian_peak


    return unique_freq, summed_intensities, freq_range, spectrum



def SQE_to_2d_spectrum(pathSQE_params, output): 
    """
    Bins simulated SQE data to a finer grid, applies smoothing, 
    and then rebins to match experimental resolution.

    Parameters:
    - pathSQE_params: dict, contains experimental binning information
    - output: ndarray, simulated data (fine resolution)

    Returns:
    - BinnedSQE: ndarray, rebinned SQE to match experimental binning
    """
    from scipy.ndimage import gaussian_filter

    E_min = float(pathSQE_params['E bins'].split(',')[0])
    E_max = float(pathSQE_params['E bins'].split(',')[2])
    E_step = float(pathSQE_params['E bins'].split(',')[1])
    
    nql_fine = output.shape[1]  # Fine Q resolution from simulation
    finer_E_step = E_step / 4   # Simulated data has finer E step                           # THIS AND BELOW 2
    nql = nql_fine // 4         # Experimental Q resolution (factor of 2 binning)
    
    # Generate energy bin edges
    E_bin_edges = np.arange(E_min, E_max + E_step, E_step)
    ne_exp = len(E_bin_edges) - 1

    # Fine binning storage
    FineBinnedSQE = np.zeros((nql_fine, len(np.arange(E_min, E_max, finer_E_step))))
    
    # Step 1: Bin energy values finely
    for ih in range(nql_fine):  # qpoints
        for j in range(len(output[0, 0, :])):  # branches
            energy_val = output[0][ih][j]
            hist, _ = np.histogram([energy_val], bins=np.arange(E_min, E_max + finer_E_step, finer_E_step), weights=[output[1][ih][j]])
            FineBinnedSQE[ih, :] += hist

    # Step 2: Apply Gaussian smoothing
    FineBinnedSQE = gaussian_filter(FineBinnedSQE, sigma=(0.5, pathSQE_params['energy blurring sigma']*4))                               # THIS ONE

    # Step 3: Rebin Q dimension (fine → coarse) using NumPy reshape and sum
    CoarseBinnedSQE = FineBinnedSQE.reshape(nql, 4, -1).sum(axis=1)  # Sum pairs of adjacent fine Q bins             # THIS ALSO 2

    # Step 4: Rebin E dimension (fine → coarse) using NumPy histogram
    BinnedSQE = np.zeros((nql, ne_exp))
    for ih in range(nql):
        hist, _ = np.histogram(np.linspace(E_min, E_max, FineBinnedSQE.shape[1]), bins=E_bin_edges, weights=CoarseBinnedSQE[ih, :])
        BinnedSQE[ih, :] = hist

    return BinnedSQE


def SQE_to_2d_spectrum_advancedRes(pathSQE_params, output): 
    """
    Bins simulated SQE data to a finer grid, applies smoothing, 
    and then rebins to match experimental resolution.

    Parameters:
    - pathSQE_params: dict, contains experimental binning information
    - output: ndarray, simulated data (fine resolution)

    Returns:
    - BinnedSQE: ndarray, rebinned SQE to match experimental binning
    """

    res_params = pathSQE_params['resolution blurring']
    ElasticFWHM = res_params[0]
    QResolution = res_params[1]/(2.35*pathSQE_params['qdim0 step size'])
    ResolutionType = "instrument"
    InstrumentName = res_params[2]
    IncidentEnergy = pathSQE_params['T and Ei conditions'][0][1]

    E_min = float(pathSQE_params['E bins'].split(',')[0])
    E_max = float(pathSQE_params['E bins'].split(',')[2])
    E_step = float(pathSQE_params['E bins'].split(',')[1])
    
    nql_fine = output.shape[1]  # Fine Q resolution from simulation
    finer_E_step = E_step / 4   # Simulated data has finer E step                           # THIS AND BELOW 2
    nql = nql_fine // 4         # Experimental Q resolution (factor of 2 binning)
    
    # Generate energy bin edges
    E_bin_edges = np.arange(E_min, E_max+E_step, E_step)
    evalues = np.arange(E_min+0.5*finer_E_step, E_max, finer_E_step)
    ne_exp = len(E_bin_edges) - 1

    # Fine binning storage
    FineBinnedSQE = np.zeros((nql_fine, len(np.arange(E_min, E_max, finer_E_step))))

    #create the resolution function object
    inst = Res.Instrument(InstrumentName,ElasticFWHM,IncidentEnergy)
    resolution = Res.Resolution(type=ResolutionType,inst=inst,sigmaq=QResolution)

    # Energy binning with convolution
    for ih in range(nql_fine):  # q-points
        for j in range(len(output[0, 0, :])):  # phonon branches
            E_phonon = output[0][ih][j]  # in meV (already converted from THz)
            intensity = output[1][ih][j]
            #print(ih, j, intensity)
            
            conv = resolution.Gauss(evalues, 0, E_phonon, 0, escale=2.4, qscale=1)
            
            FineBinnedSQE[ih, :] += intensity * conv

    
    from scipy.ndimage import gaussian_filter1d

    # Use ndimage Gaussian blur along the Q direction (axis 0)
    FineBinnedSQE = gaussian_filter1d(FineBinnedSQE, sigma=QResolution, axis=0, mode='nearest')

    # Step 3: Rebin Q dimension (fine → coarse) using NumPy reshape and sum
    CoarseBinnedSQE = FineBinnedSQE.reshape(nql, 4, -1).sum(axis=1)  # Sum pairs of adjacent fine Q bins             # THIS ALSO 2

    # Step 4: Rebin E dimension (fine → coarse) using NumPy histogram
    BinnedSQE = np.zeros((nql, ne_exp))
    for ih in range(nql):
        hist, _ = np.histogram(np.linspace(E_min, E_max, FineBinnedSQE.shape[1]), bins=E_bin_edges, weights=CoarseBinnedSQE[ih, :])
        BinnedSQE[ih, :] = hist

    return BinnedSQE



def fold_symmetric_simulations(all_sims_in_BZ):
    """
    Folds all symmetrically equivalent simulated data by summing and normalizing.

    Parameters:
    - all_sims_in_BZ: List of lists of binned SQE arrays.
                      Outer list indexes path segments, inner lists contain equivalent simulations.

    Returns:
    - folded_results: List of tuples, each containing:
        (summed_data, summed_norm, final_folded_data) for each segment.
    """
    folded_results = []

    for segment_sims in all_sims_in_BZ:
        # Ensure the segment contains at least one simulation
        if not segment_sims:
            continue  
        
        # Initialize summed arrays with zeros
        summed_data = np.zeros_like(segment_sims[0])
        summed_norm = np.zeros_like(segment_sims[0])

        for sim in segment_sims:
            norm = np.where(~np.isnan(sim), 1, 0)  # Norm is 1 where data is valid, 0 otherwise
            sim_filled = np.nan_to_num(sim)  # Replace NaNs with zero for summation
            
            summed_data += sim_filled
            summed_norm += norm
        
        # Avoid division by zero by setting invalid norms to NaN
        summed_norm[summed_norm == 0] = np.nan
        
        # Compute final folded data
        final_folded_data = summed_data / summed_norm

        folded_results.append((np.nan_to_num(summed_data, nan=0.0), np.nan_to_num(summed_norm, nan=0.0), final_folded_data))

    return folded_results
