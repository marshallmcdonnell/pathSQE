########################################################################################################
# Define dictionary input file for running pathSQE
# Author: Aiden Sable. August 2025.
########################################################################################################


import numpy as np

def define_pathSQE_params(**kwargs):
    pathSQE_params={
    # Sample and Q point info
    'sample':'Sn S', # separate elements with spaces; don't include stoichiometric subscripts, e.g. for Sr3As2, do 'Sr As'; if desired, some specific isotopes can be specified, e.g. instead of 'Fe Si', could do '56Fe 28Si'
    'space group':'Pnma', # don't worry about spaces
    'T and Ei conditions':[(295,17)], # (Temperature, Indicent energy) for each dataset in a list. MDEs must follow consistent naming scheme you specify in define_data.py
    'use seeKpath path':False, # need to fix, POSCAR required
    'user defined Qpoints':{'path':[ ('Y', 'Gamma'), ('Gamma','X'), ('X','Gamma2'), ('Gamma2','X'), ('X','Gamma'), ('Gamma','Z') ],#, ('K', 'Gamma'), ('Gamma', 'L'), ('L', 'W'), ('W', 'X') ], # e.g. [ ('X', 'U'), ('K', 'Gamma') ], total piece-wise path along which to process 2D inelastic slices, names must match strings in 'point_coords'
                            '1d_points':['Gamma','X'], # e.g. ['Gamma', 'X'], all 1d points at which to make 1d S(E) cuts, names must match strings in 'point_coords'
                            'point_coords':{ 'Gamma': [0.0, 0.0, 0.0], 'X': [0.5, 0.0, 0.0], 'Y': [0.0, 0.5, 0.0], 'Z': [0.0, 0.0, 0.5], 'Gamma2': [1, 0.0, 0.0] }}, # Q point definitions in Q basis associated with MDE's UB matrix
    'primitive to mantid transformation matrix': np.array([ [1,0,0],[0,1,0],[0,0,1] ]), # matrix defining transformation from primtive reciprocal lattice basis to Q basis associated with MDE's UB matrix
    
    # per-slice binning and integration ranges
    'E bins':'0,0.2,17', # min, step, max for energy binning
    'qdim0 step size':0.025, # if 2d 'path', this is step along qdim0 in rlu; if 1d 'point_coords', this is +/- integration range for qdim0 in rlu
    'qdim1 integration range':0.05, # +/- integration range for qdim1 in rlu
    'qdim2 integration range':0.05, # +/- integration range for qdim2 in rlu

    # which data to process and how
    'all symmetric 2d slices':True, # whether to process all symmetrically equivalent 2d inelastic slices (also includes bz folding)
    'all symmetric 1d cuts':False, # whether to process all symmetrically equivalent 1d cuts
    'BZ to process':0.2, # float or list of lists; if float between [0,1], then BZ with fractional qE coverage above provided threshold are processed (e.g. 0.5); if list of lists, data in each listed BZ will be processed (e.g. [[1,1,1], [2,0,0], [4,2,0]])
    'slice filter functions':[], # list of string(s), for automatic slice evaluation and filtering, string(s) must match filter function(s) specified in filter_functions.py
    
    # simulation details - POSCAR and FORCE_CONSTANTS required
    'perform simulations':True, # whether to perform analogous S(Q,E) simulations
    'supercell dimensions':[2,4,4], # supercell dimensions used to generate FORCE_CONSTANTS
    'use experimental coverage mask':True, # whether to impose experimental Q,E coverage mask on simulations
    'resolution blurring':(0.65,0.05), # (E FWHM in meV, Q FWHM in rlu, optional instrument string for SNS); if instrument, E FWHM should be approx. elastic line FWHM
    
    # saving and outputs
    'output directory':'../295K_17meV_foldWsims/', # absolute filepath where to save everything (will be created if doesn't exist)
    'save individual slices':True,
    'save reports':True,

    # misc.
    'u_vec':np.array([1,0,0]), # nominal - used for plots
    'v_vec':np.array([0,0,1]), # ^^
    'rhomb in hex notation':False # rare - if sample is rhombohedral but Mantid UB matrix uses hexagonal representation
    }

    return pathSQE_params

