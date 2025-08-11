########################################################################################################
# Define dictionary input file for running pathSQE
# Author: Aiden Sable. August 2025.
########################################################################################################


import numpy as np

def define_pathSQE_params(**kwargs):
    pathSQE_params={
    # Sample and Q point info
    'sample':'Ge', # separate elements with spaces; don't include stoichiometric subscripts, e.g. for Sr3As2, do 'Sr As'; if desired, some specific isotopes can be specified, e.g. instead of 'Fe Si', could do '56Fe 28Si'
    'space group':'Fd-3m', # don't worry about spaces
    'T and Ei conditions':[(5,40)], # (Temperature, Indicent energy) for each dataset in a list. MDEs must follow consistent naming scheme you specify in define_data.py
    'use seeKpath path':False, # POSCAR required
    'user defined Qpoints': {
        'path': [
            ('0s', '0e'),
            ('1s', '1e'),
            ('2s', '2e'),
            ('3s', '3e'),
            ('4s', '4e'),
            ('5s', '5e'),
            ('6s', '6e'),
            ('7s', '7e'),
            ('8s', '8e'),
            ('9s', '9e'),
            ('10s', '10e'),
            ('11s', '11e'),
            ('12s', '12e'),
            ('13s', '13e'),
            ('14s', '14e')
        ],
        '1d_points': [],
        'point_coords': {
            '0s': [0.00, 0.00, 0.00], '0e': [1.00, 0.00, 0.00],
            '1s': [0.05, 0.05, 0.00], '1e': [1.00, 0.05, 0.00],
            '2s': [0.10, 0.10, 0.00], '2e': [1.00, 0.10, 0.00],
            '3s': [0.15, 0.15, 0.00], '3e': [1.00, 0.15, 0.00],
            '4s': [0.20, 0.20, 0.00], '4e': [1.00, 0.20, 0.00],
            '5s': [0.25, 0.25, 0.00], '5e': [1.00, 0.25, 0.00],
            '6s': [0.30, 0.30, 0.00], '6e': [1.00, 0.30, 0.00],
            '7s': [0.35, 0.35, 0.00], '7e': [1.00, 0.35, 0.00],
            '8s': [0.40, 0.40, 0.00], '8e': [1.00, 0.40, 0.00],
            '9s': [0.45, 0.45, 0.00], '9e': [1.00, 0.45, 0.00],
            '10s': [0.50, 0.50, 0.00], '10e': [1.00, 0.50, 0.00],
            '11s': [0.55, 0.55, 0.00], '11e': [0.95, 0.55, 0.00],
            '12s': [0.60, 0.60, 0.00], '12e': [0.90, 0.60, 0.00],
            '13s': [0.65, 0.65, 0.00], '13e': [0.85, 0.65, 0.00],
            '14s': [0.70, 0.70, 0.00], '14e': [0.80, 0.70, 0.00]
        }
    },
    'primitive to mantid transformation matrix': np.array([[-1,1,1],[1,-1,1],[1,1,-1]]), # matrix defining transformation from primtive reciprocal lattice basis to Q basis associated with MDE's UB matrix
    
    # per-slice binning and integration ranges
    'E bins':'0,0.5,40', # min, step, max for energy binning
    'qdim0 step size':0.025, # if 2d 'path', this is step along qdim0 in rlu; if 1d 'point_coords', this is +/- integration range for qdim0 in rlu
    'qdim1 integration range':0.025, # +/- integration range for qdim1 in rlu
    'qdim2 integration range':0.025, # +/- integration range for qdim2 in rlu

    # which data to process and how
    'all symmetric 2d slices':True, # whether to process all symmetrically equivalent 2d inelastic slices (also includes bz folding)
    'all symmetric 1d cuts':False, # whether to process all symmetrically equivalent 1d cuts
    'BZ to process':0.89, # float or list of lists; if float between [0,1], then BZ with fractional qE coverage above provided threshold are processed (e.g. 0.5); if list of lists, data in each listed BZ will be processed (e.g. [[1,1,1], [2,0,0], [4,2,0]])
    'slice filter functions':[], # list of string(s), for automatic slice evaluation and filtering, string(s) must match filter function(s) specified in filter_functions.py
    
    # simulation details - POSCAR and FORCE_CONSTANTS required
    'perform simulations':False, # whether to perform analogous S(Q,E) simulations
    'supercell dimensions':[4,4,4], # supercell dimensions used to generate FORCE_CONSTANTS
    'use experimental coverage mask':True, # whether to impose experimental Q,E coverage mask on simulations
    'resolution blurring':(0.8,0.05,'ARCS'), # (E FWHM in meV, Q FWHM in rlu, optional instrument string for SNS); if instrument, E FWHM should be approx. elastic line FWHM
    
    # saving and outputs
    'output directory':'../Ge_5K_40meV_fold_BZ_plane/', # where to save everything (will be created if doesn't exist)
    'save individual slices':False,
    'save reports':True,

    # misc.
    'u_vec':np.array([1,0,0]), # nominal - used for plots
    'v_vec':np.array([0,1,0]), # ^^
    'rhomb in hex notation':False # rare - if sample is rhombohedral but Mantid UB matrix uses hexagonal representation
    }

    return pathSQE_params




