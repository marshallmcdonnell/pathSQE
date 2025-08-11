import os
import numpy as np

########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace 
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(T_Ei_conditions, **kwargs):
    shared_folder='/SNS/ARCS/IPTS-13861/shared/jen/reduction/'
    raw_data_folder='/SNS/ARCS/IPTS-13861/data/'
    mde_folder='~'

    data_set_list=[]
    for T, E in T_Ei_conditions:
        data_set={'Runs':list(range(69587, 69788+1)), #not accurate
              'BackgroundRuns':None,   # Options: None;list of runs that are added together
              'RawDataFolder':raw_data_folder,      # Options:raw_data_folder string
              'MdeFolder':mde_folder,               # Options:mde_folder string
              'MdeName':'merged_mde_MnO_{}meV_{}K_unpol'.format(E, T),   # Options:mde_name string
              'BackgroundMdeName':'',      # Options:None;bg_mde_name string
              'MaskingDataFile':None,         # Options:None;data_file_name; 
              'NormalizationDataFile':None,   # Options:None;data_file_name
              'SampleLogVariables':{'OmegaMotorName':'CCR16Rot','Temperature':T},   # found name by doing hdfview on *_event.nxs and checking under entry->DASlogs
              'UBSetup':{'a':5.66,'b':5.66,'c':5.66,'alpha':90,'beta':90,'gamma':90,'u':'1,0,0','v':'0,1,0'}, #
               #Data reduction options
              'Ei':None,                             # Options: None;Ei_somehow_determined
              'T0': None,                           # Options: None;T0_determined_from_mantid
              'BadPulsesThreshold':None,            # Options: None;bg_pulses_threshold value
              'TimeIndepBackgroundWindow':None,     # Options: None;[Tib_min,Tib_max]
              'E_min':None,                         # Options: None;Emin
              'E_max':None,                         # Options: None;Emax
              'AdditionalDimensions':None,          # Options: None;list of triplets ("name", min, max)
              'EfCorrectionFunction':None,          # Options:None;'HYSPEC_default_correction';Custom_Ef_Correction_Function_Name
              }
    data_set_list.append(data_set)


    return data_set_list

