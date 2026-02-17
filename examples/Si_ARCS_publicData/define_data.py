import os
import numpy as np


########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(T_Ei_conditions, **kwargs):
    shared_folder = "/SNS/ARCS/2013_2_18_CAL/shared/silicon/"
    raw_data_folder = "/SNS/ARCS/2013_2_18_CAL/data/"
    mde_folder = "~"  # where to save or search for mde

    data_set_list = []
    for T, E in T_Ei_conditions:
        data_set = {
            "Runs": list(
                range(37769, 37859 + 1)
            ),  # 5K: list(range(43977, 44083+1));  300K: list(range(37769, 37859+1))
            "BackgroundRuns": None,  # Options: None;list of runs that are added together
            "RawDataFolder": raw_data_folder,  # Options:raw_data_folder string
            "MdeFolder": mde_folder,  # Options:mde_folder string
            "MdeName": "Si_{}K_{}meV_UBcorrected_ub".format(T, E),
            "BackgroundMdeName": "",  # Options:None;bg_mde_name string
            "MaskingDataFile": shared_folder
            + "van37349_white_upstatic.nxs",  # Options:None;data_file_name
            "NormalizationDataFile": shared_folder
            + "van37349_white_upstatic.nxs",  # Options:None;data_file_name
            "SampleLogVariables": {
                "OmegaMotorName": "CCR12Rot",
                "Temperature": T,
            },  # 5K: 'CCR16Rot'; 300K: 'CCR12Rot'
            #'UBSetup':{'a':5.45,'b':5.45,'c':5.45,'alpha':90,'beta':90,'gamma':90, 'u':'-2.87720207, -3.57259257,  -2.94292212', 'v':'4.04505632, -3.62506976, 0.44596929'}, # For Si 5K:
            "UBSetup": {
                "ub": "0.12685824,-0.13193672,-0.01291855,0.06613873,0.07848529,-0.15209513,-0.11489062,-0.10049857,-0.10182034"
            },  # For Si 300K:
            # Data reduction options
            "Ei": None,  # Options: None;Ei_somehow_determined
            "T0": None,  # Options: None;T0_determined_from_mantid
            "BadPulsesThreshold": None,  # Options: None;bg_pulses_threshold value
            "TimeIndepBackgroundWindow": None,  # Options: None;[Tib_min,Tib_max]
            "E_min": None,  # Options: None;Emin
            "E_max": None,  # Options: None;Emax
            "AdditionalDimensions": None,  # Options: None;list of triplets ("name", min, max)
            "EfCorrectionFunction": None,  # Options:None;'HYSPEC_default_correction';Custom_Ef_Correction_Function_Name
        }
        data_set_list.append(data_set)

    return data_set_list
