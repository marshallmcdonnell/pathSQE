import os
import numpy as np


########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(T_Ei_conditions, **kwargs):
    shared_folder = "/SNS/ARCS/IPTS-5307/shared/mantid_reduction/"
    raw_data_folder = "/SNS/ARCS/IPTS-5307/data/"
    mde_folder = "/SNS/ARCS/IPTS-21211/shared/Aiden/merged_mde/"

    data_set_list = []
    for T, E in T_Ei_conditions:
        data_set = {
            "Runs": None,
            "BackgroundRuns": None,  # range(297325,297337)Options: None;list of runs that are added together
            "RawDataFolder": raw_data_folder,  # Options:raw_data_folder string
            "MdeFolder": mde_folder,  # Options:mde_folder string
            "MdeName": "FeSi_{}K_{}meV_UBcorrected".format(T, E),
            "BackgroundMdeName": None,  # Options:None;bg_mde_name string
            "MaskingDataFile": shared_folder
            + "van_mask_Aiden.nxs",  # Options:None;data_file_name;
            "NormalizationDataFile": shared_folder
            + "van_mask_Aiden.nxs",  # Options:None;data_file_name
            "SampleLogVariables": {
                "OmegaMotorName": "ccr12rot",
                "Temperature": T,
            },  # Options:None;LogVariableName;number
            # did ub correction for each temp and resaved mde but initial UB estimate is here:
            #'UBSetup':{'a':4.449,'b':4.449,'c':4.449,'alpha':90,'beta':90,'gamma':90,'u':'-0.0562,-0.0797,1','v':'0.9809,1,0.1348'}, #by hand
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
