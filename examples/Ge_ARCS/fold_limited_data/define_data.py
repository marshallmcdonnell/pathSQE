import os
import numpy as np


########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(T_Ei_conditions, **kwargs):
    shared_folder = "/SNS/ARCS/IPTS-13861/shared/jen/reduction/"
    raw_data_folder = "/SNS/ARCS/IPTS-13861/data/"
    mde_folder = "/SNS/ARCS/IPTS-13861/shared/Aiden/comprehensive/individAngle_MDEs/"

    data_set_list = []
    for T, E in T_Ei_conditions:
        data_set = {
            "Runs": 69945,  # going to sub this in using driver
            "BackgroundRuns": None,
            "RawDataFolder": raw_data_folder,
            "MdeFolder": mde_folder,
            "MdeName": "Ge_run69945",  # also subbing in using driver
            "BackgroundMdeName": "",
            "MaskingDataFile": shared_folder + "van67625.nxs",
            "NormalizationDataFile": shared_folder + "van67625.nxs",
            "SampleLogVariables": {
                "OmegaMotorName": "CCR16Rot",
                "Temperature": T,
            },  # hdfview on event.nxs, under entry->DASlogs
            "UBSetup": {
                "a": 5.66,
                "b": 5.66,
                "c": 5.66,
                "alpha": 90,
                "beta": 90,
                "gamma": 90,
                "u": "1,0,0",
                "v": "0,1,0",
            },
            # Data reduction options
            "Ei": None,
            "T0": None,
            "BadPulsesThreshold": None,
            "TimeIndepBackgroundWindow": None,
            "E_min": None,
            "E_max": None,
            "AdditionalDimensions": None,
            "EfCorrectionFunction": None,
        }
    data_set_list.append(data_set)

    return data_set_list
