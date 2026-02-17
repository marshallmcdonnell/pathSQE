import os
import numpy as np


########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(T_Ei_conditions, **kwargs):
    shared_folder = "/SNS/CNCS/IPTS-16687/shared/autoreduce/"
    raw_data_folder = "/SNS/CNCS/IPTS-16687/data/"
    mde_folder = "/SNS/CNCS/IPTS-16687/shared/Aiden/merged_mde/"

    data_set_list = []
    for T, E in T_Ei_conditions:
        data_set = {
            "Runs": list(
                range(180300, 180391 + 1)
            ),  # (15K,40meV):84040, 84159+1; (500K,40meV):79193, 79283+1; (750K,40meV):79092, 79192+1
            "BackgroundRuns": None,  # Options: None;list of runs that are added together
            "RawDataFolder": raw_data_folder,  # Options:raw_data_folder string
            "MdeFolder": mde_folder,  # Options:mde_folder string
            "MdeName": "SnS_{}K_{}meV_UBcorrected".format(T, E),
            "BackgroundMdeName": None,  # Options:None;bg_mde_name string
            "MaskingDataFile": shared_folder
            + "van180160.nxs",  # Options:None;data_file_name;
            "NormalizationDataFile": shared_folder
            + "van180160.nxs",  # Options:None;data_file_name
            "SampleLogVariables": {
                "OmegaMotorName": "CCR10G2Rot",
                "Temperature": T,
            },  # Options:None;LogVariableName;number
            # initial UB estimate  from jen is here:
            #'UBSetup':{'a':11.2,'b':3.987,'c':4.334,'alpha':90,'beta':90,'gamma':90,'u':'1,0,0','v':'0,0,1'},
            # refined UB matrix  -1.09721e-17,0.0131267,-0.230418,0.00467285,-0.250128,-0.0120591,-0.0891634,-0.0131087,-0.000631992
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
