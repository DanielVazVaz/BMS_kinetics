#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.append("./")
import os
from utils.utils import BMS_instance

nsteps = 5000
nsaves = 5000

LOCAL_FOLDER = r"Local"

# Haber_Bosch_experiments
Haber_BoschNS = BMS_instance(experiment = "Haber_Bosch", noutputs = 1, chosen_output = 1, scaling = None,
                             data_path = os.path.join(LOCAL_FOLDER, r"Data\Haber_Bosch\data.xlsx"))
Haber_BoschNS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves, save_folder_path=os.path.join(LOCAL_FOLDER,r"Results"))
Haber_BoschTS = BMS_instance(experiment = "Haber_Bosch", noutputs = 1, chosen_output = 1, scaling = "both",
                             data_path = os.path.join(LOCAL_FOLDER, r"Data\Haber_Bosch\data.xlsx"))
Haber_BoschTS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves, save_folder_path=os.path.join(LOCAL_FOLDER,r"Results"))

# Ammonia Oxidation experiments
for ammonia_output in range(2):
    Ammonia_NS  = BMS_instance(experiment = "Ammonia_Oxidation", noutputs = 2, chosen_output = ammonia_output + 1, scaling = None,
                               data_path = os.path.join(LOCAL_FOLDER, r"Data\Ammonia_Oxidation\data.xlsx"))
    Ammonia_NS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves, save_folder_path=os.path.join(LOCAL_FOLDER,r"Results"))
    Ammonia_TS  = BMS_instance(experiment = "Ammonia_Oxidation", noutputs = 2, chosen_output = ammonia_output + 1, scaling = "both",
                               data_path = os.path.join(LOCAL_FOLDER, r"Data\Ammonia_Oxidation\data.xlsx"))
    Ammonia_TS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves, save_folder_path=os.path.join(LOCAL_FOLDER,r"Results"))


# Propene Ammoxidation experiments
for propene_output in range(5):
    Propene_NS  = BMS_instance(experiment = "Propene_Ammoxidation", noutputs = 5, chosen_output = propene_output + 1, scaling = None,
                               data_path = os.path.join(LOCAL_FOLDER, r"Data\Propene_Ammoxidation\data.xlsx"))
    Propene_NS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves, save_folder_path=os.path.join(LOCAL_FOLDER,r"Results"))
    Propene_TS  = BMS_instance(experiment = "Propene_Ammoxidation", noutputs = 5, chosen_output = propene_output + 1, scaling = "both",
                               data_path = os.path.join(LOCAL_FOLDER, r"Data\Propene_Ammoxidation\data.xlsx"))
    Propene_TS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves, save_folder_path=os.path.join(LOCAL_FOLDER,r"Results"))