#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from utils.utils import BMS_instance
import sys
sys.path.append("./")
from utils.utils import BMS_instance

nsteps = 5000
nsaves = 5000

# Haber_Bosch_experiments
Haber_BoschNS = BMS_instance(experiment = "Haber_Bosch", noutputs = 1, chosen_output = 1, scaling = None)
Haber_BoschNS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves)
Haber_BoschTS = BMS_instance(experiment = "Haber_Bosch", noutputs = 1, chosen_output = 1, scaling = "both")
Haber_BoschTS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves)

# Ammonia Oxidation experiments
for ammonia_output in range(2):
    Ammonia_NS  = BMS_instance(experiment = "Ammonia_Oxidation", noutputs = 2, chosen_output = ammonia_output + 1, scaling = None)
    Ammonia_NS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves)
    Ammonia_TS  = BMS_instance(experiment = "Ammonia_Oxidation", noutputs = 2, chosen_output = ammonia_output + 1, scaling = "both")
    Ammonia_TS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves)


# Propene Ammoxidation experiments
for propene_output in range(5):
    Propene_NS  = BMS_instance(experiment = "Propene_Ammoxidation", noutputs = 5, chosen_output = propene_output + 1, scaling = None)
    Propene_NS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves)
    Propene_TS  = BMS_instance(experiment = "Propene_Ammoxidation", noutputs = 5, chosen_output = propene_output + 1, scaling = "both")
    Propene_TS.run_BMS(mcmcsteps = nsteps, save_distance = nsaves)