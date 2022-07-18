#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from utils.utils import BMS_instance
import sys
sys.path.append("./")
from utils.utils import BMS_instance

# Haber_Bosch_experiments
#Haber_BoschNS = BMS_instance(experiment = "Haber_Bosch", noutputs = 1, chosen_output = 1, scaling = None)
#Haber_BoschNS.run_BMS(mcmcsteps = 5000, save_distance = 500)
Haber_BoschTS = BMS_instance(experiment = "Haber_Bosch", noutputs = 1, chosen_output = 1, scaling = "both")
Haber_BoschTS.run_BMS(mcmcsteps = 4500, save_distance = 500)