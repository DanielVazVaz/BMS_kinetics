#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
BMS_FOLDER = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS"
sys.path.append(BMS_FOLDER)
sys.path.append(BMS_FOLDER + "\\" + "Prior")
import pandas as pd
import numpy as np
from copy import deepcopy
import mcmc
import parallel
from fit_prior import read_prior_par
import warnings
warnings.filterwarnings('ignore')
from tqdm.auto import tqdm

class PriorError(Exception):
    pass
class ScalingError(Exception):
    pass

class BMS_instance:
    def __init__(self, experiment: str, noutputs:int = 1, chosen_output:int = 1):
        """
        Load the data of an experiment
        """
        folder_data = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Data"
        file_data = folder_data + "\\" + experiment + "\\data.xlsx"
        self.train = pd.read_excel(file_data, sheet_name="train", index_col = 0)
        self.test  = pd.read_excel(file_data, sheet_name="test", index_col = 0)
        self.noutputs = noutputs
        self.ninputs  = self.train.shape[1] - self.noutputs
        self.init_prior()
        self.init_tree(chosen_output=1)
        
    def init_prior(self):
        prior_folder = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS\Prior"
        prior_files  = os.listdir(prior_folder)
        self.valid_priors = [i for i in prior_files if ".nv{0}.".format(self.ninputs) in i]
        self.chosen_prior = self.valid_priors[-1]
        if not self.chosen_prior:
            raise PriorError
        self.npar = int(re.search(r"np(\d+)", self.chosen_prior).group(1))
        self.prior_par = read_prior_par(prior_folder + "\\" + self.chosen_prior)
    
    def init_tree(self, chosen_output = 1, scaling = None):
        self.x  = self.train.iloc[:,:self.ninputs].copy()
        self.y  = self.train.iloc[:, self.ninputs + chosen_output - 1].copy()
        self.y  = pd.Series(list(self.y))
        self.x.columns = ["x" + str(i+1) for i in range(len(self.x.columns))]
        if scaling == "inputs":
            self.x = (self.x - self.x.mean())/self.x.std()
        elif scaling == "outputs":
            self.y = (self.y - self.y.mean())/self.y.std()
        elif scaling == "both":
            self.x = (self.x - self.x.mean())/self.x.std()
            self.y = (self.y - self.y.mean())/self.y.std()
        elif not scaling:
            pass
        else:
            raise ScalingError
        
        Ts = [1] + [1.04**k for k in range(1, 20)]
        self.pms = parallel.Parallel(
            Ts,
            variables= self.x.columns.tolist(),
            parameters=['a%d' % i for i in range(self.npar)],
            x=self.x, y=self.y,
            prior_par=self.prior_par,
        )
        
    def run_BMS(self, mcmcsteps = 232):
        self.description_lengths, self.mdl, self.mdl_model = [], np.inf, None
        pbar = tqdm(range(mcmcsteps), desc = "Running BMS: ")
        for i in pbar:
            pbar.set_description("Running BMS: [{0}/{1}]: ".format(i+1, mcmcsteps))
            # MCMC update
            self.pms.mcmc_step() # MCMC step within each T
            self.pms.tree_swap() # Attempt to swap two randomly selected consecutive temps
            # Add the description length to the trace
            self.description_lengths.append(self.pms.t1.E)
            # Check if this is the MDL expression so far
            # Thinning here
            if self.pms.t1.E < self.mdl:
                self.mdl, self.mdl_model = self.pms.t1.E, deepcopy(self.pms.t1)
            # Update the progress bar
            #f.value += 1
            #f.description = 'Run:{0}'.format(i)
            
if __name__ == "__main__":
    a = BMS_instance("Haber_Bosch")     
    a.run_BMS()
    print(a.mdl_model)