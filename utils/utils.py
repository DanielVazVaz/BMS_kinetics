#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
from .config import BMS_FOLDER
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
import pickle
import matplotlib.pyplot as plt

class PriorError(Exception):
    pass
class ScalingError(Exception):
    pass

class BMS_instance:
    def __init__(self, experiment: str = None, noutputs:int = 1, chosen_output:int = 1, scaling = None):
        """
        Initialize the instance. 
        
        Inputs:
            - experiment    : String indicating the name of the folder where the data.xlsx file is.
            - noutputs      : Integer indicating the number of outputs in this experiment.
            - chosen_output : Integer indicating the chosen output for this BMS instance.
            - scaling       : String with the options "inputs", "outputs", and "both", which indicates to which
                              section of the data a z-score scaling should be performed. 
        
        Called methods:
            - init_prior()
            - init_tree(scaling, chosen_output)

        """
        self.experiment = experiment
        folder_data = r"Local\Data"
        file_data = folder_data + "\\" + experiment + "\\data.xlsx"
        self.train = pd.read_excel(file_data, sheet_name="train", index_col = 0)
        self.test  = pd.read_excel(file_data, sheet_name="test", index_col = 0)
        self.noutputs = noutputs
        self.ninputs  = self.train.shape[1] - self.noutputs
        self.init_prior()
        self.init_tree(scaling = scaling, chosen_output=chosen_output)
        self.chosen_output = chosen_output
        
    @staticmethod
    def load(load):
        """
        Static method to load pickled saved models. 
        
        Inputs:
            - load: Raw string that indicates the full path to the pickle (.pkl) file.
        """
        with open(load, "rb") as input_file:
                return pickle.load(input_file)
        
    def init_prior(self):
        """
        Initialize the prior considered for the BMS. It chooses the last element of a valid list of priors regarding
        the dimensionality of the input.
        """
        prior_folder = r"Local\BMS\Prior"
        prior_files  = os.listdir(prior_folder)
        self.valid_priors = [i for i in prior_files if ".nv{0}.".format(self.ninputs) in i]
        self.chosen_prior = self.valid_priors[-1]
        if not self.chosen_prior:
            raise PriorError
        self.npar = int(re.search(r"np(\d+)", self.chosen_prior).group(1))
        self.prior_par = read_prior_par(prior_folder + "\\" + self.chosen_prior)
    
    def init_tree(self, chosen_output = 1, scaling = None):
        """
        Initializes the parallel BMS tree for the chosen output. Also applies z-score scaling if the option is selected.
        
        Inputs:
            - chosen_output : Integer indicating the chosen output for this BMS instance.
            - scaling       : String with the options "inputs", "outputs", and "both", which indicates to which
                              section of the data a z-score scaling should be performed.
        """
        self.x_test  = self.test.iloc[:,:self.ninputs].copy()
        self.y_test  = self.test.iloc[:, self.ninputs + chosen_output - 1].copy()
        self.y_test  = pd.Series(list(self.y_test))
        self.x_test.columns = ["x" + str(i+1) for i in range(len(self.x_test.columns))]
        self.x  = self.train.iloc[:,:self.ninputs].copy()
        self.y  = self.train.iloc[:, self.ninputs + chosen_output - 1].copy()
        self.y  = pd.Series(list(self.y))
        self.x.columns = ["x" + str(i+1) for i in range(len(self.x.columns))]
        self.scaling = scaling
        x_scaling_mean = self.x.mean()
        y_scaling_mean = self.y.mean()
        x_scaling_std  = self.x.std()
        y_scaling_std  = self.y.std()
        self.scaling_param = {"mean": [x_scaling_mean, y_scaling_mean], "std": [x_scaling_std, y_scaling_std]}
        if scaling == "inputs":
            self.x = (self.x - x_scaling_mean)/x_scaling_std
            self.x_test = (self.x_test - x_scaling_mean)/x_scaling_std
        elif scaling == "outputs":
            self.y = (self.y - y_scaling_mean)/y_scaling_std
            self.y_test = (self.y_test - y_scaling_mean)/y_scaling_std
        elif scaling == "both":
            self.x = (self.x - x_scaling_mean)/x_scaling_std
            self.x_test = (self.x_test - x_scaling_mean)/x_scaling_std
            self.y = (self.y - y_scaling_mean)/y_scaling_std
            self.y_test = (self.y_test - y_scaling_mean)/y_scaling_std
        elif not scaling:
            self.scaling = "None"
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
        
    def run_BMS(self, mcmcsteps = 232, save_distance = 100):
        """
        Runs the BMS for a number of mcmcsteps. Also saves the resultant model in a .pkl each save_distance points.
        
        Inputs:
            - mcmcsteps     : Integer that indicates the number of markov chain monte carlo steps to perform.
            - save_distance : Integer that indicates how many steps occur between saving an instance of the model to a pkl. file.
        """
        self.description_lengths, self.mdl, self.mdl_model = [], np.inf, None
        pbar = tqdm(range(mcmcsteps), desc = "Running BMS: ")
        result_folder = r"Local\Results"
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
            # Save pickle
            if (i+1)%save_distance == 0:
                with open(result_folder + "\\" + r'{2}_Out{3}_Scale{0}_{1}.pkl'.format(self.scaling, (i+1), self.experiment, self.chosen_output), 'wb') as outp:
                    pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)
                    
    def plot_dlength(self):
        """
        Plot the description length graph, as well as information about the best found model.
        """
        ## WE check some stats:
        print('Best model:\t', self.mdl_model)
        print('Desc. length:\t', self.mdl)
        print("BIC:\t", self.mdl_model.bic)
        print("Par:\t", self.mdl_model.par_values["d0"])
        ## We display the DL
        plt.figure(figsize=(15, 5))
        plt.plot(self.description_lengths)
        
        plt.xlabel('MCMC step', fontsize=14)
        plt.ylabel('Description length', fontsize=14)
        
        plt.title('MDL model: $%s$' % self.mdl_model.latex())
        plt.show()
        
    def plot_rsquare(self):
        """
        Plot the predicted vs. real graph, as well as information about the best found model.
        """
        print('Best model:\t', self.mdl_model)
        print('Desc. length:\t', self.mdl)
        print("BIC:\t", self.mdl_model.bic)
        print("Par:\t", self.mdl_model.par_values["d0"])
        
        predicted_train = self.mdl_model.predict(self.x)
        predicted_test  = self.mdl_model.predict(self.x_test)
        SSR_train = ((self.y - predicted_train)**2).sum()
        SST_train = ((self.y - self.y.mean())**2).sum()
        R2_train  = 1 - SSR_train/SST_train 
        SSR_test  = ((self.y_test - predicted_test)**2).sum()
        SST_test  = ((self.y_test - self.y_test.mean())**2).sum()
        R2_test   = 1 - SSR_test/SST_test
        MSE_train = ((self.y - predicted_train)**2).sum()/self.y.shape[0]
        MSE_test  = ((self.y_test - predicted_test)**2).sum()/self.y_test.shape[0]
        plt.figure(figsize=(6, 4), dpi = 200)
        ax = plt.gca()
        plt.scatter(predicted_train, self.y, color = "blue", label = "training")
        plt.scatter(predicted_test, self.y_test, color = "red", label = "test")
        lims = [np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
                np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
                ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        plt.xlabel('MDL model predictions', fontsize=14)
        plt.ylabel('Actual values', fontsize=14)
        plt.title(f"$R^2_{{train}} = {R2_train:.4f}\qquad R^2_{{test}} = {R2_test:.4f}$\n$MSE_{{train}} = {MSE_train:.3e} \qquad  MSE_{{test}} = {MSE_test:.3e}$")
        plt.legend(loc = "best")
        plt.show()
        
            
            
if __name__ == "__main__":
    a = BMS_instance("Haber_Bosch")     
    a.run_BMS(300)
    print(a.mdl_model)