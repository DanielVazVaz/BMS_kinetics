#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
#BMS_FOLDER = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS"
#sys.path.append(BMS_FOLDER)
#sys.path.append(BMS_FOLDER + "\\" + "Prior")
import pandas as pd
import numpy as np
from copy import deepcopy
import parallel
from fit_prior import read_prior_par
import warnings
warnings.filterwarnings('ignore')
from tqdm.auto import tqdm
import pickle
import matplotlib.pyplot as plt

OPS = {
    'sin': 1,
    'cos': 1,
    'tan': 1,
    'exp': 1,
    'log': 1,
    'sinh' : 1,
    'cosh' : 1,
    'tanh' : 1,
    'pow2' : 1,
    'pow3' : 1,
    'abs'  : 1,
    'sqrt' : 1,
    'fac' : 1,
    '-' : 1,
    '+' : 2,
    '*' : 2,
    '/' : 2,
    '**' : 2,
}

class PriorError(Exception):
    pass
class ScalingError(Exception):
    pass

class BMS_instance:
    def __init__(self, experiment: str = None, noutputs:int = 1, chosen_output:int = 1, scaling:str = None, data_path:str = r"data.xslx", ops:dict = OPS,
                 prior_folder:str = r".", npar:int = None):
        """Initialize an instance of hte BMS class

        Args:
            experiment (str, optional): Name of the experiment if needed. Defaults to None.
            noutputs (int, optional): Number of outputs in the data. Defaults to 1.
            chosen_output (int, optional): Position of the chosen output, from 1 to noutputs. Defaults to 1.
            scaling (str, optional): Perform zscore scaling for inputs, outputs, both, or None. Defaults to None.
            data_path (str, optional): Path to the xlsx file. Defaults to r"data.xslx".
            ops (dict, optional): Dictionary with the valid operations and the children of each operation. Defaults to OPS.
            prior_folder (str, optional): Path to the prior folder. Defaults to r".".
            npar (int, optional): Number of parameters to consider. Defaults to None.
        """        
        self.experiment = experiment
        file_data = data_path
        self.prior_folder = prior_folder
        self.ops = ops
        self.train = pd.read_excel(file_data, sheet_name="train", index_col = 0)
        self.test  = pd.read_excel(file_data, sheet_name="test", index_col = 0)
        self.noutputs = noutputs
        self.ninputs  = self.train.shape[1] - self.noutputs
        self.init_prior(npar)
        self.init_tree(scaling = scaling, chosen_output=chosen_output)
        self.chosen_output = chosen_output
        
    @staticmethod
    def load(load):
        """Static method to load pickled saved models.

        Args:
            load (str): Raw string that indicates the full path to the pickle (.pkl) file.

        Returns:
            BMS_instance: Loaded BMS instance
        """        
        with open(load, "rb") as input_file:
                return pickle.load(input_file)
        
    def init_prior(self, npar:int):
        """Initialize the prior considered for the BMS. It chooses the last element of a valid list of priors regarding
        the dimensionality of the input if no integer is passed.

        Args:
            npar (int): Numnber of parameters to consider

        Raises:
            PriorError: In case there is no priors available
        """    
        prior_folder = self.prior_folder
        prior_files  = os.listdir(prior_folder)
        self.valid_priors = [i for i in prior_files if ".nv{0}.".format(self.ninputs) in i]
        if npar:
            self.chosen_prior = [i for i in self.valid_priors if ".np{0}.".format(npar) in i][-1]
        else:
            self.chosen_prior = self.valid_priors[-1]
        if not self.chosen_prior:
            raise PriorError
        self.npar = int(re.search(r"np(\d+)", self.chosen_prior).group(1))
        self.prior_par = read_prior_par(prior_folder + "\\" + self.chosen_prior)
    
    def init_tree(self, chosen_output:int = 1, scaling:str = None):
        """Initializes the parallel BMS tree for the chosen output. Also applies z-score scaling if the option is selected.

        Args:
            chosen_output (int, optional): Position of the chosen output, from 1 to noutputs. Defaults to 1.
            scaling (str, optional): Perform zscore scaling for inputs, outputs, both, or None. Defaults to None.

        Raises:
            ScalingError: Wrong input to scaling. Options are 'inputs', 'outputs', 'both', and None.
        """        
        self.x_test  = self.test.iloc[:,:self.ninputs].copy()
        self.x_test.reset_index(inplace=True, drop = True)
        self.y_test  = self.test.iloc[:, self.ninputs + chosen_output - 1].copy()
        self.y_test  = pd.Series(list(self.y_test))
        self.x_test.columns = ["x" + str(i+1) for i in range(len(self.x_test.columns))]
        self.x  = self.train.iloc[:,:self.ninputs].copy()
        self.x.reset_index(inplace=True, drop = True)
        self.y  = self.train.iloc[:, self.ninputs + chosen_output - 1].copy()
        self.y  = pd.Series(list(self.y))
        self.original_columns = self.x.columns
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
            ops = self.ops
        )
        
    def run_BMS(self, mcmcsteps:int = 232, save_distance:int = None, save_folder_path:str = r".\saved_model.pkl"):
        """Runs the BMS for a number of mcmcsteps. Also saves the resultant model in a .pkl each save_distance points.

        Args:
            mcmcsteps (int, optional): Number of mcmcsteps to run. Defaults to 232.
            save_distance (int, optional): Number of mcmcsteps before saving a .pkl file. Defaults to None.
            save_folder_path (str, optional): Path to the saved .pkl file. Defaults to r".\saved_model.pkl".
        """
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
            # Save pickle
            if save_distance:
                if (i+1)%save_distance == 0:
                    with open(save_folder_path + "\\" + r'{2}_Out{3}_Scale{0}_{1}.pkl'.format(self.scaling, (i+1), self.experiment, self.chosen_output), 'wb') as outp:
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
        
    def calculate_rsquare(self, plot:bool = False) -> dict:
        """Calculate the predicted vs. real graph, as well as information about the best found model.

        Args:
            plot (bool, optional): Indicates if the plot of predicted vs calculated will be displayed. Defaults to False.

        Returns:
            dict: Information about the error metrics
        """        
        if plot:
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
        if plot:
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
        error_dict = {"R2_train": R2_train, "R2_test": R2_test, "MSE_train": MSE_train,
                      "MSE_test": MSE_test}
        return error_dict
    
    def save_txt(self, file_path: str):
        """Saves the BMS equation in a .txt file, as well as the value of the parameters

        Args:
            file_path (str): Path to the txt file
        """
        with open(file_path, "w") as f:
            f.writelines([str(self.mdl_model), str(self.mdl_model.par_values["d0"])] )
            

if __name__ == "__main__":
    a = BMS_instance("Haber_Bosch")     
    a.run_BMS(300)
    print(a.mdl_model)